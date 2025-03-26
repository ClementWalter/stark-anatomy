import math
from collections import namedtuple
from hashlib import blake2b
from typing import List, cast

from stark_anatomy.algebra import FieldElement
from stark_anatomy.ip import ProofStream
from stark_anatomy.merkle import Merkle
from stark_anatomy.univariate import Polynomial

# Type for colinearity check commitments
Codewords = namedtuple("Codewords", ["a", "b", "c"])
MerkleProofs = namedtuple("MerkleProof", ["a", "b", "c"])
ColinearityCheck = namedtuple("ColinearityCheck", ["codewords", "proofs"])


class FriError(Exception):
    pass


class Fri:
    def __init__(
        self,
        offset: FieldElement,
        omega: FieldElement,
        domain_length: int,
        expansion_factor: int,
        num_colinearity_tests: int,
    ):
        self.offset = offset
        self.omega = omega
        self.domain_length = domain_length
        self.field = omega.field
        self.expansion_factor = expansion_factor
        self.num_colinearity_tests = num_colinearity_tests

        if self.num_rounds < 1:
            raise ValueError("cannot do FRI with less than one round")

        if self.omega**self.domain_length != 1:
            raise ValueError("omega does not have the right order")

        if self.domain_length <= self.expansion_factor:
            raise ValueError("domain length must be at least expansion factor")

        if self.domain_length <= 4 * self.num_colinearity_tests:
            raise ValueError("domain length must be at least 4 * num_colinearity_tests")

    @property
    def num_rounds(self):
        """
        Compute the number of folding rounds of FRI.

        Given the fact that each round halves the domain length, we have
        len(initial_domain) = 2**num_rounds * len(final_domain).
        """
        return math.ceil(
            math.log2(
                self.domain_length
                / max(self.expansion_factor, 4 * self.num_colinearity_tests)
            )
        )

    @staticmethod
    def sample_index(seed: bytes, size: int) -> int:
        acc = 0
        for b in seed:
            acc = (acc << 8) ^ b
        return acc % size

    @staticmethod
    def sample_indices(
        seed: bytes, size: int, reduced_size: int, number: int
    ) -> list[int]:
        """
        Sample `number` distinct indices from `size`-sized domain,
        where the indices are also distinct when reduced modulo `reduced_size`.
        """
        if number > reduced_size:
            raise ValueError(
                f"cannot sample more indices than available in last codeword; requested: {number}, available: {reduced_size}"
            )
        if number > 2 * reduced_size:
            raise ValueError("not enough entropy in indices wrt last codeword")

        indices = []
        reduced_indices = []
        counter = 0
        while len(indices) < number:
            index = Fri.sample_index(blake2b(seed + bytes(counter)).digest(), size)
            reduced_index = index % reduced_size
            counter += 1
            if reduced_index not in reduced_indices:
                indices += [index]
                reduced_indices += [reduced_index]

        return indices

    @property
    def domain(self) -> List[FieldElement]:
        """
        The domain of the FRI.

        This is the set of points at which the polynomial is evaluated.
        """
        return [self.offset * self.omega**i for i in range(self.domain_length)]

    def commit(self, codeword: List[FieldElement], proof_stream: ProofStream):
        """
        Commit the codeword to the proof stream.

        Folding is done by splitting the codeword into odd and even parts of the polynomial.

        f*(x) = 1/2(f(x) + f(-x)) + 1/2 * alpha * (f(x) - f(-x)) / x
              = 1/2( (1 + alpha / x) * f(x) + (1 - alpha / x) * f(-x) )

        Due to the fact that -1 = omega**n/2, negative values are simply values with offset
        n/2 in the domain:

        f*(omega**(2*i)) = 1/2( (1 + alpha / omega**i) * f(omega**i) + (1 - alpha / omega**i) * f(-omega**i) )
                         = 1/2( (1 + alpha / omega**i) * f(omega**i) + (1 + alpha / omega**(n/2 + i) * f(omega**(n/2 + i)) )
        """
        if len(self.domain) != len(codeword):
            raise ValueError("initial codeword length does not match domain length")

        domain = self.domain
        codewords = []

        # for each round
        for _ in range(self.num_rounds):

            # add to list of codewords
            codewords += [codeword]

            # commit the codeword
            proof_stream.push(Merkle.commit(codeword))

            # get challenge
            alpha = self.field.sample(proof_stream.prover_fiat_shamir())

            # compute weights
            weights = [
                (alpha / omega + 1) * word for omega, word in zip(domain, codeword)
            ]

            # split and fold
            domain = domain[::2]
            n = len(domain)
            codeword = [
                (positive + negative) / 2
                for positive, negative in zip(weights[:n], weights[n:])
            ]

        # commit the last codeword
        proof_stream.push(codewords[-1])

        return codewords

    def query(
        self,
        current_codeword: List[FieldElement],
        next_codeword: List[FieldElement],
        c_indices: List[int],
        proof_stream: ProofStream,
    ):
        """
        Commit whatever is required to allow the verifier to perform the colinearity checks
        for the current codeword and the next codeword.

        A = (omega**i, f(omega**i))
        B = (-omega**i, f(-omega**i)) = (omega**(n/2 + i), f(omega**(n/2 + i)))
        C = (alpha, f*(omega**i*2))

        A and C are at the same index (both are omega**i) while B is shifted by n/2
        because -1 = omega**n/2.
        """
        for i in c_indices:
            codewords = Codewords(
                a=current_codeword[i],
                b=current_codeword[i + len(current_codeword) // 2],
                c=next_codeword[i],
            )
            proofs = MerkleProofs(
                a=Merkle.open(i, current_codeword),
                b=Merkle.open(i + len(current_codeword) // 2, current_codeword),
                c=Merkle.open(i, next_codeword),
            )
            proof_stream.push(ColinearityCheck(codewords, proofs))

    def prove(self, codeword: List[FieldElement], proof_stream: ProofStream):
        """
        Commit the codeword to the proof stream and add colinearity checks.
        """
        # commit phase
        codewords = self.commit(codeword, proof_stream)

        # get indices
        top_level_indices = self.sample_indices(
            proof_stream.prover_fiat_shamir(),
            len(codewords[0]) // 2,
            len(codewords[-1]),
            self.num_colinearity_tests,
        )

        # query phase
        indices = top_level_indices
        for current_codeword, next_codeword in zip(codewords[:-1], codewords[1:]):
            self.query(current_codeword, next_codeword, indices, proof_stream)
            indices = [index % (len(next_codeword) // 2) for index in indices]

        return top_level_indices

    def verify(self, proof_stream):
        """
        Verify the proof, i.e. mainly unwind the prover's computation in `prove()`, pulling from
        instead of pushing to the proof stream.
        """

        # Unwind `commit()`
        roots = []
        alphas = []
        for _ in range(self.num_rounds):
            roots += [proof_stream.pull()]
            alphas += [self.field.sample(proof_stream.verifier_fiat_shamir())]
        last_codeword = proof_stream.pull()

        # check if it matches the given root
        if roots[-1] != Merkle.commit(last_codeword):
            raise FriError("last codeword is not well formed")

        last_omega = self.omega ** (2 ** (self.num_rounds - 1))
        last_offset = self.offset ** (2 ** (self.num_rounds - 1))

        # assert that last_omega has the right order
        if last_omega ** len(last_codeword) != 1:
            raise FriError("omega does not have right order")

        # compute interpolant
        last_domain = [last_offset * (last_omega**i) for i in range(len(last_codeword))]
        poly = Polynomial.interpolate(last_domain, last_codeword)

        # check if it is low degree
        degree = len(last_codeword) // self.expansion_factor - 1
        if poly.degree > degree:
            raise FriError(
                "last codeword does not correspond to polynomial of low enough degree"
            )

        # get indices
        top_level_indices = self.sample_indices(
            proof_stream.verifier_fiat_shamir(),
            self.domain_length >> 1,
            self.domain_length >> (self.num_rounds - 1),
            self.num_colinearity_tests,
        )

        # Unwind `query()` loop
        omega = self.omega
        offset = self.offset
        domain_length = self.domain_length
        polynomial_values = []
        # for every round, check consistency of subsequent layers
        for r in range(0, self.num_rounds - 1):
            # fold c indices
            c_indices = [index % (domain_length // 2) for index in top_level_indices]
            a_indices = c_indices
            b_indices = [index + domain_length // 2 for index in c_indices]

            # read values and check colinearity
            for s in range(self.num_colinearity_tests):
                colinearity_check = cast(ColinearityCheck, proof_stream.pull())
                # colinearity check
                ax = offset * (omega ** a_indices[s])
                bx = offset * (omega ** b_indices[s])
                cx = alphas[r]
                ay, by, cy = colinearity_check.codewords

                if not Polynomial.is_colinear([(ax, ay), (bx, by), (cx, cy)]):
                    raise FriError("colinearity check failure")

                # verify authentication paths
                proofs = colinearity_check.proofs
                if not Merkle.verify(roots[r], a_indices[s], proofs.a, ay):
                    raise FriError(
                        "merkle authentication path verification fails for a"
                    )
                if not Merkle.verify(roots[r], b_indices[s], proofs.b, by):
                    raise FriError(
                        "merkle authentication path verification fails for b"
                    )
                if not Merkle.verify(roots[r + 1], c_indices[s], proofs.c, cy):
                    raise FriError(
                        "merkle authentication path verification fails for c"
                    )

                # record top-layer values for later verification
                if r == 0:
                    polynomial_values += [(a_indices[s], ay), (b_indices[s], by)]

            omega = omega**2
            offset = offset**2
            domain_length //= 2

        return polynomial_values
