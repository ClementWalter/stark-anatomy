import math

import pytest

from stark_anatomy.algebra import Field, FieldElement
from stark_anatomy.fri import Fri, FriError
from stark_anatomy.ip import ProofStream
from stark_anatomy.univariate import Polynomial


class TestFri:
    def test_fri(self, field: Field):
        degree = 63
        expansion_factor = 4
        num_colinearity_tests = 17

        initial_codeword_length = (degree + 1) * expansion_factor
        log_codeword_length = math.ceil(math.log2(initial_codeword_length))
        assert 1 << log_codeword_length == initial_codeword_length

        omega = field.primitive_nth_root(initial_codeword_length)
        offset = field.generator()

        fri = Fri(
            offset,
            omega,
            initial_codeword_length,
            expansion_factor,
            num_colinearity_tests,
        )

        polynomial = Polynomial([FieldElement(i, field) for i in range(degree + 1)])
        domain = [omega**i for i in range(initial_codeword_length)]

        codeword = polynomial.evaluate_domain(domain)

        # test valid codeword
        proof_stream = ProofStream()
        fri.prove(codeword, proof_stream)
        points = fri.verify(proof_stream)
        assert all(polynomial.evaluate(omega**x) == y for x, y in points)

        # disturb then test for failure
        for i in range(0, degree // 3):
            codeword[i] = field.zero

        proof_stream = ProofStream()
        fri.prove(codeword, proof_stream)
        with pytest.raises(FriError):
            fri.verify(proof_stream)
