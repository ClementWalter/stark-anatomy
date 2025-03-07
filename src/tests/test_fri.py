from stark_anatomy.algebra import Field, FieldElement
from stark_anatomy.fri import Fri
from stark_anatomy.ip import ProofStream
from stark_anatomy.univariate import Polynomial


class TestFri:
    def test_fri(self, field: Field):
        degree = 63
        expansion_factor = 4
        num_colinearity_tests = 17

        initial_codeword_length = (degree + 1) * expansion_factor
        log_codeword_length = 0
        codeword_length = initial_codeword_length
        while codeword_length > 1:
            codeword_length //= 2
            log_codeword_length += 1

        assert 1 << log_codeword_length == initial_codeword_length

        omega = field.primitive_nth_root(initial_codeword_length)
        generator = field.generator()

        assert omega ^ (1 << log_codeword_length) == field.one
        assert omega ^ (1 << (log_codeword_length - 1)) != field.one

        fri = Fri(
            generator,
            omega,
            initial_codeword_length,
            expansion_factor,
            num_colinearity_tests,
        )

        polynomial = Polynomial([FieldElement(i, field) for i in range(degree + 1)])
        domain = [omega ^ i for i in range(initial_codeword_length)]

        codeword = polynomial.evaluate_domain(domain)

        # test valid codeword
        proof_stream = ProofStream()
        fri.prove(codeword, proof_stream)
        points = []
        assert fri.verify(proof_stream, points)

        assert all(polynomial.evaluate(omega ^ x) == y for x, y in points)

        # disturb then test for failure
        for i in range(0, degree // 3):
            codeword[i] = field.zero

        proof_stream = ProofStream()
        fri.prove(codeword, proof_stream)
        points = []
        assert not fri.verify(proof_stream, points)
