import os

from stark_anatomy.algebra import Field, FieldElement
from stark_anatomy.ntt import (
    fast_coset_divide,
    fast_coset_evaluate,
    fast_evaluate,
    fast_interpolate,
    fast_multiply,
    intt,
    ntt,
)
from stark_anatomy.univariate import Polynomial


class TestNtt:
    def test_ntt(self, field: Field):
        log_n = 8
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)

        coefficients = [field.sample(os.urandom(17)) for i in range(n)]
        poly = Polynomial(coefficients)

        values = ntt(primitive_root, coefficients)

        values_again = poly.evaluate_domain(
            [primitive_root ^ i for i in range(len(values))]
        )

        assert values == values_again

    def test_intt(self, field: Field):
        log_n = 7
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)

        values = [field.sample(os.urandom(1)) for i in range(n)]
        coefficients = ntt(primitive_root, values)
        values_again = intt(primitive_root, coefficients)

        assert values == values_again

    def test_multiply(self, field: Field):
        log_n = 6
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)

        for _ in range(20):
            lhs_degree = int(os.urandom(1)[0]) % (n // 2)
            rhs_degree = int(os.urandom(1)[0]) % (n // 2)

            lhs = Polynomial(
                [field.sample(os.urandom(17)) for _ in range(lhs_degree + 1)]
            )
            rhs = Polynomial(
                [field.sample(os.urandom(17)) for _ in range(rhs_degree + 1)]
            )

            fast_product = fast_multiply(lhs, rhs, primitive_root, n)
            slow_product = lhs * rhs

            assert fast_product == slow_product

    def test_divide(self, field: Field):
        log_n = 6
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)

        for _ in range(20):
            lhs_degree = int(os.urandom(1)[0]) % (n // 2)
            rhs_degree = int(os.urandom(1)[0]) % (n // 2)

            lhs = Polynomial(
                [field.sample(os.urandom(17)) for i in range(lhs_degree + 1)]
            )
            rhs = Polynomial(
                [field.sample(os.urandom(17)) for i in range(rhs_degree + 1)]
            )

            fast_product = fast_multiply(lhs, rhs, primitive_root, n)
            quotient = fast_coset_divide(
                fast_product, lhs, field.generator(), primitive_root, n
            )

            assert quotient == rhs

    def test_interpolate(self, field: Field):
        log_n = 9
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)

        for _ in range(10):
            N = sum((1 << (8 * i)) * int(os.urandom(1)[0]) for i in range(8)) % n
            if N == 0:
                continue
            print("N:", N)
            values = [field.sample(os.urandom(17)) for _ in range(N)]
            domain = [field.sample(os.urandom(17)) for _ in range(N)]
            poly = fast_interpolate(domain, values, primitive_root, n)
            print("poly degree:", poly.degree())
            values_again = fast_evaluate(poly, domain, primitive_root, n)[0:N]

            assert values == values_again

    def test_coset_evaluate(self, field: Field):
        log_n = 9
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)

        two = FieldElement(2, field)

        domain = [two * (primitive_root ^ i) for i in range(n)]

        degree = ((int(os.urandom(1)[0]) * 256 + int(os.urandom(1)[0])) % n) - 1
        coefficients = [field.sample(os.urandom(17)) for i in range(degree + 1)]
        poly = Polynomial(coefficients)

        values_fast = fast_coset_evaluate(poly, two, primitive_root, n)
        values_traditional = [poly.evaluate(d) for d in domain]

        assert all(vf == vt for (vf, vt) in zip(values_fast, values_traditional))
