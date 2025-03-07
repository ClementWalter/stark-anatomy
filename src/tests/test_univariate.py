import os

from hypothesis import given
from hypothesis import strategies as st

from stark_anatomy.algebra import Field, FieldElement
from stark_anatomy.univariate import Polynomial


class TestUnivariate:
    def test_distributivity(self):
        field = Field.default()
        zero = field.zero
        one = field.one
        two = FieldElement(2, field)
        five = FieldElement(5, field)

        a = Polynomial([one, zero, five, two])
        b = Polynomial([two, two, one])
        c = Polynomial([zero, five, two, five, five, one])

        lhs = a * (b + c)
        rhs = a * b + a * c
        assert lhs == rhs

    def test_division(self):
        field = Field.default()
        zero = field.zero
        one = field.one
        two = FieldElement(2, field)
        five = FieldElement(5, field)

        a = Polynomial([one, zero, five, two])
        b = Polynomial([two, two, one])
        c = Polynomial([zero, five, two, five, five, one])

        # a should divide a*b, quotient should be b
        quo, rem = Polynomial.divide(a * b, a)
        assert rem.is_zero()
        assert quo == b

        # b should divide a*b, quotient should be a
        quo, rem = Polynomial.divide(a * b, b)
        assert rem.is_zero()
        assert quo == a

        # c should not divide a*b
        quo, rem = Polynomial.divide(a * b, c)
        assert not rem.is_zero()
        assert quo * c + rem == a * b

    def test_interpolate(self):
        field = Field.default()
        zero = field.zero
        one = field.one
        two = FieldElement(2, field)
        five = FieldElement(5, field)

        values = [five, two, two, one, five, zero]
        domain = [FieldElement(i, field) for i in range(6)]

        poly = Polynomial.interpolate_domain(domain, values)

        assert all(poly.evaluate(d) == v for d, v in zip(domain, values))

        # evaluation in random point is nonzero with high probability
        assert poly.evaluate(FieldElement(363, field)) != zero
        assert poly.degree() == len(domain) - 1

    @given(degree=st.integers(min_value=1, max_value=255))
    def test_zerofier(self, degree: int):
        field = Field.default()

        domain = []
        while len(domain) != degree:
            new = field.sample(os.urandom(17))
            if new not in domain:
                domain += [new]

        zerofier = Polynomial.zerofier_domain(domain)

        assert zerofier.degree() == degree
        assert all(zerofier.evaluate(d) == field.zero for d in domain)

        random = field.sample(os.urandom(17))
        while random in domain:
            random = field.sample(os.urandom(17))

        assert zerofier.evaluate(random) != field.zero
