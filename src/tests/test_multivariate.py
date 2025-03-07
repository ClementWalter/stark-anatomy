from stark_anatomy.algebra import Field, FieldElement
from stark_anatomy.multivariate import MPolynomial
from stark_anatomy.univariate import Polynomial


class TestMultivariate:
    def test_evaluate(self, field: Field):
        variables = MPolynomial.variables(4, field)
        zero = field.zero
        one = field.one
        two = FieldElement(2, field)
        five = FieldElement(5, field)

        mpoly1 = (
            MPolynomial.constant(one) * variables[0]
            + MPolynomial.constant(two) * variables[1]
            + MPolynomial.constant(five) * (variables[2] ^ 3)
        )
        mpoly2 = (
            MPolynomial.constant(one) * variables[0] * variables[3]
            + MPolynomial.constant(five) * (variables[3] ^ 3)
            + MPolynomial.constant(five)
        )

        mpoly3 = mpoly1 * mpoly2

        point = [zero, five, five, two]

        eval1 = mpoly1.evaluate(point)
        eval2 = mpoly2.evaluate(point)
        eval3 = mpoly3.evaluate(point)

        assert eval1 * eval2 == eval3
        assert eval1 + eval2 == (mpoly1 + mpoly2).evaluate(point)

    def test_lift(self, field: Field):
        zero = field.zero
        one = field.one
        two = FieldElement(2, field)
        five = FieldElement(5, field)

        univariate_poly = Polynomial.interpolate_domain(
            [zero, one, two], [two, five, five]
        )
        mpoly = MPolynomial.lift(univariate_poly, 3)

        assert univariate_poly.evaluate(five) == mpoly.evaluate(
            [zero, zero, zero, five]
        )
