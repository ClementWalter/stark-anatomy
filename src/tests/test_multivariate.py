import math

from hypothesis import assume, given
from hypothesis import strategies as st

from stark_anatomy.algebra import FieldElement
from stark_anatomy.multivariate import MPolynomial
from stark_anatomy.univariate import Polynomial


class TestMultivariate:
    @given(a=..., b=..., c=...)
    def test_distributivity(self, a: MPolynomial, b: MPolynomial, c: MPolynomial):
        assert a * (b + c) == a * b + a * c == (b + c) * a

    @given(mpoly_1=..., mpoly_2=..., point=...)
    def test_evaluate(
        self,
        mpoly_1: MPolynomial,
        mpoly_2: MPolynomial,
        point: tuple[FieldElement, ...],
    ):
        assume(mpoly_1 != 0)
        assume(mpoly_2 != 0)
        assume(len(point) > 0)

        max_num_variables = max(mpoly_1.num_variables, mpoly_2.num_variables)
        point = point * math.ceil(max_num_variables / len(point))

        eval1 = mpoly_1.evaluate(point[: mpoly_1.num_variables])
        eval2 = mpoly_2.evaluate(point[: mpoly_2.num_variables])

        assert eval1 + eval2 == (mpoly_1 + mpoly_2).evaluate(point[:max_num_variables])

    @given(
        poly=...,
        mpoly=...,
        point=...,
        variable_index=st.integers(min_value=0, max_value=5),
    )
    def test_lift(
        self,
        poly: Polynomial,
        mpoly: MPolynomial,
        point: FieldElement,
        variable_index: int,
    ):
        assume(poly != 0)

        mpoly = MPolynomial.lift(poly, variable_index)

        assert poly.evaluate(point) == mpoly.evaluate((point,) * mpoly.num_variables)
