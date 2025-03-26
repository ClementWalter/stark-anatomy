from hypothesis import assume, given, settings
from hypothesis import strategies as st

from stark_anatomy.algebra import FieldElement
from stark_anatomy.univariate import Polynomial


class TestUnivariate:
    @given(a=..., b=..., c=...)
    def test_distributivity(self, a: Polynomial, b: Polynomial, c: Polynomial):
        assert a * (b + c) == a * b + a * c == (b + c) * a

    @given(a=..., b=..., c=...)
    def test_division(self, a: Polynomial, b: Polynomial, c: Polynomial):
        assume(a != 0)
        assume(b != 0)
        assume(c != 0)

        quo, rem = divmod(a * b, a)
        assert rem == 0
        assert quo == b

        quo, rem = divmod(a * b, b)
        assert rem == 0
        assert quo == a

        quo, rem = divmod(a * b, c)
        assert quo * c + rem == a * b

    @given(size=st.integers(min_value=1, max_value=10), data=st.data())
    def test_interpolate(self, size: int, data: st.DataObject):
        domain = data.draw(
            st.sets(st.from_type(FieldElement), min_size=size, max_size=size)
        )
        values = data.draw(
            st.sets(st.from_type(FieldElement), min_size=size, max_size=size)
        )

        poly = Polynomial.interpolate(list(domain), list(values))
        assert all(poly.evaluate(d) == v for d, v in zip(domain, values))

    @given(domain=st.sets(st.from_type(FieldElement), min_size=1, max_size=255))
    def test_zerofier(self, domain: set[FieldElement]):
        zerofier = Polynomial.zerofier_domain(list(domain))

        assert zerofier.degree == len(domain)
        assert all(zerofier.evaluate(d) == 0 for d in domain)
