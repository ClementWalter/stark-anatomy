from hypothesis import assume, given, settings
from hypothesis import strategies as st

from stark_anatomy.algebra import Field, FieldElement
from stark_anatomy.ntt import (
    fast_coset_divide,
    fast_coset_evaluate,
    fast_evaluate,
    fast_interpolate,
    fast_multiply,
    fast_zerofier,
    intt,
    ntt,
)
from stark_anatomy.univariate import Polynomial


class TestNtt:
    @given(log_n=st.integers(min_value=1, max_value=8), data=st.data())
    def test_ntt(self, field: Field, log_n: int, data: st.DataObject):
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)
        domain = data.draw(st.lists(st.from_type(FieldElement), min_size=n, max_size=n))

        assert ntt(primitive_root, domain) == [
            Polynomial(domain).evaluate(primitive_root**i) for i in range(n)
        ]

    @given(log_n=st.integers(min_value=1, max_value=8), data=st.data())
    def test_intt(self, field: Field, log_n: int, data: st.DataObject):
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)
        domain = data.draw(st.lists(st.from_type(FieldElement), min_size=n, max_size=n))

        assert domain == intt(primitive_root, ntt(primitive_root, domain))

    @given(lhs=..., rhs=...)
    def test_fast_multiply(self, field: Field, lhs: Polynomial, rhs: Polynomial):
        assume(lhs != 0 and rhs != 0)
        n = 2 ** (lhs.degree + rhs.degree + 1)
        primitive_root = field.primitive_nth_root(n)

        assert fast_multiply(lhs, rhs, primitive_root, n) == lhs * rhs

    @given(log_n=st.integers(min_value=1, max_value=6), data=st.data())
    def test_fast_zerofier(self, field: Field, log_n: int, data: st.DataObject):
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)
        domain = data.draw(
            st.lists(st.from_type(FieldElement), min_size=n - 1, max_size=n - 1)
        )

        assert fast_zerofier(domain, primitive_root, n) == Polynomial.zerofier_domain(
            domain
        )

    @given(log_n=st.integers(min_value=1, max_value=4), poly=..., data=st.data())
    def test_fast_evaluate(
        self, field: Field, log_n: int, poly: Polynomial, data: st.DataObject
    ):
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)
        domain = data.draw(st.lists(st.from_type(FieldElement), min_size=n, max_size=n))

        assert fast_evaluate(poly, domain, primitive_root, n) == [
            poly.evaluate(d) for d in domain
        ]

    @given(log_n=st.integers(min_value=2, max_value=4), data=st.data())
    def test_fast_interpolate(self, field: Field, log_n: int, data: st.DataObject):
        n = 1 << log_n
        primitive_root = field.primitive_nth_root(n)
        domain = list(
            # Domain items must be unique to interpolate
            data.draw(st.sets(st.from_type(FieldElement), min_size=n, max_size=n))
        )
        values = data.draw(st.lists(st.from_type(FieldElement), min_size=n, max_size=n))

        assert fast_interpolate(
            domain, values, primitive_root, n
        ) == Polynomial.interpolate(domain, values)

    @given(log_n=st.integers(min_value=2, max_value=4), offset=..., data=st.data())
    def test_fast_coset_evaluate(
        self, field: Field, data: st.DataObject, log_n: int, offset: FieldElement
    ):
        assume(offset != 0)
        n = 1 << log_n
        generator = field.primitive_nth_root(n)
        domain = data.draw(st.lists(st.from_type(FieldElement), min_size=n, max_size=n))

        poly = Polynomial(domain)

        assert fast_coset_evaluate(poly, offset, generator, n) == [
            poly.scale(offset).evaluate(generator**i) for i in range(n)
        ]

    @given(
        lhs_degree=st.integers(min_value=2, max_value=4),
        rhs_degree=st.integers(min_value=2, max_value=4),
        data=st.data(),
    )
    def test_fast_coset_divide(
        self, field: Field, data: st.DataObject, lhs_degree: int, rhs_degree: int
    ):
        n = 1 << (lhs_degree + rhs_degree + 1)
        primitive_root = field.primitive_nth_root(n)

        rhs = Polynomial(
            data.draw(
                st.lists(
                    st.from_type(FieldElement),
                    min_size=rhs_degree + 1,
                    max_size=rhs_degree + 1,
                )
            )
        )
        assume(rhs != 0)

        lhs = (
            Polynomial(
                data.draw(
                    st.lists(
                        st.from_type(FieldElement),
                        min_size=lhs_degree + 1,
                        max_size=lhs_degree + 1,
                    )
                )
            )
            * rhs
        )

        assert lhs / rhs == fast_coset_divide(
            lhs, rhs, field.generator(), primitive_root, n
        )
