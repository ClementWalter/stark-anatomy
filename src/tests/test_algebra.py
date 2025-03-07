import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from stark_anatomy.algebra import Field, FieldElement


class TestField:
    def test_field_properties(self, field: Field):
        zero = field.zero
        one = field.one

        # Zero is the additive identity
        assert zero + zero == zero
        assert one + zero == one
        assert zero + one == one

        # One is the multiplicative identity
        assert one * one == one
        assert zero * one == zero
        assert one * zero == zero

        # Zero has no multiplicative inverse
        with pytest.raises(ValueError):
            zero.inverse()

    @given(value=...)
    def test_field_element_creation(self, value: int, field: Field):
        element = FieldElement(value, field)
        assert element.value % field.p == value % field.p
        assert element.field == field

    @given(a=..., b=...)
    def test_field_addition(self, a: int, b: int, field: Field):
        elem_a = FieldElement(a, field)
        elem_b = FieldElement(b, field)

        # Commutativity
        assert elem_a + elem_b == elem_b + elem_a

        # Associativity
        elem_c = FieldElement(5, field)
        assert (elem_a + elem_b) + elem_c == elem_a + (elem_b + elem_c)

        # Additive inverse
        assert elem_a + (-elem_a) == field.zero
        assert elem_b + (-elem_b) == field.zero

    @given(a=..., b=...)
    def test_field_multiplication(self, a: int, b: int, field: Field):
        elem_a = FieldElement(a, field)
        elem_b = FieldElement(b, field)

        # Commutativity
        assert elem_a * elem_b == elem_b * elem_a

        # Associativity
        elem_c = FieldElement(7, field)
        assert (elem_a * elem_b) * elem_c == elem_a * (elem_b * elem_c)

        # Distributivity
        assert elem_a * (elem_b + elem_c) == (elem_a * elem_b) + (elem_a * elem_c)

    @given(a=...)
    def test_field_inverse(self, a: int, field: Field):
        assume(a % field.p != 0)
        elem_a = FieldElement(a, field)
        inverse_a = elem_a.inverse()

        assert elem_a * inverse_a == field.one
        assert inverse_a * elem_a == field.one

    @given(a=..., exponent=st.integers(min_value=0, max_value=10000))
    def test_field_exponentiation(self, a: int, exponent: int, field: Field):
        assume(a % field.p != 0)
        elem_a = FieldElement(a, field)

        assert (elem_a**exponent).value % field.p == pow(a, exponent, field.p)

    def test_field_default(self, field: Field):
        assert field.p == Field.PRIME
        assert isinstance(field.zero, FieldElement)
        assert isinstance(field.one, FieldElement)
        assert field.zero.value == 0
        assert field.one.value == 1

    def test_field_generator(self, field: Field):
        generator = field.generator()
        assert isinstance(generator, FieldElement)
        assert generator.field == field

        # Test that generator is not 0 or 1
        assert generator != field.zero
        assert generator != field.one

    @given(n=st.integers(min_value=1, max_value=100))
    def test_primitive_nth_root(self, n: int, field: Field):
        n = 1 << n  # Make n a power of 2
        root = field.primitive_nth_root(n)
        assert root**n == field.one

    @given(byte_array=st.binary(min_size=1, max_size=100))
    def test_sample(self, field: Field, byte_array: bytes):
        element = field.sample(byte_array)
        assert isinstance(element, FieldElement)
        assert element.field == field
        assert 0 <= element.value < field.p
