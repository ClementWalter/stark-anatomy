from hypothesis import assume, given
from hypothesis import strategies as st

from stark_anatomy.algebra import Field, FieldElement


class TestField:
    def test_field_properties(self, field: Field):
        assert field.zero == 0
        assert field.one == 1

    @given(a=..., b=..., c=...)
    def test_field_addition(self, a: FieldElement, b: FieldElement, c: FieldElement):
        # Commutativity
        assert a + b == b + a

        # Associativity
        assert (a + b) + c == a + (b + c)

        # Additive inverse
        assert a + (-a) == 0
        assert b + (-b) == 0

    @given(a=..., b=..., c=...)
    def test_field_multiplication(
        self, a: FieldElement, b: FieldElement, c: FieldElement
    ):
        # Commutativity
        assert a * b == b * a

        # Associativity
        assert (a * b) * c == a * (b * c)

        # Distributivity
        assert a * (b + c) == (a * b) + (a * c)

    @given(a=...)
    def test_field_inverse(self, a: FieldElement):
        assume(a != 0)

        assert a * (1 / a) == 1
        assert (1 / a) * a == 1

    @given(a=..., exponent=st.integers(min_value=0, max_value=10000))
    def test_field_exponentiation(self, a: FieldElement, exponent: int):
        assume(a != 0)

        assert a**exponent == pow(a.value, exponent, a.field.p)

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
