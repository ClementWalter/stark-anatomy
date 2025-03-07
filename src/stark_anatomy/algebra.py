def xgcd(x, y):
    """
    Extended Euclidean Algorithm, see https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    """
    old_r, r = (x, y)
    old_s, s = (1, 0)
    old_t, t = (0, 1)

    while r != 0:
        quotient = old_r // r
        old_r, r = (r, old_r - quotient * r)
        old_s, s = (s, old_s - quotient * s)
        old_t, t = (t, old_t - quotient * t)

    return old_s, old_t, old_r  # a, b, g


class FieldElement:
    def __init__(self, value: int, field: "Field"):
        self.value = value
        self.field = field

    def __add__(self, right: "FieldElement") -> "FieldElement":
        return self.field.add(self, right)

    def __mul__(self, right: "FieldElement") -> "FieldElement":
        return self.field.multiply(self, right)

    def __sub__(self, right: "FieldElement") -> "FieldElement":
        return self.field.subtract(self, right)

    def __truediv__(self, right: "FieldElement") -> "FieldElement":
        return self.field.divide(self, right)

    def __neg__(self) -> "FieldElement":
        return self.field.negate(self)

    def inverse(self) -> "FieldElement":
        return self.field.inverse(self)

    def __pow__(self, exponent: int) -> "FieldElement":
        acc = self.field.one
        val = self
        for i in bin(exponent)[2:]:
            acc *= acc
            if i != "0":
                acc *= val
        return acc

    def __eq__(self, other: "FieldElement") -> bool:
        return self.field.eq(self, other)

    def __neq__(self, other: "FieldElement") -> bool:
        return self.field.neq(self, other)

    def __str__(self) -> str:
        return str(self.value)

    def __bytes__(self) -> bytes:
        return bytes(str(self).encode())

    def is_zero(self) -> bool:
        return self.field.eq(self, self.field.zero)


class Field:
    PRIME = 1 + 407 * (1 << 119)

    def __init__(self, p: int):
        self.p = p

    @property
    def zero(self) -> FieldElement:
        return FieldElement(0, self)

    @property
    def one(self) -> FieldElement:
        return FieldElement(1, self)

    def multiply(self, left: FieldElement, right: FieldElement) -> FieldElement:
        return FieldElement((left.value * right.value) % self.p, self)

    def eq(self, left: FieldElement, right: FieldElement) -> bool:
        return left.value % self.p == right.value % self.p

    def neq(self, left: FieldElement, right: FieldElement) -> bool:
        return left.value % self.p != right.value % self.p

    def add(self, left: FieldElement, right: FieldElement) -> FieldElement:
        return FieldElement((left.value + right.value) % self.p, self)

    def subtract(self, left: FieldElement, right: FieldElement) -> FieldElement:
        return FieldElement((self.p + left.value - right.value) % self.p, self)

    def negate(self, operand: FieldElement) -> FieldElement:
        return FieldElement((self.p - operand.value) % self.p, self)

    def inverse(self, operand: FieldElement) -> FieldElement:
        if operand.is_zero():
            raise ValueError("Cannot invert zero")
        a, *_ = xgcd(operand.value, self.p)
        return FieldElement(((a % self.p) + self.p) % self.p, self)

    def divide(self, left: FieldElement, right: FieldElement) -> FieldElement:
        if right.is_zero():
            raise ValueError("Cannot divide by zero")
        a, *_ = xgcd(right.value, self.p)
        return FieldElement(left.value * a % self.p, self)

    @staticmethod
    def default() -> "Field":
        return Field(Field.PRIME)

    def generator(self) -> FieldElement:
        assert (
            self.p == Field.PRIME
        ), "Do not know generator for other fields beyond 1+407*2^119"
        return FieldElement(85408008396924667383611388730472331217, self)

    def primitive_nth_root(self, n: int) -> FieldElement:
        assert self.p == Field.PRIME, "Unknown field, can't return root of unity."
        assert (
            n <= 1 << 119 and (n & (n - 1)) == 0
        ), "Field does not have nth root of unity where n > 2^119 or not power of two."
        root = FieldElement(85408008396924667383611388730472331217, self)
        order = 1 << 119
        while order != n:
            root = root**2
            order = order // 2
        return root

    def sample(self, byte_array: bytes) -> FieldElement:
        acc = 0
        for b in byte_array:
            acc = (acc << 8) ^ int(b)
        return FieldElement(acc % self.p, self)
