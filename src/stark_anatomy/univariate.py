from typing import List, Sequence, Tuple, Union

from stark_anatomy.algebra import FieldElement


class Polynomial:
    def __init__(self, coefficients: list[FieldElement]):
        self.coefficients = coefficients.copy()
        self.trim()

    def __neg__(self) -> "Polynomial":
        return Polynomial([-c for c in self.coefficients])

    def __eq__(self, other: Union["Polynomial", int, FieldElement]) -> bool:
        if isinstance(other, int):
            if not self.coefficients:
                return other == 0
            other = FieldElement(other, self.coefficients[0].field)
        if isinstance(other, FieldElement):
            other = Polynomial([other])
        if len(self.coefficients) != len(other.coefficients):
            return False

        return all(
            self.coefficients[i] == other.coefficients[i]
            for i in range(len(self.coefficients))
        )

    def __ne__(self, other: "Polynomial") -> bool:
        return not self.__eq__(other)

    def __add__(self, other: Union["Polynomial", int, FieldElement]) -> "Polynomial":
        if isinstance(other, int):
            other = FieldElement(other, self.coefficients[0].field)
        if isinstance(other, FieldElement):
            other = Polynomial([other]).trim()

        if self == 0:
            return other
        elif other == 0:
            return self

        field = self.coefficients[0].field
        coeffs = [field.zero] * max(len(self.coefficients), len(other.coefficients))
        for i in range(len(self.coefficients)):
            coeffs[i] = coeffs[i] + self.coefficients[i]
        for i in range(len(other.coefficients)):
            coeffs[i] = coeffs[i] + other.coefficients[i]

        return Polynomial(coeffs).trim()

    def __sub__(self, other: Union["Polynomial", int, FieldElement]) -> "Polynomial":
        return self.__add__(-other)

    def __mul__(self, other: Union["Polynomial", int, FieldElement]) -> "Polynomial":
        if isinstance(other, int):
            other = FieldElement(other, self.coefficients[0].field)
        if isinstance(other, FieldElement):
            other = Polynomial([other])

        if self == 0 or other == 0:
            return Polynomial([])

        zero = self.coefficients[0].field.zero
        buf = [zero] * (len(self.coefficients) + len(other.coefficients) - 1)
        for i in range(len(self.coefficients)):
            if self.coefficients[i] == 0:
                continue
            for j in range(len(other.coefficients)):
                buf[i + j] = buf[i + j] + self.coefficients[i] * other.coefficients[j]
        return Polynomial(buf).trim()

    def __divmod__(
        self, other: Union["Polynomial", int, FieldElement]
    ) -> Tuple["Polynomial", "Polynomial"]:
        if isinstance(other, int):
            other = FieldElement(other, self.coefficients[0].field)
        if isinstance(other, FieldElement):
            other = Polynomial([other])

        if other.degree == -1:
            raise ValueError("divide by zero")
        if self.degree < other.degree:
            return Polynomial([]), self

        zero = other.coefficients[0].field.zero
        remainder = Polynomial([n for n in self.coefficients])
        quotient_coefficients = [zero for _ in range(self.degree - other.degree + 1)]
        for _ in range(self.degree - other.degree + 1):
            if remainder.degree < other.degree:
                break
            coefficient = remainder.leading_coefficient / other.leading_coefficient
            shift = remainder.degree - other.degree
            subtractee = Polynomial([zero] * shift + [coefficient]) * other
            quotient_coefficients[shift] = coefficient
            remainder = remainder - subtractee
        quotient = Polynomial(quotient_coefficients)
        return quotient.trim(), remainder.trim()

    def __truediv__(
        self, other: Union["Polynomial", int, FieldElement]
    ) -> "Polynomial":
        if isinstance(other, int):
            other = FieldElement(other, self.coefficients[0].field)
        if isinstance(other, FieldElement):
            other = Polynomial([other])

        quo, rem = self.__divmod__(other)
        assert (
            rem == 0
        ), "cannot perform polynomial division because remainder is not zero"
        return quo.trim()

    def __mod__(self, other: Union["Polynomial", int, FieldElement]) -> "Polynomial":
        if isinstance(other, int):
            other = FieldElement(other, self.coefficients[0].field)
        if isinstance(other, FieldElement):
            other = Polynomial([other])

        _, rem = self.__divmod__(other)
        return rem.trim()

    def __pow__(self, exponent: int) -> "Polynomial":
        if self == 0:
            return Polynomial([])
        if exponent == 0:
            return Polynomial([self.coefficients[0].field.one])
        acc = Polynomial([self.coefficients[0].field.one])
        for b in bin(exponent)[2:]:
            acc = acc * acc
            if b == "1":
                acc = acc * self
        return acc.trim()

    def __repr__(self) -> str:
        return "[" + ",".join(str(s) for s in self.coefficients) + "]"

    @property
    def leading_coefficient(self) -> FieldElement:
        return self.coefficients[-1]

    @property
    def degree(self) -> int:
        """
        Returns -1 if the polynomial is the zero polynomial.
        """
        return len(self.coefficients) - 1

    def trim(self) -> "Polynomial":
        """
        Trims the polynomial by removing trailing zero coefficients.
        """
        while self.coefficients and self.coefficients[-1] == 0:
            self.coefficients.pop()
        return self

    def evaluate(self, point: Union[FieldElement, int]) -> FieldElement:
        if isinstance(point, int):
            point = FieldElement(point, self.coefficients[0].field)

        xi = point.field.one
        value = point.field.zero
        for c in self.coefficients:
            value = value + c * xi
            xi = xi * point
        return value

    def evaluate_domain(
        self, domain: Sequence[Union[FieldElement, int]]
    ) -> List[FieldElement]:
        return [self.evaluate(d) for d in domain]

    def scale(self, factor: Union[int, FieldElement]) -> "Polynomial":
        """
        Scales the polynomial by _factor_, i.e. actually creates a new polynomial Q such that
        Q(x) = P(factor * x) for all x.
        """
        return Polynomial([c * factor**i for i, c in enumerate(self.coefficients)])

    @staticmethod
    def interpolate(
        domain: List[FieldElement], values: List[FieldElement]
    ) -> "Polynomial":
        if len(domain) != len(values):
            raise ValueError(
                "number of elements in domain does not match number of values -- cannot interpolate"
            )
        if len(domain) == 0:
            raise ValueError("cannot interpolate between zero points")

        field = domain[0].field
        x = Polynomial([field.zero, field.one])
        acc = Polynomial([])
        for i in range(len(domain)):
            prod = Polynomial([values[i]])
            for j in range(len(domain)):
                if j == i:
                    continue
                prod = prod * (x - domain[j]) / (domain[i] - domain[j])
            acc = acc + prod
        return acc.trim()

    @staticmethod
    def zerofier_domain(domain: List[FieldElement]) -> "Polynomial":
        field = domain[0].field
        x = Polynomial([field.zero, field.one])
        acc = Polynomial([field.one])
        for d in domain:
            acc *= x - Polynomial([d])
        return acc.trim()

    @staticmethod
    def is_colinear(points: List[Tuple[FieldElement, FieldElement]]) -> bool:
        domain = [p[0] for p in points]
        values = [p[1] for p in points]
        polynomial = Polynomial.interpolate(domain, values)
        return polynomial.degree == 1
