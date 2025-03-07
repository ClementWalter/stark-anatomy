from typing import Union

from stark_anatomy.algebra import Field, FieldElement
from stark_anatomy.univariate import Polynomial


class MPolynomial:
    """
    Multivariate polynomials are represented as dictionaries with exponent vectors
    as keys and coefficients as values. E.g.:
    f(x,y,z) = 17 + 2xy + 42z - 19x^6*y^3*z^12 is represented as:
    {
        (0,0,0) => 17,
        (1,1,0) => 2,
        (0,0,1) => 42,
        (6,3,12) => -19,
    }
    """

    MPolynomialDict = dict[tuple[int, ...], FieldElement]

    def __init__(self, dictionary: MPolynomialDict):
        if dictionary and len({len(k) for k in dictionary.keys()}) != 1:
            raise ValueError("All keys must have the same length.")

        self.dictionary = dictionary

    def __neg__(self) -> "MPolynomial":
        return MPolynomial({k: -v for k, v in self.dictionary.items()})

    def __eq__(self, other: Union["MPolynomial", int, FieldElement]) -> bool:
        if isinstance(other, int):
            if not self.dictionary:
                return other == 0
            other = FieldElement(other, list(self.dictionary.values())[0].field)
        if isinstance(other, FieldElement):
            other = MPolynomial({(0,): other})
        if set(self.dictionary.keys()) != set(other.dictionary.keys()):
            return False

        return all(
            self.dictionary[k] == other.dictionary[k] for k in self.dictionary.keys()
        )

    def __ne__(self, other: "MPolynomial") -> bool:
        return not self.__eq__(other)

    def __add__(self, other: Union["MPolynomial", int, FieldElement]) -> "MPolynomial":
        if isinstance(other, int):
            other = FieldElement(other, list(self.dictionary.values())[0].field)
        if isinstance(other, FieldElement):
            other = MPolynomial({(0,): other})

        if self == 0:
            return other
        elif other == 0:
            return self

        dictionary = dict()
        num_variables = max(self.num_variables, other.num_variables)
        for k, v in self.dictionary.items():
            pad = list(k) + [0] * (num_variables - len(k))
            pad = tuple(pad)
            dictionary[pad] = v
        for k, v in other.dictionary.items():
            pad = list(k) + [0] * (num_variables - len(k))
            pad = tuple(pad)
            if pad in dictionary.keys():
                dictionary[pad] = dictionary[pad] + v
            else:
                dictionary[pad] = v

        return MPolynomial(dictionary).trim()

    def __sub__(self, other: Union["MPolynomial", int, FieldElement]) -> "MPolynomial":
        return self.__add__(-other)

    def __mul__(self, other: Union["MPolynomial", int, FieldElement]) -> "MPolynomial":
        if isinstance(other, int):
            other = FieldElement(other, list(self.dictionary.values())[0].field)
        if isinstance(other, FieldElement):
            other = MPolynomial({(0,): other})

        if self == 0 or other == 0:
            return MPolynomial({})

        dictionary = dict()
        num_variables = max(self.num_variables, other.num_variables)
        for k0, v0 in self.dictionary.items():
            for k1, v1 in other.dictionary.items():
                exponent = [0] * num_variables
                for k in range(len(k0)):
                    exponent[k] += k0[k]
                for k in range(len(k1)):
                    exponent[k] += k1[k]
                exponent = tuple(exponent)
                if exponent in dictionary.keys():
                    dictionary[exponent] = dictionary[exponent] + v0 * v1
                else:
                    dictionary[exponent] = v0 * v1
        return MPolynomial(dictionary).trim()

    def __pow__(self, exponent: int) -> "MPolynomial":
        if self == 0:
            return MPolynomial({})

        field = list(self.dictionary.values())[0].field
        acc = MPolynomial({(0,) * self.num_variables: field.one})
        for b in bin(exponent)[2:]:
            acc = acc * acc
            if b == "1":
                acc = acc * self
        return acc.trim()

    def __repr__(self) -> str:
        return "{" + ", ".join(f"{k}: {v}" for k, v in self.dictionary.items()) + "}"

    @property
    def num_variables(self) -> int:
        return len(list(self.dictionary.keys())[0]) if self.dictionary else 0

    def trim(self) -> "MPolynomial":
        """
        Trims the polynomial by removing zero coefficients.
        """
        return MPolynomial({k: v for k, v in self.dictionary.items() if v != 0})

    def evaluate(self, point: tuple[FieldElement, ...]) -> FieldElement:
        if len(point) != self.num_variables:
            raise ValueError(
                "number of variables in point does not match number of variables in polynomial"
            )

        if len(point) == 0:
            raise ValueError("point must have at least one variable")

        acc = point[0].field.zero
        for k, v in self.dictionary.items():
            prod = v
            for i in range(len(k)):
                prod = prod * (point[i] ** k[i])
            acc = acc + prod
        return acc

    def evaluate_symbolic(self, point: tuple[FieldElement, ...]) -> Polynomial:
        acc = Polynomial([])
        for k, v in self.dictionary.items():
            prod = Polynomial([v])
            for i in range(len(k)):
                prod = prod * (point[i] ** k[i])
            acc = acc + prod
        return acc

    @staticmethod
    def variables(num_variables: int, field: Field) -> list["MPolynomial"]:
        """
        Creates a list of the single variables polynomials.

        For a given number of variables, returns a list of MPolynomial objects,
        where each polynomial represents a single variable in the multivariate system.

        For example, with num_variables=3, this returns [x, y, z] as multivariate polynomials, where:
        - x is the polynomial representing the first variable (1*x + 0*y + 0*z)
        - y is the polynomial representing the second variable (0*x + 1*y + 0*z)
        - z is the polynomial representing the third variable (0*x + 0*y + 1*z)

        Parameters:
            num_variables: The number of variables to create
            field: The field to use for the coefficients

        Returns:
            A list of MPolynomial objects, each representing a single variable
        """
        variables = []
        for i in range(num_variables):
            exponent = [0] * num_variables
            exponent[i] = 1
            variables = variables + [MPolynomial({tuple(exponent): field.one})]
        return variables

    @staticmethod
    def lift(polynomial: Polynomial, variable_index: int) -> "MPolynomial":
        """
        Lift a univariate polynomial to a multivariate polynomial.

        This function takes a univariate polynomial and lifts it to a multivariate polynomial
        where the single variable in the univariate polynomial becomes the `variable_index`
        variable of the multivariate polynomial.

        Parameters:
            polynomial: The univariate polynomial to lift
            variable_index: The index of the variable to lift

        Returns:
            The lifted multivariate polynomial, with `variable_index + 1` variables
        """
        if polynomial == 0:
            return MPolynomial({})

        return MPolynomial(
            {
                (0,) * variable_index + (i,): coefficient
                for i, coefficient in enumerate(polynomial.coefficients)
            }
        )

    @staticmethod
    def constant(value: FieldElement) -> "MPolynomial":
        return MPolynomial({(): value})
