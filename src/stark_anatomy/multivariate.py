from stark_anatomy.algebra import FieldElement
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
        self.dictionary = dictionary

    @staticmethod
    def zero():
        return MPolynomial(dict())

    def __add__(self, other):
        dictionary = dict()
        num_variables = max(
            [len(k) for k in self.dictionary.keys()]
            + [len(k) for k in other.dictionary.keys()]
        )
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
        return MPolynomial(dictionary)

    def __mul__(self, other):
        dictionary = dict()
        num_variables = max(
            [len(k) for k in self.dictionary.keys()]
            + [len(k) for k in other.dictionary.keys()]
        )
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
        return MPolynomial(dictionary)

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        dictionary = dict()
        for k, v in self.dictionary.items():
            dictionary[k] = -v
        return MPolynomial(dictionary)

    def __xor__(self, exponent):
        if self.is_zero():
            return MPolynomial(dict())
        field = list(self.dictionary.values())[0].field
        num_variables = len(list(self.dictionary.keys())[0])
        exp = [0] * num_variables
        acc = MPolynomial({tuple(exp): field.one})
        for b in bin(exponent)[2:]:
            acc = acc * acc
            if b == "1":
                acc = acc * self
        return acc

    @staticmethod
    def constant(element):
        return MPolynomial({tuple([0]): element})

    def is_zero(self):
        if not self.dictionary:
            return True
        else:
            for v in self.dictionary.values():
                if not v.is_zero():
                    return False
            return True

    @staticmethod
    def variables(num_variables, field):
        """
        Creates a list of multivariate polynomial variables.

        For a given number of variables, returns a list of MPolynomial objects,
        where each polynomial represents a single variable in the multivariate system.

        For example, with num_variables=3, this returns:
        [x, y, z] as multivariate polynomials, where:
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
            exponent = [0] * i + [1] + [0] * (num_variables - i - 1)
            variables = variables + [MPolynomial({tuple(exponent): field.one})]
        return variables

    def evaluate(self, point):
        acc = point[0].field.zero
        for k, v in self.dictionary.items():
            prod = v
            for i in range(len(k)):
                prod = prod * (point[i] ^ k[i])
            acc = acc + prod
        return acc

    def evaluate_symbolic(self, point):
        acc = Polynomial([])
        for k, v in self.dictionary.items():
            prod = Polynomial([v])
            for i in range(len(k)):
                prod = prod * (point[i] ^ k[i])
            acc = acc + prod
        return acc

    @staticmethod
    def lift(polynomial: Polynomial, variable_index: int):
        if polynomial.is_zero():
            return MPolynomial({})

        field = polynomial.coefficients[0].field
        variables = MPolynomial.variables(variable_index + 1, field)
        x = variables[-1]
        acc = MPolynomial({})
        for i in range(len(polynomial.coefficients)):
            acc = acc + MPolynomial.constant(polynomial.coefficients[i]) * (x ^ i)
        return acc
