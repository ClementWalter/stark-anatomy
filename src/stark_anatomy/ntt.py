from stark_anatomy.algebra import FieldElement
from stark_anatomy.univariate import Polynomial


def ntt(
    primitive_root: FieldElement, coefficients: list[FieldElement]
) -> list[FieldElement]:
    """
    Computes the Number Theoretic Transform (NTT) of a sequence of field elements.

    The NTT is a specialized version of the Discrete Fourier Transform (DFT) that works
    over finite fields. It transforms a sequence of field elements from the coefficient
    representation to the point-value representation.

    This implementation uses a recursive divide-and-conquer approach (Cooley-Tukey algorithm).

    Parameters:
        primitive_root: A primitive nth root of unity in the field, where n is the length of coefficients
        coefficients: A list of field elements to transform, must be a power of two in length

    Returns:
        A list of field elements representing the NTT of the input sequence

    Raises:
        ValueError: If the length of coefficients is not a power of two or if primitive_root is not
                   a primitive nth root of unity
    """
    if len(coefficients) & (len(coefficients) - 1) != 0:
        raise ValueError("cannot compute ntt of non-power-of-two sequence")
    if len(coefficients) <= 1:
        return coefficients

    field = coefficients[0].field

    if primitive_root ** len(coefficients) != field.one:
        raise ValueError(
            "primitive root must be nth root of unity, where n is len(coefficients)"
        )
    if primitive_root ** (len(coefficients) // 2) == field.one:
        raise ValueError(
            "primitive root is not primitive nth root of unity, where n is len(coefficients)"
        )

    half = len(coefficients) // 2

    evens = ntt(primitive_root**2, coefficients[::2])
    odds = ntt(primitive_root**2, coefficients[1::2])

    return [
        evens[i % half] + (primitive_root**i) * odds[i % half]
        for i in range(len(coefficients))
    ]


def intt(
    primitive_root: FieldElement, values: list[FieldElement]
) -> list[FieldElement]:
    """
    Computes the Inverse Number Theoretic Transform (INTT) of a sequence of field elements.

    The INTT is the inverse operation of the NTT. It transforms a sequence of field elements
    from the point-value representation to the coefficient representation.

    This implementation uses a recursive divide-and-conquer approach (Cooley-Tukey algorithm).
    """
    if len(values) & (len(values) - 1) != 0:
        raise ValueError("cannot compute intt of non-power-of-two sequence")

    if len(values) == 1:
        return values

    field = values[0].field
    n = FieldElement(len(values), field)

    transformed_values = ntt(1 / primitive_root, values)
    return [tv / n for tv in transformed_values]


def fast_multiply(
    lhs: Polynomial, rhs: Polynomial, primitive_root: FieldElement, root_order: int
) -> Polynomial:
    """
    Computes the product of two polynomials using the Fast Fourier Transform (FFT).

    This function multiplies two polynomials efficiently by converting them to their
    point-value representations, performing a Hadamard product (element-wise multiplication),
    and then converting back to the coefficient representation.
    """
    if primitive_root**root_order != primitive_root.field.one:
        raise ValueError("supplied root does not have supplied order")
    if primitive_root ** (root_order // 2) == primitive_root.field.one:
        raise ValueError("supplied root is not primitive root of supplied order")
    if root_order <= lhs.degree + rhs.degree:
        raise ValueError("supplied root order is less than the degree of the product")

    if lhs == 0 or rhs == 0:
        return Polynomial([])

    degree = lhs.degree + rhs.degree
    # Optimization: use regular multiplication for small degrees
    if degree < 8:
        return lhs * rhs

    zero = lhs.coefficients[0].field.zero
    root = primitive_root
    order = root_order
    while degree < order // 2:
        root = root**2
        order = order // 2

    lhs_coefficients = lhs.coefficients + [zero] * (order - len(lhs.coefficients))
    rhs_coefficients = rhs.coefficients + [zero] * (order - len(rhs.coefficients))

    lhs_codeword = ntt(root, lhs_coefficients)
    rhs_codeword = ntt(root, rhs_coefficients)

    hadamard_product = [lhs * rhs for (lhs, rhs) in zip(lhs_codeword, rhs_codeword)]

    product_coefficients = intt(root, hadamard_product)
    return Polynomial(product_coefficients).trim()


def fast_zerofier(
    domain: list[FieldElement], primitive_root: FieldElement, root_order: int
) -> Polynomial:
    """
    Computes the zerofier polynomial for a given domain and primitive root.

    The zerofier polynomial is used to evaluate a polynomial at a given domain.
    It is defined as the product of the polynomials (x - domain[i]) for all i in the domain.
    """
    if primitive_root**root_order != primitive_root.field.one:
        raise ValueError("supplied root does not have supplied order")
    if primitive_root ** (root_order // 2) == primitive_root.field.one:
        raise ValueError("supplied root is not primitive root of supplied order")

    if len(domain) == 0:
        return Polynomial([])

    if len(domain) == 1:
        return Polynomial([-domain[0], primitive_root.field.one])

    half = len(domain) // 2

    left = fast_zerofier(domain[:half], primitive_root, root_order)
    right = fast_zerofier(domain[half:], primitive_root, root_order)
    return fast_multiply(left, right, primitive_root, root_order)


def fast_evaluate(
    polynomial: Polynomial,
    domain: list[FieldElement],
    primitive_root: FieldElement,
    root_order: int,
) -> list[FieldElement]:
    """
    Evaluates a polynomial at a given domain using the Fast Fourier Transform (FFT).

    This function evaluates a polynomial efficiently by converting it to its
    point-value representations, performing a Hadamard product (element-wise multiplication),
    and then converting back to the coefficient representation.
    """
    if primitive_root**root_order != primitive_root.field.one:
        raise ValueError("supplied root does not have supplied order")
    if primitive_root ** (root_order // 2) == primitive_root.field.one:
        raise ValueError("supplied root is not primitive root of supplied order")

    if len(domain) == 0:
        return []

    if len(domain) == 1:
        return [polynomial.evaluate(domain[0])]

    half = len(domain) // 2

    left_zerofier = fast_zerofier(domain[:half], primitive_root, root_order)
    right_zerofier = fast_zerofier(domain[half:], primitive_root, root_order)

    left = fast_evaluate(
        polynomial % left_zerofier, domain[:half], primitive_root, root_order
    )
    right = fast_evaluate(
        polynomial % right_zerofier, domain[half:], primitive_root, root_order
    )

    return left + right


def fast_interpolate(
    domain: list[FieldElement],
    values: list[FieldElement],
    primitive_root: FieldElement,
    root_order: int,
) -> Polynomial:
    if primitive_root**root_order != primitive_root.field.one:
        raise ValueError("supplied root does not have supplied order")
    if primitive_root ** (root_order // 2) == primitive_root.field.one:
        raise ValueError("supplied root is not primitive root of supplied order")
    if len(domain) != len(values):
        raise ValueError(f"interpolate needs {len(domain)=} == {len(values)=}")
    if len(set(domain)) != len(domain):
        raise ValueError("domain must contain unique elements")

    if len(domain) == 0:
        return Polynomial([])

    if len(domain) == 1:
        return Polynomial([values[0]])

    half = len(domain) // 2

    left_domain = domain[:half]
    right_domain = domain[half:]
    left_zerofier = fast_zerofier(left_domain, primitive_root, root_order)
    right_zerofier = fast_zerofier(right_domain, primitive_root, root_order)

    left_offset = fast_evaluate(right_zerofier, left_domain, primitive_root, root_order)
    right_offset = fast_evaluate(
        left_zerofier, right_domain, primitive_root, root_order
    )

    left_targets = [v / d for (v, d) in zip(values[:half], left_offset)]
    right_targets = [v / d for (v, d) in zip(values[half:], right_offset)]

    left_interpolant = fast_interpolate(
        left_domain, left_targets, primitive_root, root_order
    )
    right_interpolant = fast_interpolate(
        right_domain, right_targets, primitive_root, root_order
    )

    return left_interpolant * right_zerofier + right_interpolant * left_zerofier


def fast_coset_evaluate(
    polynomial: Polynomial,
    offset: FieldElement,
    generator: FieldElement,
    order: int,
) -> list[FieldElement]:
    """
    Evaluates a polynomial at a given offset and generator using the Number Theoretic Transform (NTT).
    """
    if generator**order != generator.field.one:
        raise ValueError("supplied generator does not have supplied order")
    if generator ** (order // 2) == generator.field.one:
        raise ValueError("supplied generator is not primitive root of supplied order")

    return ntt(
        generator,
        polynomial.scale(offset).coefficients
        + [offset.field.zero] * (order - len(polynomial.coefficients)),
    )


def fast_coset_divide(
    lhs: Polynomial,
    rhs: Polynomial,
    offset: FieldElement,
    primitive_root: FieldElement,
    root_order: int,
):
    if primitive_root**root_order != primitive_root.field.one:
        raise ValueError("supplied root does not have supplied order")
    if primitive_root ** (root_order // 2) == primitive_root.field.one:
        raise ValueError("supplied root is not primitive root of supplied order")
    if rhs == 0:
        raise ValueError("cannot divide by zero polynomial")

    if lhs == 0:
        return Polynomial([])

    if rhs.degree > lhs.degree:
        raise ValueError("cannot divide by polynomial of larger degree")

    root = primitive_root
    order = root_order
    degree = max(lhs.degree, rhs.degree)

    # Optimization: use regular division for small degrees
    if degree < 8:
        return lhs / rhs

    while degree < order // 2:
        root = root**2
        order = order // 2

    lhs_coefficients = lhs.scale(offset).coefficients + [offset.field.zero] * (
        order - len(lhs.coefficients)
    )
    rhs_coefficients = rhs.scale(offset).coefficients + [offset.field.zero] * (
        order - len(rhs.coefficients)
    )

    lhs_codeword = ntt(root, lhs_coefficients)
    rhs_codeword = ntt(root, rhs_coefficients)

    quotient_codeword = [lhs / rhs for (lhs, rhs) in zip(lhs_codeword, rhs_codeword)]
    scaled_quotient_coefficients = intt(root, quotient_codeword)
    scaled_quotient = Polynomial(
        scaled_quotient_coefficients[: (lhs.degree - rhs.degree + 1)]
    ).trim()

    return scaled_quotient.scale(1 / offset)
