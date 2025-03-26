import pytest
from hypothesis import strategies as st
from hypothesis.strategies import register_type_strategy

from stark_anatomy.algebra import Field, FieldElement
from stark_anatomy.ip import ProofStream
from stark_anatomy.multivariate import MPolynomial
from stark_anatomy.univariate import Polynomial


@pytest.fixture(scope="session")
def field():
    return Field.default()


register_type_strategy(
    FieldElement,
    st.integers(min_value=0, max_value=Field.PRIME - 1).map(
        lambda value: FieldElement(value, Field.default())
    ),
)

register_type_strategy(
    Polynomial,
    st.lists(st.from_type(FieldElement), min_size=0, max_size=32).map(
        lambda coefficients: Polynomial(coefficients).trim()
    ),
)

register_type_strategy(
    MPolynomial,
    st.integers(min_value=1, max_value=5).flatmap(
        lambda num_variables: st.dictionaries(
            st.tuples(
                *[st.integers(min_value=0, max_value=32) for _ in range(num_variables)]
            ),
            st.from_type(FieldElement),
            min_size=0,
            max_size=32,
        ).map(lambda coefficients: MPolynomial(coefficients).trim())
    ),
)

register_type_strategy(
    ProofStream,
    st.lists(st.binary()).map(lambda objects: ProofStream(objects)),
)
