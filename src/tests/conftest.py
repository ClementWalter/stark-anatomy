import pytest

from stark_anatomy.algebra import Field


@pytest.fixture(scope="session")
def field():
    return Field.default()
