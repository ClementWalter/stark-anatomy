import os
import random

import pytest

from stark_anatomy.algebra import Field
from stark_anatomy.rescue_prime import RescuePrime
from stark_anatomy.stark import Stark


class TestStark:
    def test_stark(self, field: Field):
        expansion_factor = 4
        num_colinearity_checks = 2
        security_level = 2

        rp = RescuePrime()
        output_element = field.sample(bytes(b"0xdeadbeef"))

        for _ in range(20):
            input_element = output_element
            output_element = rp.hash(input_element)
            num_cycles = rp.N + 1
            state_width = rp.m

            stark = Stark(
                field,
                expansion_factor,
                num_colinearity_checks,
                security_level,
                state_width,
                num_cycles,
            )

            # prove
            trace = rp.trace(input_element)
            air = rp.transition_constraints(stark.omicron)
            boundary = rp.boundary_constraints(output_element)
            proof = stark.prove(trace, air, boundary)

            # verify
            assert stark.verify(proof, air, boundary)

            output_element_ = output_element + field.one
            boundary_ = rp.boundary_constraints(output_element_)
            assert not stark.verify(proof, air, boundary_)

        # verify with false witness
        cycle = random.randint(0, len(trace) - 1)
        register = random.randint(0, state_width - 1)
        error = field.sample(os.urandom(17))

        trace[cycle][register] = trace[cycle][register] + error

        with pytest.raises(
            AssertionError,
            match="cannot perform polynomial division because remainder is not zero",
        ):
            stark.prove(trace, air, boundary)
