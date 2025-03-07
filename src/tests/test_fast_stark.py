import os

from stark_anatomy.algebra import Field
from stark_anatomy.fast_stark import FastStark
from stark_anatomy.rescue_prime import RescuePrime


class TestFastStark:
    def test_fast_stark(self, field: Field):
        expansion_factor = 4
        num_colinearity_checks = 2
        security_level = 2

        rp = RescuePrime()
        output_element = field.sample(bytes(b"0xdeadbeef"))

        for _ in range(0, 20):
            input_element = output_element
            output_element = rp.hash(input_element)
            num_cycles = rp.N + 1
            state_width = rp.m

            stark = FastStark(
                field,
                expansion_factor,
                num_colinearity_checks,
                security_level,
                state_width,
                num_cycles,
            )
            (
                transition_zerofier,
                transition_zerofier_codeword,
                transition_zerofier_root,
            ) = stark.preprocess()

            # prove
            trace = rp.trace(input_element)
            air = rp.transition_constraints(stark.omicron)
            boundary = rp.boundary_constraints(output_element)
            proof = stark.prove(
                trace, air, boundary, transition_zerofier, transition_zerofier_codeword
            )

            # verify
            verdict = stark.verify(proof, air, boundary, transition_zerofier_root)

            assert verdict, "valid stark proof fails to verify"

            # verify false claim
            output_element_ = output_element + field.one
            boundary_ = rp.boundary_constraints(output_element_)
            verdict = stark.verify(proof, air, boundary_, transition_zerofier_root)

            assert not verdict, "invalid stark proof verifies"

            # prove with false witness
            cycle = 1 + (int(os.urandom(1)[0]) % len(trace) - 1)
            register = int(os.urandom(1)[0]) % state_width
            error = field.sample(os.urandom(17))

            trace[cycle][register] = trace[cycle][register] + error

            proof = stark.prove(
                trace, air, boundary, transition_zerofier, transition_zerofier_codeword
            )

            verdict = stark.verify(proof, air, boundary, transition_zerofier_root)
            assert not verdict, "STARK produced from false witness verifies :("
