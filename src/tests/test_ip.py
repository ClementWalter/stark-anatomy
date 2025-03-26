from hypothesis import given
from hypothesis import strategies as st

from stark_anatomy.ip import ProofStream


class TestIp:
    @given(proof_stream=...)
    def test_serialize(self, proof_stream: ProofStream):
        assert (
            ProofStream.deserialize(proof_stream.serialize()).objects
            == proof_stream.objects
        )

    @given(proof_stream=...)
    def test_prover_fiat_shamir_should_return_last_fiat_shamir(
        self, proof_stream: ProofStream
    ):
        assert proof_stream.prover_fiat_shamir() == proof_stream._fiat_shamir(
            index=len(proof_stream.objects)
        )

    @given(proof_stream=..., data=st.data())
    def test_verifier_fiat_shamir_should_return_hash_at_index(
        self, proof_stream: ProofStream, data: st.DataObject
    ):
        index = data.draw(st.integers(min_value=0, max_value=len(proof_stream.objects)))
        proof_stream.read_index = index
        assert proof_stream.verifier_fiat_shamir() == proof_stream._fiat_shamir(
            index=index
        )

    @given(num_bytes=st.integers(min_value=1, max_value=32))
    def test_fiat_shamir_should_return_num_bytes(self, num_bytes: int):
        proof_stream = ProofStream()
        assert len(proof_stream.prover_fiat_shamir(num_bytes)) == num_bytes
