from stark_anatomy.ip import ProofStream


class TestIp:
    def test_serialize(self):
        proof1 = ProofStream()
        proof1.push(1)
        proof1.push({1: "1"})
        proof1.push([1])
        proof1.push(2)

        serialized = proof1.serialize()
        proof2 = ProofStream.deserialize(serialized)

        assert proof1.pull() == proof2.pull()
        assert proof1.pull() == proof2.pull()
        assert proof1.pull() == proof2.pull()
        assert proof1.pull() == 2
        assert proof2.pull() == 2
        assert proof1.prover_fiat_shamir() == proof2.prover_fiat_shamir()
