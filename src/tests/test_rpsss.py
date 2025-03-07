from stark_anatomy.fast_rpsss import FastRPSSS
from stark_anatomy.rpsss import RPSSS


class TestRPSSS:
    def test_rpsss(self):
        rpsss = RPSSS()

        sk, pk = rpsss.keygen()

        doc = b"Hello, world!"
        sig = rpsss.sign(sk, doc)

        assert rpsss.verify(pk, doc, sig)
        not_doc = b"Byebye."
        assert not rpsss.verify(pk, not_doc, sig)

    def test_fast_rpsss(self):
        rpsss = FastRPSSS()

        sk, pk = rpsss.keygen()

        doc = b"Hello, world!"
        sig = rpsss.sign(sk, doc)

        assert rpsss.verify(pk, doc, sig)
        not_doc = b"Byebye."
        assert not rpsss.verify(pk, not_doc, sig)
