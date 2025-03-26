import pickle as pickle
from hashlib import shake_256
from typing import Any, List, Optional


class ProofStream:
    def __init__(self, objects: Optional[List[Any]] = None):
        self.objects: List[Any] = objects or []
        self.read_index = len(self.objects)

    def push(self, obj: Any):
        self.objects += [obj]

    def pull(self) -> Any:
        if self.read_index >= len(self.objects):
            raise ValueError("ProofStream: cannot pull object; queue empty.")

        obj = self.objects[self.read_index]
        self.read_index += 1
        return obj

    def serialize(self, index: Optional[int] = None) -> bytes:
        if index is None:
            return pickle.dumps(self.objects)
        else:
            return pickle.dumps(self.objects[:index])

    def _fiat_shamir(self, num_bytes: int = 32, index: Optional[int] = None) -> bytes:
        """
        Generate a challenge from the proof stream.

        If index is provided, the challenge is generated from the part of the
        proof stream up to (but not including) the given index.
        """
        return shake_256(self.serialize(index)).digest(num_bytes)

    def prover_fiat_shamir(self, num_bytes: int = 32) -> bytes:
        """
        Fiat-Shamir transform for prover.

        The prover always uses the full proof stream to generate the challenge.
        """
        return self._fiat_shamir(num_bytes)

    def verifier_fiat_shamir(self, num_bytes: int = 32) -> bytes:
        """
        Fiat-Shamir transform for verifier.

        The verifier only uses the part of the proof stream that has been revealed
        so far.
        """
        return self._fiat_shamir(num_bytes, self.read_index)

    @staticmethod
    def deserialize(bb: bytes) -> "ProofStream":
        ps = ProofStream()
        ps.objects = pickle.loads(bb)
        return ps
