from hashlib import blake2b
from typing import NewType, Protocol, Sequence, Union


class IntoBytes(Protocol):
    def __bytes__(self) -> bytes: ...


BytesConvertible = Union[IntoBytes, bytes]
Bytes64 = NewType("Bytes64", bytes)


class Merkle:
    @staticmethod
    def H(data: bytes) -> Bytes64:
        return Bytes64(blake2b(data).digest())

    @staticmethod
    def _commit(leafs: list[Bytes64]) -> Bytes64:
        """
        Commit a list of leafs.

        The leafs must be of even length.
        """
        if len(leafs) & (len(leafs) - 1) != 0:
            raise ValueError("length must be power of two")

        if len(leafs) == 1:
            return leafs[0]

        return Merkle.H(
            Merkle._commit(leafs[: len(leafs) // 2])
            + Merkle._commit(leafs[len(leafs) // 2 :])
        )

    @staticmethod
    def commit(data: Sequence[BytesConvertible]) -> Bytes64:
        """
        Commit a list of objects.

        The objects must implement the `__bytes__` method, and the initial leaf value is
        computed by hashing the bytes of the object.
        """
        return Merkle._commit([Merkle.H(bytes(da)) for da in data])

    @staticmethod
    def _open(index: int, leafs: list[Bytes64]) -> list[Bytes64]:
        """
        Open a leaf at a given index, i.e. return the inclusion proof for the leaf at `index`.
        """
        if len(leafs) & (len(leafs) - 1) != 0:
            raise ValueError("length must be power of two")
        if index < 0 or index >= len(leafs):
            raise ValueError("index not in range")

        if len(leafs) == 2:
            return [leafs[1 - index]]

        half = len(leafs) // 2
        if index < half:
            return Merkle._open(index, leafs[:half]) + [Merkle._commit(leafs[half:])]
        else:
            return Merkle._open(index - half, leafs[half:]) + [
                Merkle._commit(leafs[:half])
            ]

    @staticmethod
    def open(index: int, data: Sequence[BytesConvertible]) -> list[Bytes64]:
        return Merkle._open(index, [Merkle.H(bytes(da)) for da in data])

    @staticmethod
    def _verify(root: Bytes64, index: int, path: list[Bytes64], leaf: Bytes64) -> bool:
        """
        Verify that the given leaf is included in the Merkle tree with the given root
        at the given index, using the given inclusion proof.
        """
        if index < 0 or index >= (1 << len(path)):
            raise ValueError("index not in range")

        if index % 2 == 0:
            leaf = Merkle.H(leaf + path[0])
        else:
            leaf = Merkle.H(path[0] + leaf)

        if len(path) == 1:
            return root == leaf

        return Merkle._verify(root, index >> 1, path[1:], leaf)

    @staticmethod
    def verify(
        root: Bytes64, index: int, path: list[Bytes64], item: BytesConvertible
    ) -> bool:
        return Merkle._verify(root, index, path, Merkle.H(bytes(item)))
