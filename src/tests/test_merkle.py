from collections import namedtuple

from hypothesis import assume, given
from hypothesis import strategies as st

from stark_anatomy.merkle import Bytes64, Merkle

MerkleTestData = namedtuple("MerkleTestData", ["data_array", "leafs", "root", "index"])


@st.composite
def data_array(draw: st.DrawFn) -> MerkleTestData:
    log_n = draw(st.integers(min_value=1, max_value=6))
    n = 1 << log_n
    item = st.binary(min_size=0, max_size=64)
    data_array = draw(st.sets(item, min_size=n, max_size=n))
    index = draw(st.integers(min_value=0, max_value=n - 1))
    leafs = [Merkle.H(item) for item in data_array]
    root = Merkle._commit(leafs)
    return MerkleTestData(list(data_array), leafs, root, index)


class TestMerkle:
    @given(data=data_array())
    def test_commit(self, data: MerkleTestData):
        assert Merkle.commit(data.data_array) == Merkle._commit(data.leafs)

    @given(data=data_array())
    def test_open(self, data: MerkleTestData):
        assert Merkle.open(data.index, data.data_array) == Merkle._open(
            data.index, data.leafs
        )

    @given(data=data_array())
    def test_verify(self, data: MerkleTestData):
        assert Merkle.verify(
            data.root,
            data.index,
            Merkle.open(data.index, data.data_array),
            data.data_array[data.index],
        ) == Merkle._verify(
            data.root,
            data.index,
            Merkle.open(data.index, data.data_array),
            data.leafs[data.index],
        )

    @given(data=data_array())
    def test_open_verify_consistency(self, data: MerkleTestData):
        assert Merkle._verify(
            data.root,
            data.index,
            Merkle._open(data.index, data.leafs),
            data.leafs[data.index],
        )

    @given(data=data_array(), leaf=st.binary(min_size=64, max_size=64))
    def test_verify_should_fail_wrong_leaf(self, data: MerkleTestData, leaf: Bytes64):
        assume(leaf != data.leafs[data.index])
        assert not Merkle._verify(
            data.root,
            data.index,
            Merkle._open(data.index, data.leafs),
            leaf,
        )

    @given(data=data_array(), index=...)
    def test_verify_should_fail_wrong_index(self, data: MerkleTestData, index: int):
        index = index % len(data.data_array)
        assume(index != data.index)
        assert not Merkle._verify(
            data.root,
            index,
            Merkle._open(index, data.leafs),
            data.leafs[data.index],
        )

    @given(data=data_array(), root=st.binary(min_size=64, max_size=64))
    def test_verify_should_fail_wrong_root(self, data: MerkleTestData, root: Bytes64):
        assume(root != data.root)
        assert not Merkle._verify(
            root,
            data.index,
            Merkle._open(data.index, data.leafs),
            data.leafs[data.index],
        )

    @given(data=data_array(), index=..., node=st.binary(min_size=64, max_size=64))
    def test_verify_should_fail_wrong_path(
        self, data: MerkleTestData, index: int, node: Bytes64
    ):
        path = Merkle._open(data.index, data.leafs)
        path[index % len(path)] = node
        assert not Merkle._verify(
            data.root,
            data.index,
            path,
            data.leafs[data.index],
        )
