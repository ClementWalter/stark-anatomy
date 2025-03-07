import os

from stark_anatomy.merkle import Merkle


class TestMerkle:
    def test_merkle(self):
        n = 64
        leafs = [os.urandom(int(os.urandom(1)[0])) for i in range(n)]
        root = Merkle.commit_(leafs)

        # opening any leaf should work
        for i in range(n):
            path = Merkle.open_(i, leafs)
            assert Merkle.verify_(root, i, path, leafs[i])

        # opening non-leafs should not work
        for i in range(n):
            path = Merkle.open_(i, leafs)
            assert not Merkle.verify_(root, i, path, os.urandom(51))

        # opening wrong leafs should not work
        for i in range(n):
            path = Merkle.open_(i, leafs)
            j = (i + 1 + (int(os.urandom(1)[0] % (n - 1)))) % n
            assert not Merkle.verify_(root, i, path, leafs[j])

        # opening leafs with the wrong index should not work
        for i in range(n):
            path = Merkle.open_(i, leafs)
            j = (i + 1 + (int(os.urandom(1)[0] % (n - 1)))) % n
            assert not Merkle.verify_(root, j, path, leafs[i])

        # opening leafs to a false root should not work
        for i in range(n):
            path = Merkle.open_(i, leafs)
            assert not Merkle.verify_(os.urandom(32), i, path, leafs[i])

        # opening leafs with even one falsehood in the path should not work
        for i in range(n):
            path = Merkle.open_(i, leafs)
            for j in range(len(path)):
                fake_path = path[0:j] + [os.urandom(32)] + path[j + 1 :]
                assert not Merkle.verify_(root, i, fake_path, leafs[i])

        # opening leafs to a different root should not work
        fake_root = Merkle.commit_([os.urandom(32) for i in range(n)])
        for i in range(n):
            path = Merkle.open_(i, leafs)
            assert not Merkle.verify_(fake_root, i, path, leafs[i])
