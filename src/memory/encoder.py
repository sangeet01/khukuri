"""
encoder.py
----------
Thin wrapper around LimitNumen's character n-gram hash encoder.
Swappable: swap this file to use any other encoder -- interface stays the same.

Interface contract:
    encoder = NumenEncoder(dim=4096)
    vec = encoder.encode("some text")          # -> np.ndarray (dim,) float32, L2-normed
    vecs = encoder.encode_batch(["a", "b"])    # -> np.ndarray (N, dim) float32, L2-normed
"""

import numpy as np


class NumenEncoder:
    """
    Training-free character n-gram hash encoder.
    Polynomial rolling hash (deterministic, no CRC32 Python loop).
    Drop-in: replace __init__ and _encode_one to swap backend.
    """

    def __init__(self, dim: int = 4096, ns: tuple = (3, 4, 5), p: int = 31):
        self.dim = dim
        self.ns  = ns
        # precompute power tables once
        self._pows = {
            n: np.array([pow(p, i, dim * 65537) for i in range(n)], dtype=np.int64)
            for n in ns
        }

    # -- public ----------------------------------------------------------------

    def encode(self, text: str) -> np.ndarray:
        vec = self._encode_one(text)
        return self._normalise(vec)

    def encode_batch(self, texts: list[str]) -> np.ndarray:
        out = np.zeros((len(texts), self.dim), dtype=np.float32)
        for i, t in enumerate(texts):
            out[i] = self._encode_one(t)
        # batch normalise
        norms = np.linalg.norm(out, axis=1, keepdims=True)
        norms = np.where(norms == 0, 1.0, norms)
        return (out / norms).astype(np.float32)

    # -- internal --------------------------------------------------------------

    def _encode_one(self, text: str) -> np.ndarray:
        b = np.frombuffer(text.lower().encode(), dtype=np.uint8).astype(np.int64)
        vec = np.zeros(self.dim, dtype=np.float32)
        for n in self.ns:
            if len(b) < n:
                continue
            wins = np.lib.stride_tricks.as_strided(
                b,
                shape=(len(b) - n + 1, n),
                strides=(b.strides[0], b.strides[0]),
            )
            np.add.at(vec, (wins @ self._pows[n]) % self.dim, 1.0)
        return vec

    @staticmethod
    def _normalise(vec: np.ndarray) -> np.ndarray:
        norm = np.linalg.norm(vec)
        return (vec / norm).astype(np.float32) if norm > 0 else vec.astype(np.float32)
