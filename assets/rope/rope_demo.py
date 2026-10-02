"""NumPy reference implementation of RoPE; run this file to verify its identities."""
import numpy as np


def rope(x, positions, base=10_000.0):
    """Rotate (..., tokens, head_dim) using adjacent coordinate pairs."""
    d = x.shape[-1]
    if d % 2:
        raise ValueError('head_dim must be even')
    positions = np.asarray(positions)
    if positions.shape != (x.shape[-2],):
        raise ValueError('One position is required per token')
    freq = base ** (-np.arange(0, d, 2, dtype=float) / d)
    angles = positions[:, None] * freq[None, :]
    pairs = x.reshape(*x.shape[:-1], d // 2, 2)
    a, b = pairs[..., 0], pairs[..., 1]
    co, si = np.cos(angles), np.sin(angles)
    return np.stack((a * co - b * si, a * si + b * co), axis=-1).reshape(x.shape)


def verify():
    rng = np.random.default_rng(7)
    q, k = rng.normal(size=(2, 2, 3, 5, 8))
    pos = np.arange(5)
    qr, kr = rope(q, pos), rope(k, pos)
    np.testing.assert_allclose(np.linalg.norm(qr, axis=-1), np.linalg.norm(q, axis=-1))
    scores = qr @ kr.swapaxes(-1, -2)
    shifted = rope(q, pos + 11) @ rope(k, pos + 11).swapaxes(-1, -2)
    np.testing.assert_allclose(scores, shifted, atol=1e-12)
    for m in range(5):
        for n in range(5):
            rhs = np.sum(q[..., m:m+1, :] * rope(k[..., n:n+1, :], np.array([n-m])), axis=-1)
            np.testing.assert_allclose(scores[..., m, n], rhs[..., 0], atol=1e-12)
    # Cached decoding must use the token's original absolute position.
    np.testing.assert_allclose(rope(q[..., 4:5, :], np.array([4])), qr[..., 4:5, :])
    # Split-half storage is the adjacent layout permuted to [all first, all second].
    d = q.shape[-1]
    permutation = np.r_[np.arange(0, d, 2), np.arange(1, d, 2)]
    split = q[..., permutation]
    a, b = np.split(split, 2, axis=-1)
    angle = pos[:, None] * 10_000.0 ** (-np.arange(0, d, 2) / d)
    split_rotated = np.concatenate((a*np.cos(angle)-b*np.sin(angle), a*np.sin(angle)+b*np.cos(angle)), axis=-1)
    np.testing.assert_allclose(split_rotated, qr[..., permutation], atol=1e-12)
    print('Passed: norms, relative-position identity, shared shift, cache positions, pairing permutation.')


if __name__ == '__main__':
    verify()
