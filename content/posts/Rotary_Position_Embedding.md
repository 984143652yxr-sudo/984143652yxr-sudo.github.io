---
{
  "slug": "rotary-position-embedding",
  "title": "Rotary position embeddings: from rotations to attention",
  "date": "2026-10-02",
  "date_label": "October 2, 2026",
  "status": "LLM study note",
  "reading_time": "12 min read",
  "summary": "A geometric and algebraic derivation of RoPE, followed by coordinate pairing, attention integration, and numerical checks."
}
---

## 1. The question: how does attention represent position?

Attention compares the query of one token with the keys of other tokens. Their dot products measure compatibility. **Rotary position embeddings (RoPE) rotate queries and keys according to their positions, so their dot product contains a relative-position rotation.** The construction was introduced in [RoFormer](https://arxiv.org/abs/2104.09864).

Consider one attention head. The input at position $m$ is a hidden-state vector $x_m$. Linear projections produce

$$
q_m=W_qx_m,\qquad k_m=W_kx_m,\qquad v_m=W_vx_m.
$$

| Symbol | Meaning |
|---|---|
| $m,n$ | Token positions, indexed from zero |
| $x_m\in\mathbb R^D$ | Hidden state at position $m$ |
| $q_m,k_m,v_m\in\mathbb R^d$ | Query, key, and value for one head |
| $d$ | Head dimension; even in the full-head construction below |
| $R_m$ | Rotation determined by position $m$ |
| $\widetilde q_m,\widetilde k_m$ | Position-rotated query and key |
| $\theta_i$ | Rotation frequency of coordinate pair $i$, in radians per token |

RoPE enters attention through the following sequence:

$$
\begin{aligned}
x_m&\longrightarrow(q_m,k_m,v_m),\\
(q_m,k_m)&\longrightarrow(R_mq_m,R_mk_m),\\
s_{mn}&=\frac{\widetilde q_m^\top\widetilde k_n}{\sqrt d},\\
a_{mn}&=\frac{\exp(s_{mn})}{\sum_{\ell\le m}\exp(s_{m\ell})},\quad n\le m,\\
o_m&=\sum_{n\le m}a_{mn}v_n.
\end{aligned}
$$

Here the last two lines describe **causal** attention: token $m$ can attend to positions up to $m$. Standard RoPE rotates queries and keys; values carry the information being averaged. Each score depends on token content and relative position. The softmax also depends on all keys available to that query.

## 2. One coordinate pair: a rotation in the plane

Begin with $d=2$. Define the rotation matrix

$$
R(\phi)=
\begin{pmatrix}
\cos\phi&-\sin\phi\\
\sin\phi&\cos\phi
\end{pmatrix}.
$$

For frequency $\theta$, position $m$ gives angle $m\theta$:

$$
\widetilde q_m=R(m\theta)q_m,
\qquad
\widetilde k_n=R(n\theta)k_n.
$$

The key identity follows from $R(\phi)^\top=R(-\phi)$ and addition of rotation angles:

$$
\boxed{
\widetilde q_m^\top\widetilde k_n
=q_m^\top R(m\theta)^\top R(n\theta)k_n
=q_m^\top R((n-m)\theta)k_n.
}
$$

Each token receives its own position rotation. The comparison between two tokens contains the **signed relative offset $n-m$**. For fixed unrotated queries and keys, shifting both positions by the same amount leaves this score unchanged.

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/rope/relative-rotation.svg"><img src="{{base}}/assets/rope/relative-rotation.svg" alt="Query and key unit vectors at positions 1 and 3, then 4 and 6. Both pairs have a 60 degree relative angle and dot product 0.5."></a>
<figcaption>A two-dimensional example. The common shift changes both absolute angles while preserving their difference. Click the figure to enlarge.</figcaption>
</figure>

For general $q=(q_1,q_2)$ and $k=(k_1,k_2)$, writing $\delta=(n-m)\theta$ makes the content-position interaction explicit:

$$
q^\top R(\delta)k
=(q_1k_1+q_2k_2)\cos\delta
+(q_2k_1-q_1k_2)\sin\delta.
$$

The two coefficients depend on content. The cosine and sine depend on position. In particular, the sign of the relative offset can matter.

### The same derivation with complex numbers

Identify the real vector $(q_1,q_2)$ with $z_q=q_1+\mathrm i q_2$. Multiplication by $e^{\mathrm i m\theta}$ performs the rotation. The real-vector dot product is the real part of a conjugate product:

$$
\widetilde q_m^\top\widetilde k_n
=\operatorname{Re}\!\left[
(z_qe^{\mathrm i m\theta})\overline{(z_ke^{\mathrm i n\theta})}
\right]
=\operatorname{Re}\!\left[z_q\overline{z_k}e^{\mathrm i(m-n)\theta}\right].
$$

The complex formula uses $m-n$ because the key is conjugated. It gives exactly the same real scalar as the matrix formula above.

## 3. A full attention head: many rotation frequencies

For even $d$, divide the coordinates into $d/2$ pairs. In the **adjacent-pair convention**, these are $(0,1),(2,3),\ldots,(d-2,d-1)$. Define

$$
\theta_i=b^{-2i/d},\qquad i=0,\ldots,d/2-1,
$$

where $b=10{,}000$ is the original illustrative base. Model configurations may use other bases or modify the frequency schedule.

The complete rotation is block diagonal:

$$
R_m=\operatorname{diag}\!\left(
R(m\theta_0),R(m\theta_1),\ldots,R(m\theta_{d/2-1})
\right).
$$

All the two-dimensional identities apply block by block:

$$
\widetilde q_m^\top\widetilde k_n=q_m^\top R_{n-m}k_n,
\qquad
\|R_mq_m\|_2=\|q_m\|_2.
$$

The score sums contributions from all coordinate pairs. Early pairs rotate rapidly; later pairs change more slowly with token position. The period of pair $i$ is $2\pi/\theta_i$ tokens. For example, $d=8$ and $b=10{,}000$ give frequencies $(1,0.1,0.01,0.001)$.

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/rope/frequencies.svg"><img src="{{base}}/assets/rope/frequencies.svg" alt="Multiple rotation frequencies and an oscillating dot product showing that larger distances do not guarantee smaller scores."></a>
<figcaption>Top: three frequency components. Bottom: a fixed unit query and key at frequency π/6; the score returns to 1 after 12 tokens.</figcaption>
</figure>

### What can be concluded about distance?

For the simple choice $q=k=(1,0)$, the dot product is $\cos((n-m)\theta)$. It oscillates. Thus an individual RoPE score can increase at a larger distance. The paper's discussion of long-term decay concerns the aggregate multi-frequency construction and its analysis; it does not imply monotonic decay for every query-key pair.

RoPE can be evaluated at positions beyond those seen during training. Whether a trained model uses those positions reliably is a separate empirical question. Frequency scaling, the training context length, and the model's learned weights all affect long-context behavior.

The shared-shift identity also holds with **fixed input queries and keys**. Hidden states in deeper layers already reflect context, attention masks, and previous computations; the identity alone does not establish translation invariance of the whole language model.

## 4. Implementation: rotate pairs without building a matrix

For one pair $(a,b)$ and angle $\phi$, compute

$$
(a,b)\longmapsto
(a\cos\phi-b\sin\phi,\ a\sin\phi+b\cos\phi).
$$

The following NumPy reference accepts shape `(..., tokens, head_dim)`. Leading dimensions can represent batches and heads. An explicit `positions` argument also supports decoding with a key-value cache.

```python
import numpy as np


def rope(x, positions, base=10_000.0):
    d = x.shape[-1]
    if d % 2:
        raise ValueError("head_dim must be even")
    positions = np.asarray(positions)
    if positions.shape != (x.shape[-2],):
        raise ValueError("One position is required per token")

    freq = base ** (-np.arange(0, d, 2, dtype=float) / d)
    angles = positions[:, None] * freq[None, :]
    pairs = x.reshape(*x.shape[:-1], d // 2, 2)
    a, b = pairs[..., 0], pairs[..., 1]
    co, si = np.cos(angles), np.sin(angles)
    return np.stack(
        (a * co - b * si, a * si + b * co), axis=-1
    ).reshape(x.shape)
```

This is an algebraic reference in floating-point NumPy. A training implementation additionally manages device placement, tensor dtype, and cached sine/cosine tables.

### Why do some implementations split the vector in half?

Another convention stores the first coordinate of every pair together, followed by the second coordinate of every pair. For $d=4$:

| Layout | Stored vector | Pairs being rotated |
|---|---|---|
| Adjacent | $(a,b,c,d)$ | $(a,b)$ and $(c,d)$ |
| Split-half | $(a,c,b,d)$ | $(a,b)$ and $(c,d)$ |

With split-half storage $x=[x_1,x_2]$, the expression

$$
\operatorname{rotateHalf}(x)=[-x_2,x_1],
\qquad
\widetilde x=x\odot\cos\Phi+
\operatorname{rotateHalf}(x)\odot\sin\Phi
$$

uses the angle vector $\Phi=(m\theta_0,\ldots,m\theta_{d/2-1},m\theta_0,\ldots,m\theta_{d/2-1})$. This is the convention used by the split-half code in [Raschka's Llama conversion notebook](https://github.com/rasbt/LLMs-from-scratch/blob/main/ch05/07_gpt_to_llama/converting-gpt-to-llama2.ipynb).

The two conventions are related by a coordinate permutation. Queries, keys, and their projection weights must use a consistent layout. When loading pretrained weights, the implementation must match the checkpoint's convention. [Meta's original Llama implementation](https://github.com/meta-llama/llama/blob/main/llama/model.py) represents adjacent pairs as complex numbers.

## 5. Where RoPE sits in causal attention

After linear projections and reshaping, suppose `q`, `k`, and `v` each have shape `(batch, heads, T, d)`. For a complete causal sequence, the reference calculation is:

```python
positions = np.arange(q.shape[-2])
qr = rope(q, positions)
kr = rope(k, positions)

scores = qr @ kr.swapaxes(-1, -2) / np.sqrt(q.shape[-1])
T = q.shape[-2]
future = np.triu(np.ones((T, T), dtype=bool), k=1)
scores = np.where(future, -np.inf, scores)

# Stable softmax over keys, separately for each query.
weights = np.exp(scores - scores.max(axis=-1, keepdims=True))
weights /= weights.sum(axis=-1, keepdims=True)
output = weights @ v
```

The head outputs are subsequently concatenated and passed through the output projection. This example assumes equal query/key sequence lengths and the same number of query and key heads.

### Positions during cached decoding

After processing positions $0,\ldots,t-1$, the next token has position $t$ even when its tensor has sequence length one. Rotate its new query and key using `positions=np.array([t])`. Compare that query against cached rotated keys and the new rotated key. Previously rotated keys retain their original positions. The corresponding cached values remain unrotated.

A useful distinction is **tensor index within the current chunk** versus **position within the sequence**. The frequency calculation needs the latter. Cache layouts, padding, and chunked decoding require consistent position IDs and attention masks.

## 6. Numerical checks

The [downloadable reference script]({{base}}/assets/rope/rope_demo.py) includes the implementation and executable checks. Run it with Python and NumPy:

```bash
python rope_demo.py
```

The checks cover:

| Property | Check |
|---|---|
| Norm preservation | $\|R_mq\|_2=\|q\|_2$ |
| Relative-position identity | $(R_mq)^\top(R_nk)=q^\top R_{n-m}k$ |
| Common shift | Scores agree for positions $(m,n)$ and $(m+s,n+s)$ |
| Cached token | Rotating a single token at its original position matches the full-sequence result |
| Pairing convention | Permuting adjacent-pair output agrees with rotation in split-half storage |

These checks pass for the supplied deterministic random example. They verify the reference algebra and indexing; trained-model quality requires evaluation on the intended tasks and context lengths.

## References

1. Su et al. [RoFormer: Enhanced Transformer with Rotary Position Embedding](https://arxiv.org/abs/2104.09864). Original formulation, complex-number derivation, and multi-frequency construction.
2. Meta. [Llama model implementation](https://github.com/meta-llama/llama/blob/main/llama/model.py). Complex-valued rotation and position handling during cached inference.
3. Sebastian Raschka. [Converting GPT to Llama 2](https://github.com/rasbt/LLMs-from-scratch/blob/main/ch05/07_gpt_to_llama/converting-gpt-to-llama2.ipynb). A code-oriented reference for attention and the split-half implementation.
