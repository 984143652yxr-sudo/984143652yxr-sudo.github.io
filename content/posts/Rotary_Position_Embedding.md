##RoPE (Rotary Position Embedding) is a way to use absolute position rotation matrix to capture the relationship between different words in a sentence.
It is a position embedding method popular in Llama and GLM based Language model.

In standard self-attention model, for a given position m and the embedding of the word X_m, and another word x_n at position n, the query and key
could be obtained through linear transformation 
q_m = f_q(x_m,m),k_m = f_k(x_m,m),q_n = f_q(x_n,n), k_n = f_k(x_n,n)
The aim of RoPe is to find a transformation f_1 and f_k, so that the inner product of q_m and k_n, which captures the attention score, depend purely on
x_m, x_n and their absolute distance (m-n). This means we want to find a function g, sucg that 
<f_q(x_m,m),f_k(x_n,n)> = f(x_m,x_n,m-n)

ReFormer utilize the complex number in 2d PLANE and Eular formula to achieve the goal.
Suppose we have a word embedding of dimension 2, then the query in the complex plane would be q=[q^{(1)},q^{(2)}] with q = q^{(1)}+iq^{(2)}.

To incooperate the position m, we could multiply a transformation factor e^{im} to capture the position at m, which in reality would be 
f_q(x_m,m) = q_m e^{im\theta} and f_k(x_n,n) = k_n e^{intheta}

When we do inner product in complex space for the two rotated vector 

<f_q(x_m,m),f_k(x_n,n)> = Re((q_m e^{imtheta})(k_n e^{intheta})^{*}) = Re(q_m k_n^{*} e^{i(m-n)theta})
we could see the inner product contains the e^{i(m-n)\theta} capture the relative position (m-n), this means only doing transformation on their own position,
after the dot product in complex plain, we get the relative position of two vectors.

# What is it like in d-Dimensions in real case

In real world, the attention head usually has d = 64/128, the real implementation of RoPE is to selerate the high dimensional vector by two, 
then we get  d/2 subspace each with dimension 2, and we did transformation on each space,
By defining the diagonal transformation matrix R^d_{\Theta,m}:
R^d_{\Theta,m} = cos(m\theta_1) -sin(m\theta_1) 0. 0   
                 sin(m\theta_1)  cos(m\theta_1). 0. 0
                 0.              0              cos(m\theta_2) -sin(m\theta_2)
                 0.              0              sin(m\theta_2)  cos(m\theta_2)
Then q_m = R_{\Theta,m}^d(W_q x_m), as we see each subspace has different transformation degree theta_i = 10000^{-2(i-1)/d} for i \in [1,2,..,d/2]

As we see with |m-n| increase, the self-attention score <q_m,k_n> will decrease, which relates nartural human language pattern. Words with higher distance, usually has less relationship. 
DUring the computation, the position is computed through absolute information, while in attention computation, the relative position is preserved by the mathematical computation.

In reality, we did rotated = [-x2,x1]*sin + [x1,x2]*cos to get the final score to save computation without constructing the large matrix. 

def precompute_rope_params(head_dim, theta_base=10_000, context_length=4096):
    assert head_dim % 2 == 0, "Embedding dimension must be even"

    # Compute the inverse frequencies
    inv_freq = 1.0 / (theta_base ** (torch.arange(0, head_dim, 2)[: (head_dim // 2)].float() / head_dim))

    # Generate position indices
    positions = torch.arange(context_length)

    # Compute the angles
    angles = positions.unsqueeze(1) * inv_freq.unsqueeze(0)  # Shape: (context_length, head_dim // 2)

    # Expand angles to match the head_dim
    angles = torch.cat([angles, angles], dim=1)  # Shape: (context_length, head_dim)

    # Precompute sine and cosine
    cos = torch.cos(angles)
    sin = torch.sin(angles)

    return cos, sin

def compute_rope(x, cos, sin):
    # x: (batch_size, num_heads, seq_len, head_dim)
    batch_size, num_heads, seq_len, head_dim = x.shape
    assert head_dim % 2 == 0, "Head dimension must be even"

    # Split x into first half and second half
    x1 = x[..., : head_dim // 2]  # First half
    x2 = x[..., head_dim // 2 :]  # Second half

    # Adjust sin and cos shapes
    cos = cos[:seq_len, :].unsqueeze(0).unsqueeze(0)  # Shape: (1, 1, seq_len, head_dim)
    sin = sin[:seq_len, :].unsqueeze(0).unsqueeze(0)

    # Apply the rotary transformation
    rotated = torch.cat((-x2, x1), dim=-1)
    x_rotated = (x * cos) + (rotated * sin)

    return x_rotated.to(dtype=x.dtype)




class MultiHeadAttention(nn.Module):
    def __init__(self, d_in, d_out, context_length, num_heads, dtype=None):  # ,dropout, num_heads, qkv_bias=False):
        super().__init__()
        assert d_out % num_heads == 0, "d_out must be divisible by n_heads"

        self.d_out = d_out
        self.num_heads = num_heads
        self.head_dim = d_out // num_heads  # Reduce the projection dim to match desired output dim

        ################################### NEW ###################################
        # Set bias=False and dtype=dtype for all linear layers below
        ###########################################################################
        self.W_query = nn.Linear(d_in, d_out, bias=False, dtype=dtype)
        self.W_key = nn.Linear(d_in, d_out, bias=False, dtype=dtype)
        self.W_value = nn.Linear(d_in, d_out, bias=False, dtype=dtype)
        self.out_proj = nn.Linear(d_out, d_out, bias=False, dtype=dtype)  # Linear layer to combine head outputs
        # self.dropout = nn.Dropout(dropout)
        self.register_buffer("mask", torch.triu(torch.ones(context_length, context_length), diagonal=1))

        ################################### NEW ###################################
        cos, sin = precompute_rope_params(head_dim=self.head_dim, context_length=context_length)
        self.register_buffer("cos", cos)
        self.register_buffer("sin", sin)
        ###########################################################################


    def forward(self, x):

        b, num_tokens, d_in = x.shape

        keys = self.W_key(x)  # Shape: (b, num_tokens, d_out)
        queries = self.W_query(x)
        values = self.W_value(x)

        # We implicitly split the matrix by adding a `num_heads` dimension
        # Unroll last dim: (b, num_tokens, d_out) -> (b, num_tokens, num_heads, head_dim)
        keys = keys.view(b, num_tokens, self.num_heads, self.head_dim)
        values = values.view(b, num_tokens, self.num_heads, self.head_dim)
        queries = queries.view(b, num_tokens, self.num_heads, self.head_dim)

        # Transpose: (b, num_tokens, num_heads, head_dim) -> (b, num_heads, num_tokens, head_dim)
        keys = keys.transpose(1, 2)
        queries = queries.transpose(1, 2)
        values = values.transpose(1, 2)

        ################################### NEW ###################################
        keys = compute_rope(keys, self.cos, self.sin)
        queries = compute_rope(queries, self.cos, self.sin)
        ###########################################################################

        # Compute scaled dot-product attention (aka self-attention) with a causal mask
        attn_scores = queries @ keys.transpose(2, 3)  # Dot product for each head

        # Original mask truncated to the number of tokens and converted to boolean
        mask_bool = self.mask.bool()[:num_tokens, :num_tokens]

        # Use the mask to fill attention scores
        attn_scores.masked_fill_(mask_bool, -torch.inf)

        attn_weights = torch.softmax(attn_scores / keys.shape[-1]**0.5, dim=-1)
        # attn_weights = self.dropout(attn_weights)

        # Shape: (b, num_tokens, num_heads, head_dim)
        context_vec = (attn_weights @ values).transpose(1, 2)

        # Combine heads, where self.d_out = self.num_heads * self.head_dim
        context_vec = context_vec.reshape(b, num_tokens, self.d_out)
        context_vec = self.out_proj(context_vec)  # optional projection

        return context_vec
