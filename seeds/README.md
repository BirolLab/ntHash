# Linearity of hashing operations

An operation $f$ is linear over XOR if: $f(a \oplus b) = f(a) \oplus f(b)$. Both `rol` and `srol` are linear, as they move the same bits around for both $a$ and $b$.

Therefore, ntHash and ntHash2 are linear hash functions: $h(S) = f^{k-1}(s_0) \oplus f^{k-2}(s_1) \oplus \dots \oplus f(s_{k-1})$.

The hashes of two k-mers $X$ and $Y$ collide if $h(X) = h(Y)$. If we XOR both sides by $h(Y)$, we get $h(X) \oplus h(Y) = 0$. Since the hash functions are linear, $h(X) \oplus h(Y) = h(X \oplus Y)$, or the hash of the position-wise mismatches between the two k-mers.

# Searching for hash collisions through mismatches

Considering the 4-letter DNA alphabet, there are 16 possible character pairs for each position in two k-mers. Since XOR is symmetric ($a \oplus b = b \oplus a$) and self-cancelling ($a \oplus a = 0$), we only need to look at 6 unique non-zero values for the mismatches: $AC, AG, AT, CG, CT, GT$.

Currently, we assign each character to an independent random 64-bit value. The 6 differences can be controlled by exactly 3 independent basis vectors (degrees of freedom), using $A$ as the baseline:

- $d_1 = h(A) \oplus h(C)$
- $d_2 = h(A) \oplus h(G)$
- $d_3 = h(A) \oplus h(T)$

All other mismatches can be defined by combining these three. For example, $h(C) \oplus h(G) = d_1 \oplus d_2$.

These basis vectors represent the 6 possible differences in position 0. To extend to the $i^{\text{th}}$ position, we have to replace each base with the rotated seed ($f$ can be `rol` or `srol`):

- $d_{i,1} = f^{k-i-1}(h(A)) \oplus f^{k-i-1}(h(C)) = f^{k-i-1}(h(A) \oplus h(C)) = f^{k-i-1}(d_1)$
- Similarly: $d_{i,2} = f^{k-i-1}(d_2)$
- Similarly: $d_{i,3} = f^{k-i-1}(d_3)$

Therefore, each extra position adds 3 degrees of freedom to represent the total possible mismatches in the k-mer.

Since we are operating in the 64-bit space, our maximum capacity for basis vectors is 64. We can build every possible 64-bit value by XORing some of the following:

- $e_0 = 0000 \dots 0001$
- $e_1 = 0000 \dots 0010$
- $e_2 = 0000 \dots 0100$
- $\dots$
- $e_{63} = 1000 \dots 0000$

Because we need $3k$ vectors for the mismatches in the $k$ positions in a k-mer and the maximum capacity is 64, it is **impossible** to have a collision-free ntHash with 4 random seeds for $k \geq \lceil 64/3 \rceil = 22$.

# Theoretically-perfect hashing with a seed constraint

Ideally, we want to be collision free for all $k \leq 32$, since the total number of k-mers is less than $2^{64}$ in that range.

Ragnar proposes a nice solution [here](https://curiouscoding.nl/posts/nthash/#proving-perfection). If we take one of the seeds (for example, $T$) and make it the XOR of the other three, we get $h(T) = h(A) \oplus h(C) \oplus h(G)$, which now gives us:

- $d^*_1 = h(A) \oplus h(C)$
- $d^*_2 = h(A) \oplus h(G)$
- $d^*_3 = h(A) \oplus h(T) = h(A) \oplus h(A) \oplus h(C) \oplus h(G) = h(C) \oplus h(G) = d^*_1 \oplus d^*_2$

With this new seed design, we only need two basis vectors per position (rotations of $d^*_1$ and $d^*_2$), leading to a **mathematically possible** collision-free ntHash for $k \leq 64 / 2 = 32$.

# Designing seeds to fix systematic collisions

We now need to extend $d^*_1$ and $d^*_2$ to all positions in a k-mer ($k \leq 32$). Let's say we've chosen three seeds for $A, C, G$ and subsequently calculated the seed for $T$. This gives us $d^*_1$ and $d^*_2$ for position 0; we now have to verify that no rotation of these two basis vectors causes a collision in other positions.

For $k \leq 32$, we have up to 31 possible rotations of $d^*_1$ and $d^*_2$, resulting in $2^{64}$ combinations to check and make sure the XOR of no non-empty combination returns 0. This is intractable to loop over and validate with brute force.

In linear algebra, we say a set of vectors are linearly independent if no combination of them creates 0, and the only way to reach 0 is to combine none of them. Our problem now translates into finding a set of 64 vectors (32 for the rotations of $d^*_1$ and 32 for $d^*_2$) that are linearly independent. We do this by building a $64 \times 64$ matrix where the rows are the rotations of $d^*_1$ and $d^*_2$, and the columns are the bits. If the rank of this matrix is **exactly 64**, the vectors are linearly independent, and we have successfully found the seeds that make ntHash collision-free for $k \leq 32$.

This is implemented in `tests.cpp`, which checks the seeds in `internal.hpp`. The `generate_seeds.py` script generates and validates optimal seeds which can be piped into `make_internal_hpp.py` to update `internal.hpp`.
