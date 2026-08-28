import functools
import itertools
import operator
import os
import random
import typing

import pynthash


def generate(num_chars: int = 4) -> tuple[list[int], list[int]]:
    """Generate basis vectors for an alphabet size of `num_chars`."""
    num_basis = (num_chars - 1).bit_length()
    basis_vectors = [random.getrandbits(64) for _ in range(num_basis)]
    base_seed = random.getrandbits(64)
    seeds = []
    for i in range(num_chars):
        selection = (v for j, v in enumerate(basis_vectors) if (i >> j) & 1)
        seeds.append(functools.reduce(operator.xor, selection, base_seed))
    return seeds, basis_vectors


def build_matrix(
    basis_vectors: list[int],
    roll_function: typing.Callable[[int, int], int],
) -> list[int]:
    """
    Builds the binary transformation matrix.
    Each row is represented directly as a 64-bit integer.
    """
    d = len(basis_vectors)
    max_rotations = 64 // (d.bit_length() if d > 0 else 1)
    rv_iter = itertools.product(range(max_rotations), basis_vectors)
    return [roll_function(v, r) if r > 0 else v for r, v in rv_iter]


def get_rank(matrix: list[int]) -> int:
    """
    Computes the rank of a binary matrix over GF(2) using Gaussian elimination.
    matrix is a list of 64-bit integers representing rows.
    """
    rank = 0
    for col in range(64):
        mask = 1 << (63 - col)
        pivot_idx = -1
        for row in range(rank, len(matrix)):
            if matrix[row] & mask:
                pivot_idx = row
                break
        if pivot_idx == -1:
            continue
        matrix[rank], matrix[pivot_idx] = matrix[pivot_idx], matrix[rank]
        for row in range(len(matrix)):
            if row != rank and (matrix[row] & mask):
                matrix[row] ^= matrix[rank]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def generate_optimized_seeds(num_seeds: int, rng_seed: int) -> list[int]:
    random.seed(rng_seed)
    d = (num_seeds - 1).bit_length()
    expected_rank = d * (64 // d)
    seeds, rank, rotl_rank = [], 0, 0
    while rank < expected_rank or rotl_rank < expected_rank:
        seeds, basis_vectors = generate(num_seeds)
        rank = get_rank(build_matrix(basis_vectors, pynthash.roll_next))
        rotl_rank = get_rank(build_matrix(basis_vectors, pynthash.rotl))
    return seeds


if __name__ == "__main__":
    print(os.linesep.join(map(hex, generate_optimized_seeds(4, 42))))
