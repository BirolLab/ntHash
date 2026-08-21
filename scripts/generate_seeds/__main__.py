import itertools
import os
import random
import typing

import pynthash


def generate(num_seeds: int) -> tuple[list[int], list[int]]:
    """Generates random 64-bit seeds and their relative basis vectors."""
    # Generate num_seeds random 64-bit unsigned integers
    seeds = [random.getrandbits(64) for _ in range(num_seeds)]
    # Basis vectors are the differences (XOR) relative to the first seed
    basis_vectors = [seeds[0] ^ seeds[i] for i in range(1, num_seeds)]
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
    seeds, rank = [], 0
    while rank < expected_rank:
        seeds, basis_vectors = generate(num_seeds)
        matrix = build_matrix(basis_vectors, pynthash.roll)
        rank = get_rank(matrix)
    return seeds


if __name__ == "__main__":
    print(os.linesep.join(map(hex, generate_optimized_seeds(4, 42))))
