import argparse
import dataclasses
import functools
import operator
import os
import random
import typing


@dataclasses.dataclass
class ProgramArguments:
    alphabet_size: int
    rng: typing.Optional[int]
    use_rol: bool

    @staticmethod
    def parse() -> "ProgramArguments":
        parser = argparse.ArgumentParser(
            description="Generate a collision-free seed set",
        )
        parser.add_argument(
            "-n",
            dest="alphabet_size",
            help="character alphabet size (number of seeds to generate)",
            type=int,
            default=4,
        )
        parser.add_argument(
            "-r",
            dest="rng",
            help="pseudo-random number generator seed",
            type=int,
            default=42,
        )
        parser.add_argument(
            "--rol",
            dest="use_rol",
            help="use rol instead of srol",
            action=argparse.BooleanOptionalAction,
        )
        return ProgramArguments(**vars(parser.parse_args()))


def generate(num_chars: int) -> tuple[list[int], list[int]]:
    """Generate basis vectors for an alphabet size of `num_chars`."""
    num_basis = (num_chars - 1).bit_length()
    basis_vectors = [random.getrandbits(64) for _ in range(num_basis)]
    base_seed = random.getrandbits(64)
    seeds = []
    for i in range(num_chars):
        selection = (v for j, v in enumerate(basis_vectors) if (i >> j) & 1)
        seeds.append(functools.reduce(operator.xor, selection, base_seed))
    return seeds, basis_vectors


def srol(val: int, shift: int) -> int:
    """ntHash2 Split Rotate Left (srol)."""
    if shift == 0:
        return val
    bits_high = 31
    bits_low = 33
    mask_low = (1 << bits_low) - 1
    mask_high = (1 << bits_high) - 1
    val_low = val & mask_low
    val_high = (val >> bits_low) & mask_high
    shift_low = shift % bits_low
    shift_high = shift % bits_high
    rot_low = ((val_low << shift_low) & mask_low) | (val_low >> (bits_low - shift_low))
    rot_high = ((val_high << shift_high) & mask_high) | (
        val_high >> (bits_high - shift_high)
    )
    return (rot_high << bits_low) | rot_low


def rol(val: int, shift: int, bits: int = 64) -> int:
    """Cyclic left shift for a 64-bit integer."""
    if shift == 0:
        return val
    mask = (1 << bits) - 1
    return ((val << shift) & mask) | (val >> (bits - shift))


def build_matrix(
    basis_vectors: list[int],
    rot_function: typing.Callable[[int, int], int],
) -> list[int]:
    max_k = 64 // len(basis_vectors)
    matrix = []
    for seed in basis_vectors:
        matrix.extend([rot_function(seed, shift) for shift in range(max_k)])
    return matrix


def get_rank(matrix: list[int]) -> int:
    """Calculates the rank of a GF(2) matrix represented as a list of 64-bit integers."""
    rank = 0
    num_rows = len(matrix)
    for col in range(63, -1, -1):
        if rank >= num_rows:
            break
        pivot = -1
        for row in range(rank, num_rows):
            if (matrix[row] >> col) & 1:
                pivot = row
                break
        if pivot == -1:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for row in range(rank + 1, num_rows):
            if (matrix[row] >> col) & 1:
                matrix[row] ^= matrix[rank]
        rank += 1
    return rank


def main() -> None:
    args = ProgramArguments.parse()
    random.seed(args.rng)
    d = (args.alphabet_size - 1).bit_length()
    expected_rank = d * (64 // d)
    seeds, rank = [], 0
    while rank < expected_rank:
        seeds, basis_vectors = generate(args.alphabet_size)
        matrix = build_matrix(basis_vectors, rol if args.use_rol else srol)
        rank = get_rank(matrix)
    print(os.linesep.join(map(str, seeds)))


if __name__ == "__main__":
    main()
