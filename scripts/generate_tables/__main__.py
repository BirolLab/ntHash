import argparse
import dataclasses
import functools
import operator
import os
import random
import typing


@dataclasses.dataclass
class CommandLineArguments:
    num_chars: int
    rng_seed: int

    @staticmethod
    def parse() -> "CommandLineArguments":
        parser = argparse.ArgumentParser(
            description="Generate a collision-free seed set and update internals.hpp",
        )
        parser.add_argument(
            "-n",
            dest="num_chars",
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
        return CommandLineArguments(**vars(parser.parse_args()))


def main() -> None:
    args = CommandLineArguments.parse()
    random.seed(args.rng_seed)
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
