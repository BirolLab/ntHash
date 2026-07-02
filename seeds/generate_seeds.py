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
        )
        return ProgramArguments(**vars(parser.parse_args()))


def generate(num_chars: int) -> list[int]:
    num_basis = (num_chars - 1).bit_length()
    basis_vectors = [random.getrandbits(64) for _ in range(num_basis)]
    base_seed = random.getrandbits(64)
    seeds = []
    for i in range(num_chars):
        selection = (v for j, v in enumerate(basis_vectors) if (i >> j) & 1)
        seeds.append(functools.reduce(operator.xor, selection, base_seed))
    return seeds


def main() -> None:
    args = ProgramArguments.parse()
    random.seed(args.rng)
    print(os.linesep.join(map(str, generate(args.alphabet_size))))


if __name__ == "__main__":
    main()
