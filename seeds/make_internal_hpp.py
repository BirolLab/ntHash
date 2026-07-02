import functools
import itertools
import operator
import pathlib
import string
import sys
import typing

DEFAULT_TEMPLATE_PATH = pathlib.Path(__file__).parent / "internal.hpp.in"


def rol(val, r_bits, max_bits):
    r_bits = r_bits % max_bits
    return ((val << r_bits) | (val >> (max_bits - r_bits))) & ((1 << max_bits) - 1)


def get_substitutions(seeds: list[int]) -> typing.Dict[str, str]:
    subs = {}
    l31, r33 = {}, {}
    for seed, base in zip(seeds, "ACGT"):
        l31[base] = [rol(seed >> 33, i, 31) << 33 for i in range(31)]
        r33[base] = [rol(seed & ((1 << 33) - 1), i, 33) for i in range(33)]
        subs[f"SEED_{base}"] = f"{seed:#x}"
        subs[f"{base}31L"] = f"{{ {','.join(f'{a:#x}' for a in l31[base])} }}"
        subs[f"{base}33R"] = f"{{ {','.join(f'{a:#x}' for a in r33[base])} }}"
    for n, key in enumerate(("DI", "TRI", "TETRA"), start=2):
        tab = itertools.product("ACGT", repeat=n)
        tab = [[l31[b][i] ^ r33[b][i] for i, b in enumerate(reversed(s))] for s in tab]
        tab = [functools.reduce(operator.xor, x) for x in tab]
        subs[f"{key}MER_TAB"] = f"{{ {','.join(f'{a}U' for a in tab)} }}"
    return subs


def main():
    seeds = [int(input()) for _ in range(4)]
    hpp_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_TEMPLATE_PATH
    template = string.Template(hpp_path.read_text())
    print(template.substitute(**get_substitutions(seeds)), end="")


if __name__ == "__main__":
    main()
