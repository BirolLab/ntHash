#!/usr/bin/env bash

diff -s <(echo "4362857412768957556
3572411708064410444
2319985823310095140
2978368046464386134" | python make_internal_hpp.py | clang-format --style=file) known.hpp