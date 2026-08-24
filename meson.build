project(
    'ntHash',
    'cpp',
    version: '3.0.0',
    license: 'MIT',
    default_options: [
        'cpp_std=c++17',
        'werror=true',
        'warning_level=3',
        'optimization=3',
    ],
)

python3 = find_program('python3', required: true)
clang_format = find_program('clang-format', required: true)

nthash_headers = files(
  'src/internal/internal.hpp',
  'src/internal/tables.hpp',
  'src/internal/utils.hpp',
  'src/kmer/blind.hpp',
  'src/kmer/nthash.hpp',
  'src/kmer/utils.hpp',
  'src/seed/blind.hpp',
  'src/seed/seed.hpp',
  'src/seed/utils.hpp',
  'src/nthash.hpp',
)

single_header = custom_target(
  'single_header',
  input: nthash_headers,
  output: 'nthash.hpp.tmp',
  command: [python3, files('scripts/amalgamate.py'), '@INPUT@'],
  capture: true,
)

formatted_header = custom_target(
  'formatted_header',
  input: single_header,
  output: 'nthash.hpp',
  command: [clang_format, '--style=file', '@INPUT@'],
  capture: true,           
  install: true,           
  install_dir: get_option('includedir'),
  build_by_default: true,
)

nthash_dep = declare_dependency(
  include_directories: include_directories('.'),
  sources: formatted_header
)

clang_tidy = find_program('clang-tidy', required: false)
if clang_tidy.found()
    run_target(
      'clang-tidy',
      command: [clang_tidy, '-p=build', nthash_headers]
    )
endif

test_exec = executable(
    'tests',
    'tests/tests.cpp',
    dependencies: [nthash_dep],
)
test('tests', test_exec)

doxygen = find_program('doxygen', required: false)
if doxygen.found()
    custom_target(
      'docs',
      input: 'Doxyfile',
      depends: formatted_header,
      output: 'docs',
      command: [doxygen, '@INPUT@'],
      env: {
        'DOXY_NTHASH_INPUT': meson.current_build_dir() / 'nthash.hpp',
        'DOXY_NTHASH_OUTPUT': meson.project_source_root() / 'docs'
      },
      build_by_default: false,
    )
else
    warning('Documentation will not be updated without doxygen')
endif

subdir('scripts')

