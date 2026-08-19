py_mod = import('python')
py3 = py_mod.find_installation('python3', required: true)
pybind11_dep = dependency('pybind11', required: true)

pynthash = py3.extension_module(
    'pynthash',
    'bindings.cpp',
    dependencies: pybind11_dep,
    include_directories: include_dirs + ['../../src'],
)
