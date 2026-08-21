py_mod = import('python')
py3 = py_mod.find_installation('python3', required: true)
pybind11_dep = dependency('pybind11', required: true)

pynthash = py3.extension_module(
    'pynthash',
    'bindings.cpp',
    dependencies: [pybind11_dep, nthash_dep],
)

foreach src : ['__init__.py', '__main__.py']
    configure_file(
        input: src,
        output: src,
        copy: true
    )
endforeach