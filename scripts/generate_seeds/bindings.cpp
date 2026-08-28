#include <pybind11/pybind11.h>

#include "nthash.hpp"

PYBIND11_MODULE(pynthash, m, pybind11::mod_gil_not_used())
{
  m.doc() = "ntHash Python bindings for seed generation";
  m.def("rotl",
        pybind11::overload_cast<uint64_t, unsigned>(&nthash::internal::rotl),
        "Bitwise rotation function",
        pybind11::arg("x"),
        pybind11::arg("d"));
  m.def(
    "roll_next",
    pybind11::overload_cast<uint64_t, unsigned>(&nthash::internal::roll_next),
    "Base rolling function",
    pybind11::arg("x"),
    pybind11::arg("d"));
}