#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include "dartsflash/stability/stability.hpp"

namespace py = pybind11;

void pybind_stability(py::module& m)
{
    using namespace pybind11::literals;  // bring in '_a' literal
    
    // Expose Stability class
    py::class_<Stability>(m, "Stability", R"pbdoc(
            This is the base class for performing stability tests.
            )pbdoc")
        .def(py::init<FlashParams&>(), R"pbdoc(
            :param flash_params: FlashParams object
            )pbdoc", "flash_params"_a)

        .def("init", &Stability::init, R"pbdoc(
            Initialise stability algorithm at p, T, with reference phase composition x
            
            :param x: List of reference phase composition.
            :type x: list
            )pbdoc", "x"_a)

        .def("run", &Stability::run, R"pbdoc(
            Run stability test from trial phase composition initial guess Y.
            
            :param Y: List of trial phase composition
            :type Y: list
            )pbdoc", "Y"_a)
        ;
}