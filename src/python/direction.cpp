#include <pybind11/pybind11.h>

#include <moc/direction.hpp>

#include <sstream>

namespace py = pybind11;

using namespace scarabee;

void init_Direction(py::module& m) {
  py::class_<Direction, Vector>(m, "Direction")
      .def(py::init<double, double>(),
           "A vector in R2 of unit length.\n\n"
           "Arguments:\n"
           "    x  x-axis component"
           "    y  y-axis component",
           py::arg("x"), py::arg("y"))

      .def(py::init<double>(),
           "A vector in R2 of unit length.\n\n"
           "Arguments:\n"
           "    phi  angle with x-axis, in radians",
           py::arg("phi"))

      .def("__add__",
           [](const Direction& d1, const Direction& d2) { return d1 + d2; })

      .def("__add__", [](const Direction& d, const Vector& v) { return d + v; })

      .def("__radd__",
           [](const Direction& d, const Vector& v) { return v + d; })

      .def("__sub__",
           [](const Direction& d1, const Direction& d2) { return d1 - d2; })

      .def("__sub__", [](const Direction& d, const Vector& v) { return d - v; })

      .def("__rsub__",
           [](const Direction& d, const Vector& v) { return v - d; })

      .def("__mul__",
           [](const Direction& d1, const Direction& d2) { return d1 * d2; })

      .def("__mul__", [](const Direction& d, const Vector& v) { return d * v; })

      .def("__rmul__",
           [](const Direction& d, const Vector& v) { return v * d; })

      .def("__mul__", [](const Direction& u, double d) { return u * d; })

      .def("__rmul__", [](const Direction& u, double d) { return d * u; })

      .def("__truediv__", [](const Direction& u, double d) { return u / d; })

      .def("__repr__", [](const Vector& d) {
        std::stringstream mssg;
        mssg << "<<" << d.x() << "," << d.y() << ">>";
        return mssg.str();
      });
}
