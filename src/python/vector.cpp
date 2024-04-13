#include <pybind11/pybind11.h>

#include <moc/vector.hpp>

#include <sstream>

namespace py = pybind11;

using namespace scarabee;

void init_Vector(py::module& m) {
  py::class_<Vector>(m, "Vector")
      .def(py::init<double, double>(),
           "A vector in R2 of arbitrary length.\n\n"
           "Arguments:\n"
           "    x  x-axis component\n"
           "    y  y-axis component",
           py::arg("x"), py::arg("y"))

      .def("dot", &Vector::dot,
           "Computes dot product with another vector.\n\n"
           "Arguments:\n"
           "    v  other vector for dot product",
           py::arg("v"))

      .def_property_readonly("x", &Vector::x, "y component of vector")

      .def_property_readonly("y", &Vector::y, "y component of vector")

      .def_property_readonly("norm", &Vector::norm, "Length of the vector")

      .def("__neg__", &Vector::operator-)

      .def("__eq__", &Vector::operator==)

      .def("__ne__", &Vector::operator!=)

      .def("__add__",
           [](const Vector& v1, const Vector& v2) { return v1 + v2; })

      .def("__sub__",
           [](const Vector& v1, const Vector& v2) { return v1 - v2; })

      .def("__mul__",
           [](const Vector& v1, const Vector& v2) { return v1 * v2; })

      .def("__mul__", [](const Vector& v, double d) { return v * d; })

      .def("__rmul__", [](const Vector& v, double d) { return d * v; })

      .def("__truediv__", [](const Vector& v, double d) { return v / d; })

      .def("__repr__", [](const Vector& v) {
        std::stringstream mssg;
        mssg << "<" << v.x() << "," << v.y() << ">";
        return mssg.str();
      });
}
