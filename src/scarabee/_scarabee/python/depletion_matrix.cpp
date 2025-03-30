#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pytensor.hpp>

#include <data/depletion_matrix.hpp>

#include <memory>

namespace py = pybind11;

using namespace scarabee;

void init_DepletionMatrix(py::module& m) {
  py::class_<DepletionMatrix, std::shared_ptr<DepletionMatrix>>(
      m, "DepletionMatrix",
      "Represents a depletion matrix containing transfer terms due to "
      "radioactive decay and transmutation reactions.")

      .def(
          py::init<const std::vector<std::string>& /*nuclides*/>(),
          "Initializes an empty depletion matrix for the provided nuclides.\n\n"
          "Parameters\n"
          "----------\n"
          "nuclides : list of string\n"
          "    All possible nuclides in the material.\n",
          py::arg("nuclides"))

      .def_property_readonly("size", &DepletionMatrix::size,
                             "Number of nuclides in the depletion matrix.")

      .def_property_readonly(
          "nuclides", &DepletionMatrix::nuclides,
          "List of nuclides considered in the depletion matrix.")

      .def("has_nuclide", &DepletionMatrix::has_nuclide,
           "Returns True if the depletion matrix contains the specified "
           "nuclide, and False otherwise.\n\n"
           "Parameters\n"
           "----------\n"
           "nuclide : string\n"
           "    Name of the nuclide.\n\n"
           "Returns\n"
           "-------\n"
           "bool\n"
           "    True if the nuclide is present, False otherwise.\n",
           py::arg("nuclide"))

      .def("get_nuclide_index", &DepletionMatrix::get_nuclide_index,
           "Returns the index of the specified nuclide, if contained in the "
           "depletion matrix.\n\n"
           "Parameters\n"
           "----------\n"
           "nuclide : string\n"
           "    Name of the nuclide.\n\n"
           "Returns\n"
           "-------\n"
           "int\n"
           "    Index of the nuclide.\n",
           py::arg("nuclide"))

      .def("zero", &DepletionMatrix::zero,
           "Sets all entries in the depletion matrix to zero.")

      .def_property_readonly(
          "is_compressed", &DepletionMatrix::is_compressed,
          "True if the matrix is in compressed format, False otherwise.")

      .def("compress", &DepletionMatrix::compress,
           "Puts the matrix into compressed format.")

      .def(
          "exponential_product",
          [](const DepletionMatrix& m, xt::pytensor<double, 1>& N,
             bool cram48 = false) {
            std::span<double> Nspn(N.data(), N.size());
            m.exponential_product(Nspn, cram48);
          },
          "Computes the matrix product of the input vector N with the "
          "exponential of the depletion matrix. The array N is modified in "
          "place to contain the result of this product.\n\n"
          ".. math:: N' = \\exp{A} N.\n\n"
          "This products is approximated using a Chebyshev Rational "
          "Approximation (CRAM). By default, at 16th order CRAM is used, but a "
          "48th order approximation is available.\n\n"
          "Parameters\n"
          "----------\n"
          "N : ndarray\n"
          "    1D array of the number densities of the nuclides in the "
          "depletion matrix (in the same order). This array is modified in "
          "place, and will contain the result of the product after the "
          "function is complete.\n"
          "cram48 : bool\n"
          "    If True, a 48th order CRAM is used. If False, a 16th order CRAM "
          "is employed. Default value is False.\n",
          py::arg("N"), py::arg("cram48") = false)

      .def("__getitem__",
           [](const DepletionMatrix& m,
              std::pair<std::size_t, std::size_t> indx) {
             return m.value(indx.first, indx.second);
           })

      .def("__setitem__",
           [](DepletionMatrix& m, std::pair<std::size_t, std::size_t> indx,
              double v) { m.ref(indx.first, indx.second) = v; })

      .def("__imul__", &DepletionMatrix::operator*=)
      .def("__idiv__", &DepletionMatrix::operator/=)
      .def("__iadd__", &DepletionMatrix::operator+=)
      .def("__isub__", &DepletionMatrix::operator-=)
      .def("__add__", &DepletionMatrix::operator+)
      .def("__sub__", &DepletionMatrix::operator-);
}
