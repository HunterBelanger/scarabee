#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/cartesian_2d.hpp>

#include <memory>

namespace py = pybind11;

using namespace scarabee;

void init_Cartesian2D(py::module& m) {
  // Tile
  py::class_<Cartesian2D::Tile>(m, "Tile")
      .def_readwrite("c2d", &Cartesian2D::Tile::c2d)
      .def_readwrite("cell", &Cartesian2D::Tile::cell)
      .def("valid", &Cartesian2D::Tile::valid);

  // TileIndex
  py::class_<Cartesian2D::TileIndex>(m, "TileIndex")
      .def_readwrite("i", &Cartesian2D::TileIndex::i)
      .def_readwrite("j", &Cartesian2D::TileIndex::j);

  // Cartesian2D
  py::class_<Cartesian2D, std::shared_ptr<Cartesian2D>>(m, "Cartesian2D")
      .def(py::init<const std::vector<double>& /*dx*/,
                    const std::vector<double>& /*dy*/>(),
           "Creates a 2D cartesian geometry, with empty tiles.\n\n"
           "Arguments:\n"
           "    dx  List of all x widths\n"
           "    dy  List of all y heights",
           py::arg("dx"), py::arg("dy"))

      .def_property_readonly("nx", &Cartesian2D::nx,
                             "Number of cells in x direction.")

      .def_property_readonly("ny", &Cartesian2D::ny,
                             "Number of cells in y direction.")

      .def_property_readonly("dx", &Cartesian2D::dx, "Width of geometry in x.")

      .def_property_readonly("dy", &Cartesian2D::dy, "Width of geometry in y.")

      .def_property_readonly("x_min", &Cartesian2D::x_min,
                             "Minimum x coordinate.")

      .def_property_readonly("x_max", &Cartesian2D::x_max,
                             "Maximum x coordinate.")

      .def_property_readonly("y_min", &Cartesian2D::y_min,
                             "Minimum y coordinate.")

      .def_property_readonly("y_max", &Cartesian2D::y_max,
                             "Maximum y coordinate.")

      .def("tiles_valid", &Cartesian2D::tiles_valid,
           "Returns true if all tiles are populated.")

      .def("get_tile_index", &Cartesian2D::get_tile_index,
           "Returns TileIndex for given position and direction.\n\n"
           "Arguments:\n"
           "    r  position Vector\n"
           "    u  Direction",
           py::arg("r"), py::arg("u"))

      .def("get_xs", &Cartesian2D::get_xs,
           "Returns TransportXS for given position and direction.\n\n"
           "Arguments:\n"
           "    r  position Vector\n"
           "    u  Direction",
           py::arg("r"), py::arg("u"))

      .def("tile", &Cartesian2D::tile,
           "Returns the Tile for the given TileIndex.\n\n"
           "Arguments:\n"
           "    ti  TileIndex for desired Tile.",
           py::arg("ti"))

      .def("set_tile",
           py::overload_cast<const Cartesian2D::TileIndex&,
                             const std::shared_ptr<Cartesian2D>&>(
               &Cartesian2D::set_tile),
           "Fills tile with provided Cartesian2D.\n\n"
           "Arguments:\n"
           "    ti   TileIndex\n"
           "    c2d  Cartesian2D to fill Tile.",
           py::arg("ti"), py::arg("c2d"))

      .def("set_tile",
           py::overload_cast<const Cartesian2D::TileIndex&,
                             const std::shared_ptr<Cell>&>(
               &Cartesian2D::set_tile),
           "Fills tile with provided Cell.\n\n"
           "Arguments:\n"
           "    ti    TileIndex\n"
           "    cell  Cell to fill Tile.",
           py::arg("ti"), py::arg("cell"))

      .def("set_tiles", &Cartesian2D::set_tiles,
           "Sets all tiles in the geometry.\n\n"
           "Arguments:\n"
           "    fills  List of Cartesian2D or Cell instances to fill all tiles",
           py::arg("fills"));
}
