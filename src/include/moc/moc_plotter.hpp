#ifndef MOC_PLOTTER_H
#define MOC_PLOTTER_H

#include <moc/moc_driver.hpp>

#include <ImApp/imapp.hpp>

namespace scarabee {

class MOCPlotter : public ImApp::Layer {
 public:
  MOCPlotter(const MOCDriver* moc);

  void render() override final;

 private:
  void render_viewport();
  void render_controls();

  ImApp::Pixel get_color(UniqueFSR ufsr);
  ImApp::Pixel get_random_color();

  Direction get_tracking_direction() const;
  Vector get_start_position(uint64_t i) const;
  Direction get_comp_tracking_direction() const;
  Vector get_comp_start_position(uint64_t j) const;

  void render_image();

  enum ColorBy : int { Cell = 0, Material = 1 };

  // Maps for colors from plotter.hpp
  const MOCDriver* moc_;
  const Cartesian2D* geom_;
  std::map<uint32_t, ImApp::Pixel> cell_id_to_color;
  std::map<CrossSection*, ImApp::Pixel> material_id_to_color;
  ImApp::Image image;
  std::mutex create_color_mutex;
  int adjust_w_or_h;
  double height, width;   // Width of image in physical space ([cm]).
  double ox, oy;      // Plot origin
  double mx, my;      // Mouse position
  UniqueFSR mufsr;        // Mouse UniqueFSR
  std::size_t mfsri;      // Mouse FSR index
  uint32_t mcell_inst;    // Mouse cell instance
  CrossSection* mxs;      // Mouse CrossSection
  double dist_per_pixel;
  ImApp::Pixel background;
  ColorBy colorby;
  bool must_rerender;
  bool outline_boundaries;
};

}

#endif