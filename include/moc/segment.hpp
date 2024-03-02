#ifndef SEGMENT_H
#define SEGMENT_H

#include <transport_xs.hpp>
#include <xtensor/xview.hpp>

class Segment {
  public:
    Segment();

    double length() const { return length_; }
    const TransportXS& xs() const { return *xs_; }

  private:
    xt::xview<double> fsr_;
    TransportXS* xs_;
    double length_;
};

#endif
