#include <moc/surface.hpp>
#include <utils/constants.hpp>

#include <cmath>

namespace scarabee {

//=============================================================================
// XPlane
//-----------------------------------------------------------------------------

inline Surface::Side xplane_side(const double x0, const Vector& r,
                                 const Direction& u) {
  if (r.x() - x0 > SURFACE_COINCIDENT)
    return Surface::Side::Positive;
  else if (r.x() - x0 < -SURFACE_COINCIDENT)
    return Surface::Side::Negative;
  else {
    if (u.x() > 0.)
      return Surface::Side::Positive;
    else
      return Surface::Side::Negative;
  }
}

inline double xplane_distance(const double x0, const Vector& r,
                              const Direction& u) {
  const double diff = x0 - r.x();
  if (std::abs(diff) < SURFACE_COINCIDENT || u.x() == 0.)
    return INF;
  else if (diff / u.x() < 0.)
    return INF;
  else
    return diff / u.x();
}

inline double xplane_integrate_x(const double /*x0*/, const double /*xmin*/,
                                 const double /*xmax*/,
                                 const Surface::Side /*side*/) {
  // An XPlane is perpendicular to the x-axis; dx integration makes no sense.
  return 0.;
}

inline double xplane_integrate_y(const double x0, const double xmin,
                                 const double xmax,
                                 const Surface::Side /*side*/) {
  return x0 * (xmax - xmin);
}

//=============================================================================
// YPlane
//-----------------------------------------------------------------------------

inline Surface::Side yplane_side(const double y0, const Vector& r,
                                 const Direction& u) {
  if (r.y() - y0 > SURFACE_COINCIDENT)
    return Surface::Side::Positive;
  else if (r.y() - y0 < -SURFACE_COINCIDENT)
    return Surface::Side::Negative;
  else {
    if (u.y() > 0.)
      return Surface::Side::Positive;
    else
      return Surface::Side::Negative;
  }
}

inline double yplane_distance(const double y0, const Vector& r,
                              const Direction& u) {
  const double diff = y0 - r.y();
  if (std::abs(diff) < SURFACE_COINCIDENT || u.y() == 0.)
    return INF;
  else if (diff / u.y() < 0.)
    return INF;
  else
    return diff / u.y();
}

inline double yplane_integrate_x(const double y0, const double xmin,
                                 const double xmax,
                                 const Surface::Side /*side*/) {
  return y0 * (xmax - xmin);
}

inline double yplane_integrate_y(const double /*y0*/, const double /*ymin*/,
                                 const double /*ymax*/,
                                 const Surface::Side /*side*/) {
  // A YPlane is perpendicular to the y-axis; dy integration makes no sense.
  return 0.;
}

//=============================================================================
// Plane
//-----------------------------------------------------------------------------

inline Surface::Side plane_side(const double A, const double B, const double C,
                                const Vector& r, const Direction& u) {
  const double eval = A * r.x() + B * r.y() - C;
  if (eval > SURFACE_COINCIDENT)
    return Surface::Side::Positive;
  else if (eval < -SURFACE_COINCIDENT)
    return Surface::Side::Negative;
  else {
    if ((A * u.x() + B * u.y()) > 0.)
      return Surface::Side::Positive;
    else
      return Surface::Side::Negative;
  }
}

inline double plane_distance(const double A, const double B, const double C,
                             const Vector& r, const Direction& u) {
  const double num = C - A * r.x() - B * r.y();
  const double denom = A * u.x() + B * u.y();
  const double d = num / denom;
  if (std::abs(d) < SURFACE_COINCIDENT || denom == 0.)
    return INF;
  else if (d < 0.)
    return INF;
  else
    return d;
}

inline double plane_integrate_x(const double A, const double B, const double C,
                                const double xmin, const double xmax,
                                const Surface::Side /*side*/) {
  if (B == 0.) return 0.;

  auto F = [&](double x) { return (C / B) * x - 0.5 * (A / B) * x * x; };
  return F(xmax) - F(xmin);
}

inline double plane_integrate_y(const double A, const double B, const double C,
                                const double ymin, const double ymax,
                                const Surface::Side /*side*/) {
  if (A == 0.) return 0.;

  auto F = [&](double x) { return (C / A) * x - 0.5 * (B / A) * x * x; };
  return F(ymax) - F(ymin);
}

//=============================================================================
// Cylinder
//-----------------------------------------------------------------------------

inline Surface::Side cylinder_side(const double x0, const double y0,
                                   const double rc, const Vector& r,
                                   const Direction& u) {
  const double x = r.x() - x0;
  const double y = r.y() - y0;
  const double eval = y * y + x * x - rc * rc;
  if (eval > SURFACE_COINCIDENT)
    return Surface::Side::Positive;
  else if (eval < -SURFACE_COINCIDENT)
    return Surface::Side::Negative;
  else {
    Direction norm(r.x() - x0, r.y() - y0);
    if (u.dot(norm) > 0.)
      return Surface::Side::Positive;
    else
      return Surface::Side::Negative;
  }
}

inline double cylinder_distance(const double x0, const double y0,
                                const double rc, const Vector& r,
                                const Direction& u) {
  const double a = u.y() * u.y() + u.x() * u.x();
  if (a == 0.) return INF;

  const double x = r.x() - x0;
  const double y = r.y() - y0;
  const double k = y * u.y() + x * u.x();
  const double c = y * y + x * x - rc * rc;
  const double quad = k * k - a * c;

  if (quad < 0.)
    return INF;
  else if (std::abs(c) < SURFACE_COINCIDENT) {
    if (k >= 0.)
      return INF;
    else
      return (-k + std::sqrt(quad)) / a;
  } else if (c < 0.) {
    return (-k + std::sqrt(quad)) / a;
  } else {
    const double d = (-k - std::sqrt(quad)) / a;
    if (d < 0.)
      return INF;
    else
      return d;
  }
}

inline double cylinder_integrate_x(const double x0, const double y0,
                                   const double rc, const double xmin,
                                   const double xmax,
                                   const Surface::Side side) {
  const double r2 = rc * rc;

  auto F = [&](double x) {
    const double dx = x - x0;
    const double dx2 = dx * dx;
    const double sqrt_r2_dx2 = std::sqrt(r2 - dx2);

    const double sqrt_int =
        0.5 * (dx * sqrt_r2_dx2 + r2 * std::atan(dx / sqrt_r2_dx2));

    double out = x * y0;
    if (side == Surface::Side::Negative) {
      out -= sqrt_int;
    } else {
      out += sqrt_int;
    }

    return out;
  };

  const double new_xmin = std::max(xmin, x0 - rc);
  const double new_xmax = std::min(xmax, x0 + rc);

  return F(new_xmax) - F(new_xmin);
}

inline double cylinder_integrate_y(const double x0, const double y0,
                                   const double rc, const double ymin,
                                   const double ymax,
                                   const Surface::Side side) {
  const double r2 = rc * rc;

  auto F = [&](double y) {
    const double dy = y - y0;
    const double dy2 = dy * dy;
    const double sqrt_r2_dy2 = std::sqrt(r2 - dy2);

    const double sqrt_int =
        0.5 * (dy * sqrt_r2_dy2 + r2 * std::atan(dy / sqrt_r2_dy2));

    double out = y * x0;
    if (side == Surface::Side::Negative) {
      out -= sqrt_int;
    } else {
      out += sqrt_int;
    }

    return out;
  };

  const double new_ymin = std::max(ymin, y0 - rc);
  const double new_ymax = std::min(ymax, y0 + rc);

  return F(new_ymax) - F(new_ymin);
}

//=============================================================================
// BWRBox
//-----------------------------------------------------------------------------

inline Surface::Side bwr_box_side(const double xl, const double xh,
                                  const double yl, const double yh,
                                  const double Rcx, const double Rcy,
                                  const double rc, const Surface::Type type,
                                  const Vector& r, const Direction& u) {
  // First, we check if we are inside the cell at all
  const auto side_xl = xplane_side(xl, r, u);
  const auto side_xh = xplane_side(xh, r, u);
  const auto side_yl = yplane_side(yl, r, u);
  const auto side_yh = yplane_side(yh, r, u);

  const bool in_box = side_xl == Surface::Side::Positive &&
                      side_xh == Surface::Side::Negative &&
                      side_yl == Surface::Side::Positive &&
                      side_yh == Surface::Side::Negative;
  if (in_box == false) return Surface::Side::Positive;

  // Next, check if we are in the circle box
  const auto side_Rcx = xplane_side(Rcx, r, u);
  const auto side_Rcy = yplane_side(Rcy, r, u);

  // This statement is dependent on the quadrant !
  const bool in_circle_box = [&]() {
    if (type == Surface::Type::BWRCornerI) {
      return side_Rcx == Surface::Side::Positive &&
             side_xh == Surface::Side::Negative &&
             side_Rcy == Surface::Side::Positive &&
             side_yh == Surface::Side::Negative;
    } else if (type == Surface::Type::BWRCornerII) {
      return side_xl == Surface::Side::Positive &&
             side_Rcx == Surface::Side::Negative &&
             side_Rcy == Surface::Side::Positive &&
             side_yh == Surface::Side::Negative;
    } else if (type == Surface::Type::BWRCornerIII) {
      return side_xl == Surface::Side::Positive &&
             side_Rcx == Surface::Side::Negative &&
             side_yl == Surface::Side::Positive &&
             side_Rcy == Surface::Side::Negative;
    } else {
      return side_Rcx == Surface::Side::Positive &&
             side_xh == Surface::Side::Negative &&
             side_yl == Surface::Side::Positive &&
             side_Rcy == Surface::Side::Negative;
    }
  }();
  if (in_circle_box == false) return Surface::Side::Negative;

  // Now, we check if we are inside the circle, and use that to determine if we
  // are on the positive or negative side.
  return cylinder_side(Rcx, Rcy, rc, r, u);
}

inline double bwr_box_distance(const double xl, const double xh,
                               const double yl, const double yh,
                               const double Rcx, const double Rcy,
                               const double rc, const Surface::Type type,
                               const Vector& r, const Direction& u) {
  const auto orig_side = bwr_box_side(xl, xh, yl, yh, Rcx, Rcy, rc, type, r, u);

  auto next_side = orig_side;
  Vector r_next(r);

  double d_sum = 0.;
  while (next_side == orig_side) {
    const double d_xl = xplane_distance(xl, r_next, u);
    const double d_xh = xplane_distance(xh, r_next, u);
    const double d_yl = yplane_distance(yl, r_next, u);
    const double d_yh = yplane_distance(yh, r_next, u);
    const double d_circ = cylinder_distance(Rcx, Rcy, rc, r_next, u);

    const double d = std::min({d_xl, d_xh, d_yl, d_yh, d_circ});

    if (d == INF) {
      // We will never hit the surface
      return d;
    }

    d_sum += d;
    r_next = Vector(r_next.x() + d * u.x(), r_next.y() + d * u.y());
    next_side = bwr_box_side(xl, xh, yl, yh, Rcx, Rcy, rc, type, r_next, u);
  }

  return d_sum;
}

inline double bwr_box_I_integrate_x(const double xl, const double xh,
                                    const double yl, const double yh,
                                    const double rc, double xmin, double xmax,
                                    const Surface::Side side) {
  // Truncate integration bounds if needed to be in range of the surface
  if (xmin < xl) xmin = xl;
  if (xmax > xh) xmax = xh;

  const double Rcx = xh - rc;
  const double Rcy = yh - rc;

  if (side == Surface::Side::Negative) {
    // Bottom surface is flat y-plane
    return yl * (xmax - xmin);
  }

  double int_out = 0.;

  // Top surface has the curve with two parts: straight line and circle
  if (xmin < Rcx) {
    int_out += yh * (std::min(Rcx, xmax) - xmin);
  }

  // Now we grab the curved portion
  if (xmax > Rcx) {
    int_out += cylinder_integrate_x(Rcx, Rcy, rc, std::max(Rcx, xmin), xmax,
                                    Surface::Side::Positive);
  }

  return int_out;
}

inline double bwr_box_I_integrate_y(const double xl, const double xh,
                                    const double yl, const double yh,
                                    const double rc, double ymin, double ymax,
                                    const Surface::Side side) {
  // Truncate integration bounds if needed to be in range of the surface
  if (ymin < yl) ymin = yl;
  if (ymax > yh) ymax = yh;

  const double Rcx = xh - rc;
  const double Rcy = yh - rc;

  if (side == Surface::Side::Negative) {
    // Bottom surface is flat x-plane
    return xl * (ymax - ymin);
  }

  double int_out = 0.;

  // Top surface has the curve with two parts: straight line and circle
  if (ymin < Rcy) {
    int_out += xh * (std::min(Rcy, ymax) - ymin);
  }

  // Now we grab the curved portion
  if (ymax > Rcy) {
    int_out += cylinder_integrate_y(Rcx, Rcy, rc, std::max(Rcy, ymin), ymax,
                                    Surface::Side::Positive);
  }

  return int_out;
}

inline double bwr_box_II_integrate_x(const double xl, const double xh,
                                     const double yl, const double yh,
                                     const double rc, double xmin, double xmax,
                                     const Surface::Side side) {
  // Truncate integration bounds if needed to be in range of the surface
  if (xmin < xl) xmin = xl;
  if (xmax > xh) xmax = xh;

  const double Rcx = xl + rc;
  const double Rcy = yh - rc;

  if (side == Surface::Side::Negative) {
    // Bottom surface is flat y-plane
    return yl * (xmax - xmin);
  }

  double int_out = 0.;

  // Curve
  if (xmin < Rcx) {
    int_out += cylinder_integrate_x(Rcx, Rcy, rc, xmin, std::min(Rcx, xmax),
                                    Surface::Side::Positive);
  }

  // Now we grab the straight portion
  if (xmax > Rcx) {
    int_out += yh * (xmax - std::max(Rcx, xmin));
  }

  return int_out;
}

inline double bwr_box_II_integrate_y(const double xl, const double xh,
                                     const double yl, const double yh,
                                     const double rc, double ymin, double ymax,
                                     const Surface::Side side) {
  // Truncate integration bounds if needed to be in range of the surface
  if (ymin < yl) ymin = yl;
  if (ymax > yh) ymax = yh;

  const double Rcx = xl + rc;
  const double Rcy = yh - rc;

  if (side == Surface::Side::Positive) {
    // Top surface is flat y-plane
    return xh * (ymax - ymin);
  }

  double int_out = 0.;

  // Straight
  if (ymin < Rcy) {
    int_out += xl * (std::min(Rcy, ymax) - ymin);
  }

  // Curve
  if (ymax > Rcy) {
    int_out += cylinder_integrate_y(Rcx, Rcy, rc, std::max(Rcy, ymin), ymax,
                                    Surface::Side::Negative);
  }

  return int_out;
}

inline double bwr_box_III_integrate_x(const double xl, const double xh,
                                      const double yl, const double yh,
                                      const double rc, double xmin, double xmax,
                                      const Surface::Side side) {
  // Truncate integration bounds if needed to be in range of the surface
  if (xmin < xl) xmin = xl;
  if (xmax > xh) xmax = xh;

  const double Rcx = xl + rc;
  const double Rcy = yl + rc;

  if (side == Surface::Side::Positive) {
    // Top surface is flat y-plane
    return yh * (xmax - xmin);
  }

  double int_out = 0.;

  // Curve
  if (xmin < Rcx) {
    int_out += cylinder_integrate_x(Rcx, Rcy, rc, xmin, std::min(Rcx, xmax),
                                    Surface::Side::Negative);
  }

  // Now we grab the straight portion
  if (xmax > Rcx) {
    int_out += yl * (xmax - std::max(Rcx, xmin));
  }

  return int_out;
}

inline double bwr_box_III_integrate_y(const double xl, const double xh,
                                      const double yl, const double yh,
                                      const double rc, double ymin, double ymax,
                                      const Surface::Side side) {
  // Truncate integration bounds if needed to be in range of the surface
  if (ymin < yl) ymin = yl;
  if (ymax > yh) ymax = yh;

  const double Rcx = xl + rc;
  const double Rcy = yl + rc;

  if (side == Surface::Side::Positive) {
    // Top surface is flat y-plane
    return xh * (ymax - ymin);
  }

  double int_out = 0.;

  // Curve
  if (ymin < Rcy) {
    int_out += cylinder_integrate_y(Rcx, Rcy, rc, ymin, std::min(Rcy, ymax),
                                    Surface::Side::Negative);
  }

  // Straight
  if (ymax > Rcy) {
    int_out += xl * (ymax - std::max(Rcy, ymin));
  }

  return int_out;
}

inline double bwr_box_IV_integrate_x(const double xl, const double xh,
                                     const double yl, const double yh,
                                     const double rc, double xmin, double xmax,
                                     const Surface::Side side) {
  // Truncate integration bounds if needed to be in range of the surface
  if (xmin < xl) xmin = xl;
  if (xmax > xh) xmax = xh;

  const double Rcx = xh - rc;
  const double Rcy = yl + rc;

  if (side == Surface::Side::Positive) {
    // Top surface is flat y-plane
    return yh * (xmax - xmin);
  }

  double int_out = 0.;

  // Straight
  if (xmin < Rcx) {
    int_out += yl * (std::min(Rcx, xmax) - xmin);
  }

  // Curve
  if (xmax > Rcx) {
    int_out += cylinder_integrate_x(Rcx, Rcy, rc, std::max(Rcx, xmin), xmax,
                                    Surface::Side::Negative);
  }

  return int_out;
}

inline double bwr_box_IV_integrate_y(const double xl, const double xh,
                                     const double yl, const double yh,
                                     const double rc, double ymin, double ymax,
                                     const Surface::Side side) {
  // Truncate integration bounds if needed to be in range of the surface
  if (ymin < yl) ymin = yl;
  if (ymax > yh) ymax = yh;

  const double Rcx = xh - rc;
  const double Rcy = yl + rc;

  if (side == Surface::Side::Negative) {
    // Bottom surface is flat x-plane
    return xl * (ymax - ymin);
  }

  double int_out = 0.;

  // Curve
  if (ymin < Rcy) {
    int_out += cylinder_integrate_y(Rcx, Rcy, rc, ymin, std::min(Rcy, ymax),
                                    Surface::Side::Positive);
  }

  // Straight
  if (ymax > Rcy) {
    int_out += yh * (ymax - std::max(Rcy, ymin));
  }

  return int_out;
}

//=============================================================================
// Surface
//-----------------------------------------------------------------------------

Surface::Side Surface::side(const Vector& r, const Direction& u) const {
  switch (type_) {
    case Type::XPlane:
      return xplane_side(this->x0(), r, u);
      break;

    case Type::YPlane:
      return yplane_side(this->y0(), r, u);
      break;

    case Type::Plane:
      return plane_side(this->A(), this->B(), this->C(), r, u);
      break;

    case Type::Cylinder:
      return cylinder_side(this->x0(), this->y0(), this->r(), r, u);
      break;

    case Type::BWRCornerI: {
      const double Rcx = this->xh() - this->rc();
      const double Rcy = this->yh() - this->rc();
      return bwr_box_side(this->xl(), this->xh(), this->yl(), this->yh(), Rcx,
                          Rcy, this->rc(), type_, r, u);
    } break;

    case Type::BWRCornerII: {
      const double Rcx = this->xl() + this->rc();
      const double Rcy = this->yh() - this->rc();
      return bwr_box_side(this->xl(), this->xh(), this->yl(), this->yh(), Rcx,
                          Rcy, this->rc(), type_, r, u);
    } break;

    case Type::BWRCornerIII: {
      const double Rcx = this->xl() + this->rc();
      const double Rcy = this->yl() + this->rc();
      return bwr_box_side(this->xl(), this->xh(), this->yl(), this->yh(), Rcx,
                          Rcy, this->rc(), type_, r, u);
    } break;

    case Type::BWRCornerIV: {
      const double Rcx = this->xh() - this->rc();
      const double Rcy = this->yl() + this->rc();
      return bwr_box_side(this->xl(), this->xh(), this->yl(), this->yh(), Rcx,
                          Rcy, this->rc(), type_, r, u);
    } break;

    default:
      return Side::Positive;
      break;
  }
}

double Surface::distance(const Vector& r, const Direction& u) const {
  switch (type_) {
    case Type::XPlane:
      return xplane_distance(this->x0(), r, u);
      break;

    case Type::YPlane:
      return yplane_distance(this->y0(), r, u);
      break;

    case Type::Plane:
      return plane_distance(this->A(), this->B(), this->C(), r, u);
      break;

    case Type::Cylinder:
      return cylinder_distance(this->x0(), this->y0(), this->r(), r, u);
      break;

    case Type::BWRCornerI: {
      const double Rcx = this->xh() - this->rc();
      const double Rcy = this->yh() - this->rc();
      return bwr_box_distance(this->xl(), this->xh(), this->yl(), this->yh(),
                              Rcx, Rcy, this->rc(), type_, r, u);
    } break;

    case Type::BWRCornerII: {
      const double Rcx = this->xl() + this->rc();
      const double Rcy = this->yh() - this->rc();
      return bwr_box_distance(this->xl(), this->xh(), this->yl(), this->yh(),
                              Rcx, Rcy, this->rc(), type_, r, u);
    } break;

    case Type::BWRCornerIII: {
      const double Rcx = this->xl() + this->rc();
      const double Rcy = this->yl() + this->rc();
      return bwr_box_distance(this->xl(), this->xh(), this->yl(), this->yh(),
                              Rcx, Rcy, this->rc(), type_, r, u);
    } break;

    case Type::BWRCornerIV: {
      const double Rcx = this->xh() - this->rc();
      const double Rcy = this->yl() + this->rc();
      return bwr_box_distance(this->xl(), this->xh(), this->yl(), this->yh(),
                              Rcx, Rcy, this->rc(), type_, r, u);
    } break;

    default:
      return INF;
      break;
  }
}

double Surface::integrate_x(double xmin, double xmax, const Side side) const {
  double int_out = 0.;
  bool swapped = false;
  if (xmax < xmin) {
    std::swap<double>(xmin, xmax);
    swapped = true;
  }

  switch (type_) {
    case Type::XPlane:
      int_out = xplane_integrate_x(this->x0(), xmin, xmax, side);
      break;

    case Type::YPlane:
      int_out = yplane_integrate_x(this->y0(), xmin, xmax, side);
      break;

    case Type::Plane:
      int_out =
          plane_integrate_x(this->A(), this->B(), this->C(), xmin, xmax, side);
      break;

    case Type::Cylinder:
      int_out = cylinder_integrate_x(this->x0(), this->y0(), this->r(), xmin,
                                     xmax, side);
      break;

    case Type::BWRCornerI:
      int_out = bwr_box_I_integrate_x(this->xl(), this->xh(), this->yl(),
                                      this->yh(), this->rc(), xmin, xmax, side);
      break;

    case Type::BWRCornerII:
      int_out =
          bwr_box_II_integrate_x(this->xl(), this->xh(), this->yl(), this->yh(),
                                 this->rc(), xmin, xmax, side);
      break;

    case Type::BWRCornerIII:
      int_out =
          bwr_box_III_integrate_x(this->xl(), this->xh(), this->yl(),
                                  this->yh(), this->rc(), xmin, xmax, side);
      break;

    case Type::BWRCornerIV:
      int_out =
          bwr_box_IV_integrate_x(this->xl(), this->xh(), this->yl(), this->yh(),
                                 this->rc(), xmin, xmax, side);
      break;
  }

  if (swapped) int_out = -int_out;
  return int_out;
}

double Surface::integrate_y(double ymin, double ymax, const Side side) const {
  double int_out = 0.;
  bool swapped = false;
  if (ymax < ymin) {
    std::swap<double>(ymin, ymax);
    swapped = true;
  }

  switch (type_) {
    case Type::XPlane:
      return xplane_integrate_y(this->x0(), ymin, ymax, side);
      break;

    case Type::YPlane:
      return yplane_integrate_y(this->y0(), ymin, ymax, side);
      break;

    case Type::Plane:
      return plane_integrate_y(this->A(), this->B(), this->C(), ymin, ymax,
                               side);
      break;

    case Type::Cylinder:
      return cylinder_integrate_y(this->x0(), this->y0(), this->r(), ymin, ymax,
                                  side);
      break;

    case Type::BWRCornerI:
      int_out = bwr_box_I_integrate_y(this->xl(), this->xh(), this->yl(),
                                      this->yh(), this->rc(), ymin, ymax, side);
      break;

    case Type::BWRCornerII:
      int_out =
          bwr_box_II_integrate_y(this->xl(), this->xh(), this->yl(), this->yh(),
                                 this->rc(), ymin, ymax, side);
      break;

    case Type::BWRCornerIII:
      int_out =
          bwr_box_III_integrate_y(this->xl(), this->xh(), this->yl(),
                                  this->yh(), this->rc(), ymin, ymax, side);
      break;

    case Type::BWRCornerIV:
      int_out =
          bwr_box_IV_integrate_y(this->xl(), this->xh(), this->yl(), this->yh(),
                                 this->rc(), ymin, ymax, side);
      break;
  }

  if (swapped) int_out = -int_out;
  return int_out;
}

}  // namespace scarabee
