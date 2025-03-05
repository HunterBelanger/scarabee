#include <moc/surface.hpp>
#include <utils/constants.hpp>

#include <cmath>

namespace scarabee {

//=============================================================================
// Side Methods
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

//=============================================================================
// Distance Methods
//-----------------------------------------------------------------------------

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

}  // namespace scarabee
