#include <moc/surface.hpp>
#include <utils/constants.hpp>

#include <cmath>

namespace scarabee {

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

    default:
      return Side::Positive;
      break;
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

    default:
      return INF;
      break;
  }
}

}  // namespace scarabee
