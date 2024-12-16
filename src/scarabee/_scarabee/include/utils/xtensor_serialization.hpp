#ifndef SCARABEE_SERIALIZATION_H
#define SCARABEE_SERIALIZATION_H

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>

#include <xtensor/xtensor.hpp>
#include <xtensor/xarray.hpp>

namespace scarabee {

// xarray
template<class Archive, class T>
void save(Archive& arc, const xt::xarray<T>& a) {
  const std::size_t ndims = a.shape().size();
  arc(CEREAL_NVP(ndims));

  for (std::size_t i = 0; i < ndims; i++) {
    arc(a.shape()[i]);
  }

  const std::size_t size = a.size();
  arc(CEREAL_NVP(size));
  for (std::size_t i = 0; i < size; i++) {
    arc(a.flat(i));
  }
}

template<class Archive, class T>
void load(Archive& arc, xt::xarray<T>& a) {
  std::size_t ndims = 0;
  arc(CEREAL_NVP(ndims));

  std::vector<std::size_t> shape(ndims, 0);
  for (std::size_t i = 0; i < ndims; i++) {
    arc(shape[i]);
  }

  a.resize(shape);
  std::size_t size = 0;
  arc(CEREAL_NVP(size));
  for (std::size_t i = 0; i < size; i++) {
    arc(a.flat(i));
  }
}

// 1D Tensor
template<class Archive, class T>
void save(Archive& arc, const xt::xtensor<T, 1>& a) {
  const std::size_t shape_0 = a.shape()[0];
  arc(CEREAL_NVP(shape_0));
  for (std::size_t i = 0; i < shape_0; i++) {
    arc(a.flat(i));
  }
}

template<class Archive, class T>
void load(Archive& arc, xt::xtensor<T, 1>& a) {
  std::size_t shape_0 = 0;
  arc(CEREAL_NVP(shape_0));
  a.resize({shape_0});

  for (std::size_t i = 0; i < shape_0; i++) {
    arc(a.flat(i));
  }
}

// 2D Tensor
template<class Archive, class T>
void save(Archive& arc, const xt::xtensor<T, 2>& a) {
  const std::size_t shape_0 = a.shape()[0];
  const std::size_t shape_1 = a.shape()[1];
  arc(CEREAL_NVP(shape_0));
  arc(CEREAL_NVP(shape_1));
  for (std::size_t i = 0; i < shape_0*shape_1; i++) {
    arc(a.flat(i));
  }
}

template<class Archive, class T>
void load(Archive& arc, xt::xtensor<T, 2>& a) {
  std::size_t shape_0 = 0;
  std::size_t shape_1 = 0;
  arc(CEREAL_NVP(shape_0));
  arc(CEREAL_NVP(shape_1));
  a.resize({shape_0, shape_1});

  for (std::size_t i = 0; i < shape_0*shape_1; i++) {
    arc(a.flat(i));
  }
}

// 3D Tensor
template<class Archive, class T>
void save(Archive& arc, const xt::xtensor<T, 3>& a) {
  const std::size_t shape_0 = a.shape()[0];
  const std::size_t shape_1 = a.shape()[1];
  const std::size_t shape_2 = a.shape()[2];
  arc(CEREAL_NVP(shape_0));
  arc(CEREAL_NVP(shape_1));
  arc(CEREAL_NVP(shape_2));
  for (std::size_t i = 0; i < shape_0*shape_1*shape_2; i++) {
    arc(a.flat(i));
  }
}

template<class Archive, class T>
void load(Archive& arc, xt::xtensor<T, 3>& a) {
  std::size_t shape_0 = 0;
  std::size_t shape_1 = 0;
  std::size_t shape_2 = 0;
  arc(CEREAL_NVP(shape_0));
  arc(CEREAL_NVP(shape_1));
  arc(CEREAL_NVP(shape_2));
  a.resize({shape_0, shape_1, shape_2});

  for (std::size_t i = 0; i < shape_0*shape_1*shape_2; i++) {
    arc(a.flat(i));
  }
}

// svector
template<class Archive, class T>
void save(Archive& arc, const xt::svector<T>& a) {
  const std::size_t size = a.size();
  arc(CEREAL_NVP(size));
  for (std::size_t i = 0; i < size; i++) {
    arc(a[i]);
  }
}

template<class Archive, class T>
void load(Archive& arc, xt::svector<T>& a) {
  std::size_t size = 0;
  arc(CEREAL_NVP(size));
  a.resize(size);
  for (std::size_t i = 0; i < size; i++) {
    arc(a[i]);
  }
}

}

#endif