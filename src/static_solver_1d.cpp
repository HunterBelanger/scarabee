#include <solver.hpp>

StaticSolver1D::StaticSolver1D(Geom1D&& geom): geometry_(geom), indxr_(geometry_.size(), geometry_.ngroups()), L_(), F_() {
  this->load_operators();
}

void StaticSolver1D::load_operators() {
  // Get the number of elements in the problem solution. This is the number of
  // spatial regions times the number of energy groups.
  const std::size_t elems = geometry_.size()*geometry_.ngroups();

  // Resize fission operator
  F_.resize(elems, elems);

  // Make temp removal, in-scattering, and leakage operators
  Eigen::SparseMatrix<float> R_, S_, Lk_;
  R_.resize(elems, elems);
  S_.resize(elems, elems);
  Lk_.resize(elems, elems);

  // Load fission, removal, and in-scattering operators
  for (std::size_t i = 0; i < geometry_.size(); i++) {
    const MatTile& tile = geometry_(i);

    if (tile.material == nullptr) continue;

    for (std::uint32_t gin = 0; gin < geometry_.ngroups(); gin++) {
      //------------------------------------------------------------------------
      // First, we determine the leakage operators for this node
      const MatTile* tile_im1 = nullptr;
      const MatTile* tile_ip1 = nullptr;
      const float dx = geometry_.dx(i);
      const float alph_i = tile.material->D(gin) / dx;
      float alph_im1 = 0.;
      float alph_ip1 = 0.;

      Lk_.insert(indxr_(i, gin), indxr_(i, gin)) = 0.;
      if (i > 0) {
        Lk_.insert(indxr_(i-1, gin), indxr_(i-1, gin)) = 0.;
        tile_im1 = &geometry_(i-1);
      }
      if (i < (geometry_.size() - 1)) {
        Lk_.insert(indxr_(i+1, gin), indxr_(i+1, gin)) = 0.;
        tile_ip1 = &geometry_(i+1);
      }

      if (tile_im1 && tile_im1->material && tile_ip1 && tile_ip1->material) {
        // This is a normal node, with no B.C.s
        alph_im1 = tile_im1->material->D(gin) / geometry_.dx(i-1);
        alph_ip1 = tile_ip1->material->D(gin) / geometry_.dx(i+1);

        Lk_.coeffRef(indxr_(i,gin), indxr_(i+1,gin)) += (-2.f*alph_i*alph_ip1/(alph_i+alph_ip1)) / dx;
        Lk_.coeffRef(indxr_(i,gin), indxr_(i,gin)) += (2.f*alph_i + (2.f*alph_im1*alph_i/(alph_im1+alph_i)) - (2.f*alph_i*alph_i/(alph_i+alph_ip1))) / dx;
        Lk_.coeffRef(indxr_(i,gin), indxr_(i-1,gin)) += (-2.f*alph_im1 + 2.f*alph_im1*alph_im1/(alph_im1+alph_i)) / dx;
      } else if (tile_im1 && tile_im1->type == MatTile::Type::Reflection) {
        alph_ip1 = tile_ip1->material->D(gin) / geometry_.dx(i+1);
        
        Lk_.coeffRef(indxr_(i,gin), indxr_(i+1,gin)) += (-2.f*alph_i*alph_ip1/(alph_i+alph_ip1)) / dx;
        Lk_.coeffRef(indxr_(i,gin), indxr_(i,gin)) += (2.f*alph_i - (2.f*alph_i*alph_i/(alph_i+alph_ip1))) / dx;
      } else if (tile_im1 && tile_im1->type == MatTile::Type::Vacuum) {
        alph_ip1 = tile_ip1->material->D(gin) / geometry_.dx(i+1);

      } else if (tile_ip1 && tile_ip1->type == MatTile::Type::Reflection) {
        alph_im1 = tile_im1->material->D(gin) / geometry_.dx(i-1);
        
        Lk_.coeffRef(indxr_(i,gin), indxr_(i,gin)) += (2.f*alph_im1*alph_i/(alph_im1+alph_i)) / dx;
        Lk_.coeffRef(indxr_(i,gin), indxr_(i-1,gin)) += (-2.f*alph_im1 + 2.f*alph_im1*alph_im1/(alph_im1+alph_i)) / dx;
      } else if (tile_ip1 && tile_ip1->type == MatTile::Type::Vacuum) {
        alph_im1 = tile_im1->material->D(gin) / geometry_.dx(i-1);

      }
      
      //------------------------------------------------------------------------
      // Get fission yield * fission xs
      const float nuEf = tile.material->nu(gin) * tile.material->Ef(gin);

      // To calculate removal xs, we start with absorption
      float Er = tile.material->Ea(gin);

      for (std::uint32_t gout = 0; gout < geometry_.ngroups(); gout++) {
        const float Es_gin_to_gout = tile.material->Es(gin, gout);
        const float Es_gout_to_gin = tile.material->Es(gout, gin);

        if (gout != gin) {
          // Must add scattering out of our current group to removal xs
          Er += Es_gin_to_gout;

          // Must account for scattering from gout into gin in the
          // in-scattering matrix
          S_.insert(indxr_(i,gout), indxr_(i,gin)) = Es_gout_to_gin;
        }
        
        const float chi_gin_gout = tile.material->chi(gin, gout);
        F_.insert(indxr_(i,gout), indxr_(i,gin)) = nuEf * chi_gin_gout;
      }

      R_.insert(indxr_(i,gin), indxr_(i,gin)) = Er;
    }
  }

}

void StaticSolver1D::solve() {

}