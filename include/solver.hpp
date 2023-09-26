#ifndef SCARABEE_SOLVER_H
#define SCARABEE_SOLVER_H

#include <geometry.hpp>

#include <Eigen/SparseLU>

class Indexer1D {
  public:
    Indexer1D(std::size_t nbins, std::size_t ngroups): nbins_(nbins), ngroups_(ngroups) {}

    std::size_t operator()(std::size_t i, std::size_t g) const {
      return g*nbins_ + i;
    }

  private:
    std::size_t nbins_;
    std::size_t ngroups_;
};

class Solver {
  public:
    virtual ~Solver() = default;
    virtual void solve() = 0;
};

class StaticSolver1D : public Solver {
  public:
    StaticSolver1D(Geom1D&& geom);

    void solve() override final;

  private:
    Geom1D geometry_; // Geometry for the problem
    Indexer1D indxr_; // Indexer for a 1D problem
    Eigen::SparseMatrix<float> L_; // Loses operator
    Eigen::SparseMatrix<float> F_; // Fission operator

    void load_operators();
};

#endif