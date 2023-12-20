#ifndef SCARABEE_ONE_OVER_E_H
#define SCARABEE_ONE_OVER_E_H

class OneOverE {
  public:
    OneOverE(double C = 1.): C_{C} {}

    double operator()(double E) const { return C_ / E; }

    double constant() const {
      return C_;
    }

    void set_constant(double C) {
      C_ = C;
    }

  private:
    double C_;
};

#endif