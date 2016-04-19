#include "Matrix.h"
#include "twoscale.h"
#include "gauss_legendre.h"

using Vector = vector<double>;

double normf(Vector v) {
  auto s = 0.0;
  for (auto i = 0; i < v.size(); i++) s+=v[i];
  return s;
}

class Function1D {
private:
  maxlevel = 10;
  Tree tree;
  int k;
  Vector quad_x;
  Vector quad_w;
  int quad_npts;
  real_matrix hg;
  real_matrix hgT;

public:
  Function1D() {}

  ~Function1D() {}

  void init_twoscale(int k) {
    hg = TwoScaleCoeffs::instance()->hg(k); 
    hgT = transpose(hg);
  }

  void init_quadrature(int order) {
    quad_x = gauss_legendre_x(order);
    quad_w = gauss_legendre_w(order);
    quad_npts = w.size();
    real_matrix quad_phi  = zeros(npt, k);
    real_matrix quad_phiT = zeros(k, npt);
    real_matrix quad_phiw = zeros(npt, k);
    auto p = ScalingFunction::instance()->phi(quad_x, k);
    for (auto i = 0; i < npt; i++) {
      for (auto m = 0; m < k; m++) {
        quad_phi(i,m) = p[m];
        quad_phiT(m,i) = p[m];
        quad_phiw[i,m] = quad_w[i]*p[m];
      }
    }
  }

  Vector project_box(double* f, int n, int l) {
    auto s = Vector(k);
    auto h = std::pow(0.5,n);
    auto scale = std::sqrt(h);
    for (auto mu = 0; mu < npt; mu++) {
      auto x = (l + quad_x[mu]) * h;  
      auto fx = f(x);
      for (auto i = 0; i < k; i++) {
        s[i] += scale*fx*quad_phiw(mu,i);
      }
    }
    return s;
  }

  void refine(double* f, int n, int l) {
    Vector s0 = project_box(f, n+1, 2*l);
    Vector s1 = project_box(f, n+1, 2*l+1);
    Vector s(2*k);
    for (auto i = 0; i < k; i++) s[i] = s0[i];
    for (auto i = k; i < 2*k; i++) s[i+k] = s1[i];
    Vector d = hg*s;
    if (normf(d(d.begin()+k,d.begin()+2*k)) < thresh || n >= maxlevel-1) {
      tree[Key(n+1,2*l)] = s0;
      tree[Key(n+1,2*l+1)] = s1;
    } else {
      refine(f, n+1, 2*l);
      refine(f, n+1, 2*l+1);
    }
  }
};
