#include "Matrix.h"
#include "tree.h"
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
  int k = 8;
  double thresh = 1e-6;
  int maxlevel = 30;
  std::map<Key,Vector> tree;
  Vector quad_x;
  Vector quad_w;
  int quad_npts;
  real_matrix hg;
  real_matrix hgT;
  real_matrix quad_phi;
  real_matrix quad_phiT;
  real_matrix quad_phiw;

public:
  Function1D(int k, double thresh, int maxlevel = 30) 
   : k(k), thresh(thresh) {
    init_twoscale(k);
    init_quadrature(k);
  }

  ~Function1D() {}

  void init_twoscale(int k) {
    hg = TwoScaleCoeffs::instance()->hg(k); 
    hgT = transpose(hg);
  }

  void init_quadrature(int order) {
    quad_x = gauss_legendre_x(order);
    quad_w = gauss_legendre_w(order);
    quad_npts = quad_w.size();
    quad_phi  = zeros<double>(quad_npts, k);
    quad_phiT = zeros<double>(k, quad_npts);
    quad_phiw = zeros<double>(quad_npts, k);
    for (auto i = 0; i < quad_npts; i++) {
      auto p = ScalingFunction::instance()->phi(quad_x[i], k);
      for (auto m = 0; m < k; m++) {
        quad_phi(i,m) = p[m];
        quad_phiT(m,i) = p[m];
        quad_phiw[i,m] = quad_w[i]*p[m];
      }
    }
  }

  Vector project_box(double (*f)(double), int n, int l) {
    auto s = Vector(k);
    auto h = std::pow(0.5,n);
    auto scale = std::sqrt(h);
    for (auto mu = 0; mu < quad_npts; mu++) {
      auto x = (l + quad_x[mu]) * h;  
      auto fx = f(x);
      for (auto i = 0; i < k; i++) {
        s[i] += scale*fx*quad_phiw(mu,i);
      }
    }
    return s;
  }

  void refine(double (*f)(double), int n, int l) {
    printf("refine n: %d     l: %d\n", n, l);
    Vector s0 = project_box(f, n+1, 2*l);
    Vector s1 = project_box(f, n+1, 2*l+1);
    Vector s(2*k);
    for (auto i = 0; i < k; i++) s[i] = s0[i];
    for (auto i = k; i < 2*k; i++) s[i+k] = s1[i];
    Vector d = hg*s;
    if (normf(Vector(d.begin()+k,d.begin()+2*k)) < thresh || n >= maxlevel-1) {
      tree[Key(n+1,2*l)] = s0;
      tree[Key(n+1,2*l+1)] = s1;
    } else {
      refine(f, n+1, 2*l);
      refine(f, n+1, 2*l+1);
    }
  }
};
