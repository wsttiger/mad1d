#include "Matrix.h"
#include "tree.h"
#include "twoscale.h"
#include "gauss_legendre.h"

using Vector = real_vector;

double normf(Vector v) {
  auto s = 0.0;
  for (auto t : v) s+=t*t;
  return std::sqrt(s);
}

class Function1D {
private:
  bool debug = false;
  int k = 8;
  double thresh = 1e-6;
  int maxlevel = 30;
  int initiallevel = 4;
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
   : k(k), thresh(thresh), maxlevel(maxlevel) {
    init_twoscale(k);
    init_quadrature(k);
  }

  Function1D(double (*f) (double), int k, double thresh, int maxlevel = 30, int initiallevel = 4) 
   : k(k), thresh(thresh), maxlevel(maxlevel), initiallevel(initiallevel) {
    init_twoscale(k);
    init_quadrature(k);
    int ntrans = std::pow(2, initiallevel);
    for (auto l = 0; l < ntrans; l++) refine(f, initiallevel, l);
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
        quad_phiw(i,m) = quad_w[i]*p[m];
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
    if (debug) printf("\nrefine---> n: %d     l: %d\n", n, l);
    Vector s0 = project_box(f, n+1, 2*l);
    if (debug) printf("  computed s0\n");
    Vector s1 = project_box(f, n+1, 2*l+1);
    if (debug) printf("  computed s1\n");
    Vector s(2*k);
    if (debug) printf("  copied s0 and s1 to s\n");
    for (auto i = 0; i < k; i++) {
      s[i] = s0[i];
      s[i+k] = s1[i];
    }
    Vector d = hg*s;
    if (debug) {
      printf("s0: ");
      print(s0);
      printf("s1: ");
      print(s1);
      printf("s: ");
      print(s);
      printf("d: ");
      print(d);
      printf("d.slice(k,2*k-1): ");
      print(Vector(d.slice(k,2*k-1)));
      printf("d slice norm: %15.8e\n", normf(Vector(d.slice(k,2*k-1))));
    }
    if (debug) printf("\n  d = hg*s\n");
    if (normf(Vector(d.slice(k,2*k-1))) < thresh || n >= maxlevel-1) {
      tree[Key(n+1,2*l)] = s0;
      if (debug) printf("set n+1 2*l coeff (%d  %d)\n",n+1,2*l);
      tree[Key(n+1,2*l+1)] = s1;
      if (debug) printf("set n+1 2*l+1 coeff (%d  %d)\n",n+1,2*l+1);
    } else {
      refine(f, n+1, 2*l);
      refine(f, n+1, 2*l+1);
    }
  }

  double operator()(double x) {
    return eval(x, 0, 0);
  }

  double eval(double x, int n, int l) {
    assert(n < maxlevel);
    auto treep = tree.find(Key(n,l));
    if (treep != tree.end()) {
      auto p = ScalingFunction::instance()->phi(x, k);
      auto t = inner(treep->second,p)*std::sqrt(std::pow(2.0,n));
      return t;
    } else {
      auto n2 = n + 1;
      auto l2 = 2*l;
      auto x2 = 2*x; 
      if (x2 >= 1.0) {
        l2 = l2 + 1;
        x2 = x2 - 1;
      }
      return eval(x2, n2, l2);
    } 
  }

  void print_coeffs(int n, int l) {
    auto s = tree[Key(n,l)];
    printf("[%d, %d] (", n, l);
    for (auto v : s) {
      printf("%8.4f  ", v);
    }
    printf(")  %15.8e\n",normf(s));
  }

  void print_tree() {
    for (auto c : tree) {
      auto k = c.first;
      auto s = c.second; 
      printf("[%d  %d]     %15.8e\n", k.n, k.l, normf(s));
    }
  }

  // void summarize() {
  //   printf("sum coeffs:\n");
  //   for (auto c : tree) {
  //     auto key = c.first;
  //     auto s = c.second;
  //   }
  // }
};
