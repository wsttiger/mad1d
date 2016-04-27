#include "Matrix.h"
#include "tree.h"
#include "twoscale.h"
#include "gauss_legendre.h"

using Vector = real_vector;
using CoeffTree = std::map<Key,Vector>;

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
  CoeffTree stree;
  CoeffTree dtree;
  Vector quad_x;
  Vector quad_w;
  int quad_npts;
  real_matrix hg;
  real_matrix hgT;
  real_matrix quad_phi;
  real_matrix quad_phiT;
  real_matrix quad_phiw;

public:
  friend Function1D compress(const Function1D& f);
  friend Function1D reconstruct(const Function1D& f);

  Function1D(int k, double thresh, int maxlevel = 30, int initiallevel = 4) 
   : k(k), thresh(thresh), maxlevel(maxlevel), initiallevel(initiallevel) {
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
      stree[Key(n+1,2*l)] = s0;
      if (debug) printf("set n+1 2*l coeff (%d  %d)\n",n+1,2*l);
      stree[Key(n+1,2*l+1)] = s1;
      if (debug) printf("set n+1 2*l+1 coeff (%d  %d)\n",n+1,2*l+1);
    } else {
      refine(f, n+1, 2*l);
      refine(f, n+1, 2*l+1);
    }
  }

  void reconstruct_spawn(CoeffTree& stree_r, const Vector& ss, int n, int l) const {
    auto dp = dtree.find(Key(n,l));  
    if (dp != dtree.end()) {
      Vector dd = dp->second;
      Vector d(2*k);
      for (auto i = 0; i < k; i++) {
        d[i]   = ss[i];
        d[i+k] = dd[i];
      }
      auto s = hgT*d;
      Vector s0(k);
      Vector s1(k);
      for (auto i = 0; i < k; i++) {
        s0[i] = s[i];
        s1[i] = s[i+k];
      }
      reconstruct_spawn(stree_r, s0, n+1, 2*l);
      reconstruct_spawn(stree_r, s0, n+1, 2*l+1);
    } else {
      stree_r[Key(n,l)] = ss;
    }
  }
  
  Vector compress_spawn(CoeffTree& dtree_r, int n, int l) const {
    auto s0p = stree.find(Key(n+1,2*l));
    auto s1p = stree.find(Key(n+1,2*l+1));
    Vector s0, s1;
    if (s0p == stree.end()) {
      s0 = compress_spawn(dtree_r, n+1, 2*l);
    } else {
      s0 = s0p->second;
    }
    if (s1p == stree.end()) {
      s1 = compress_spawn(dtree_r, n+1, 2*l+1);
    } else {
      s1 = s1p->second;
    }
    Vector s(2*k);
    for (auto i = 0; i < k; i++) {
      s[i]   = s0[i];
      s[i+k] = s1[i];
    }
    Vector d = hg*s;
    auto sr = d.slice(0,k); 
    auto dr = d.slice(k,2*k);
    dtree_r[Key(n,l)] = dr; 
    return sr; 
  }

  double operator()(double x) {
    return eval(x, 0, 0);
  }

  double eval(double x, int n, int l) {
    assert(n < maxlevel);
    auto treep = stree.find(Key(n,l));
    if (treep != stree.end()) {
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
    printf("sum coeffs:\n");
    auto s = stree[Key(n,l)];
    printf("[%d, %d] (", n, l);
    for (auto v : s) {
      printf("%8.4f  ", v);
    }
    printf(")  %15.8e\n",normf(s));
    printf("diff coeffs:\n");
    auto d = dtree[Key(n,l)];
    printf("[%d, %d] (", n, l);
    for (auto v : d) {
      printf("%8.4f  ", v);
    }
    printf(")  %15.8e\n",normf(d));
  }

  void print_tree() {
    printf("sum coeffs:\n");
    for (auto c : stree) {
      auto k = c.first;
      auto s = c.second; 
      printf("[%d  %d]     %15.8e\n", k.n, k.l, normf(s));
    }
    printf("diff coeffs:\n");
    for (auto c : dtree) {
      auto k = c.first;
      auto d = c.second; 
      printf("[%d  %d]     %15.8e\n", k.n, k.l, normf(d));
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

// Function1D compress(const Function1D& f) {
//   Function1D r(f.k, f.thresh, f.maxlevel, f.initiallevel);
//   int sz = f.stree.size();
//   assert(sz % 2 == 0);
//   sz = sz / 2;
//   bool keepgoing = true;
//   int cnt = 0;
//   for (auto n = 0; keepgoing && n < r.maxlevel; n++) {
//     auto maxl = 1<<n;
//     for (auto l = 0; keepgoing && l < maxl; l++) {
//       auto s0p = f.stree.find(Key(n+1,2*l));
//       auto s1p = f.stree.find(Key(n+1,2*l+1));
//       if ((s0p != f.stree.end()) && (s0p != f.stree.end())) {
//         auto s0 = s0p->second;
//         auto s1 = s1p->second;
//         Vector s(2*r.k);
//         for (auto i = 0; i < r.k; i++) {
//           s[i] = s0[i];
//           s[i+r.k] = s1[i];
//         }
//         Vector d = r.hg*s;
//         r.stree[Key(n,l)] = d.slice(0,r.k);
//         r.dtree[Key(n,l)] = d.slice(r.k,2*r.k);
//         keepgoing = sz > ++cnt;
//       }
//     }
//   }  
//   return r;
// }

Function1D compress(const Function1D& f) {
  Function1D r(f.k, f.thresh, f.maxlevel, f.initiallevel);
  auto s0 = f.compress_spawn(r.dtree, 0, 0);
  r.stree[Key(0,0)] = s0;
  return r;
}

Function1D reconstruct(const Function1D& f) {
  Function1D r(f.k, f.thresh, f.maxlevel, f.initiallevel);
  const auto s0 = f.stree.find(Key(0,0))->second;
  f.reconstruct_spawn(r.stree, s0, 0, 0);
  return r;
}





