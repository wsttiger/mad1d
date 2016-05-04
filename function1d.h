#include "Matrix.h"
#include "tree.h"
#include "twoscale.h"
#include "gauss_legendre.h"

using Vector = real_vector;
using CoeffTree = std::map<Key,Vector>;

double normf(Vector v) {
  auto s = 0.0;
  //for (auto t : v) s+=t*t;
  for (auto i = 0; i < v.size(); i++) s+=v[i]*v[i];
  return std::sqrt(s);
}

class Function1D {
private:
  enum class FunctionForm {RECONSTRUCTED, COMPRESSED, UNDEFINED};

  FunctionForm form;
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
  friend Function1D operator*(const double& s, const Function1D& f);
  friend Function1D operator*(const Function1D& f, const double& s);
  friend Function1D compress(const Function1D& f);
  friend Function1D reconstruct(const Function1D& f);

  Function1D(int k, double thresh, int maxlevel = 30, int initiallevel = 4) 
   : k(k), thresh(thresh), maxlevel(maxlevel), initiallevel(initiallevel) {
    form = FunctionForm::UNDEFINED;
    init_twoscale(k);
    init_quadrature(k);
  }

  Function1D(double (*f) (double), int k, double thresh, int maxlevel = 30, int initiallevel = 4) 
   : k(k), thresh(thresh), maxlevel(maxlevel), initiallevel(initiallevel) {
    form = FunctionForm::RECONSTRUCTED;
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
    Vector s0 = project_box(f, n+1, 2*l);
    Vector s1 = project_box(f, n+1, 2*l+1);
    Vector s(2*k);
    for (auto i = 0; i < k; i++) {
      s[i] = s0[i];
      s[i+k] = s1[i];
    }
    Vector d = hg*s;
    if (normf(Vector(d.slice(k,2*k-1))) < thresh || n >= maxlevel-1) {
      stree[Key(n+1,2*l)] = s0;
      stree[Key(n+1,2*l+1)] = s1;
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
      reconstruct_spawn(stree_r, s1, n+1, 2*l+1);
    } else {
      stree_r[Key(n,l)] = ss;
    }
  }
  
  Vector compress_spawn(CoeffTree& dtree_r, int n, int l) const {
    auto s0p = stree.find(Key(n+1,2*l));
    auto s1p = stree.find(Key(n+1,2*l+1));
    Vector s0 = (s0p == stree.end()) ? compress_spawn(dtree_r, n+1, 2*l)   : s0p->second;
    Vector s1 = (s1p == stree.end()) ? compress_spawn(dtree_r, n+1, 2*l+1) : s1p->second;
    Vector s(2*k);
    for (auto i = 0; i < k; i++) {
      s[i]   = s0[i];
      s[i+k] = s1[i];
    }
    Vector d = hg*s;
    auto sr = d.slice(0,k-1); 
    auto dr = d.slice(k,2*k-1);
    if (normf(dr) > 1e5) {
      printf("%d  %d     %15.8e  %15.8e  %15.8e  %15.8e\n",n, l, normf(s0), normf(s1), normf(sr), normf(dr)); 
      printf("s0\n");
      print(s0);
      printf("\n");
      printf("s1\n");
      print(s1);
      printf("\n");
      printf("sr\n");
      print(sr);
      printf("\n");
      printf("dr\n");
      print(dr);
      printf("\n");
      printf("s\n");
      print(s);
      printf("\n");
      printf("d\n");
      print(d);
      printf("\n");
    }
    dtree_r[Key(n,l)] = dr; 
    return sr; 
  }

  double operator()(double x) {
    return eval(x, 0, 0);
  }

  Function1D operator+(const Function1D& f) const {
    // Make sure that everybody is compressed
    assert(form == FunctionForm::COMPRESSED);
    assert(f.form == FunctionForm::COMPRESSED);
    Function1D r(f.k, f.thresh, f.maxlevel, f.initiallevel);
    r.form = f.form;
    auto& dtree_r = r.dtree;
    auto& stree_r = r.stree;
    // Loop over d-coeffs in this tree and add these coeffs
    // to the d-coeffs in the f tree IF THEY EXIST 
    // then insert into the result
    for (auto c : dtree) {
      auto key = c.first;
      auto dcoeffs = copy(c.second);
      auto c2 = f.dtree.find(key);
      if (c2 != f.dtree.end()) {
        dcoeffs += c2->second;
        dtree_r[key] = dcoeffs;
      } else {
        dtree_r[key] = dcoeffs;
      }
    }
    // Loop over the remainder d-coeffs in the f tree and insert
    // into the result tree
    for (auto c : f.dtree) {
      auto key = c.first;
      auto c2 = dtree_r.find(key);
      if (c2 == dtree_r.end()) {
        dtree_r[key] = copy(c.second);
      }
    }
    // Do s0 coeffs
    auto c1 = stree.find(Key(0,0));
    auto c2 = f.stree.find(Key(0,0));
    assert(c1 != stree.end());
    assert(c2 != f.stree.end());
    stree_r[Key(0,0)] = c1->second + c2->second;
    return r;
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
};

Function1D operator*(const double& s, const Function1D& f) {
  Function1D r(f.k, f.thresh, f.maxlevel, f.initiallevel);
  r.form = f.form;
  for (auto cp : f.stree) {
    auto key = cp.first;
    auto c = cp.second; 
    auto c2 = copy(c);
    c2.scale(s);
    r.stree[key] = c2;
  }
  for (auto cp : f.dtree) {
    auto key = cp.first;
    auto c = cp.second; 
    auto c2 = copy(c);
    c2.scale(s);
    r.dtree[key] = c2; 
  }
  return r;
}

Function1D operator*(const Function1D& f, const double& s) {
  return s*f;
}

Function1D compress(const Function1D& f) {
  Function1D r(f.k, f.thresh, f.maxlevel, f.initiallevel);
  auto s0 = f.compress_spawn(r.dtree, 0, 0);
  r.stree[Key(0,0)] = s0;
  r.form = Function1D::FunctionForm::COMPRESSED;
  return r;
}

Function1D reconstruct(const Function1D& f) {
  Function1D r(f.k, f.thresh, f.maxlevel, f.initiallevel);
  const auto s0 = f.stree.find(Key(0,0))->second;
  f.reconstruct_spawn(r.stree, s0, 0, 0);
  r.form = Function1D::FunctionForm::RECONSTRUCTED;
  return r;
}


