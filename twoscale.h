#include <string>
#include <fstream>
#include "Matrix.h"

//Return the twoscale coefficients for the multiwavelets of order k.
//
//Note that a cached value is returned ... if you want to modify it
//take a copy first
class TwoScaleCoeffs {
private:
  std::vector<real_matrix> coeffs;
  static TwoScaleCoeffs* _instance;
  TwoScaleCoeffs() {
    std::ifstream fil("tscoeffs");
    if (fil.is_open()) {
      int maxk = -1;
      fil >> maxk;
      for (auto k = 1; k <= maxk; k++) {
        real_matrix matk;
        matk.create(2*k,2*k);
        int readk; double val = 0.0;
        fil >> readk;
        assert(k == readk);
        for (auto j = 0; j < 4*k*k; j++) {
          fil >> val;
          matk[j] = val; 
        }
        coeffs.push_back(matk); 
      }
    }
  }
  ~TwoScaleCoeffs() {}
public:
  static TwoScaleCoeffs* instance() {
    if (!_instance) _instance = new TwoScaleCoeffs();
    return _instance;
  }
  real_matrix hg(int k) {return coeffs[k-1];}
};
TwoScaleCoeffs* TwoScaleCoeffs::_instance = 0;

// Evaluate the Legendre polynomials up to the given order at x
// defined on [-1,1].
std::vector<double> legendre(double x,int order) {
  std::vector<double> p(order+1);
  p[0] = 1.0;
  if (order == 0) return p;
  p[1] = x;
  for (int j = 1; j < order; j++) {
    p[j+1] = j*(x*p[j] - p[j-1])/(j+1) + x*p[j];
  }
  return p;
}

class ScalingFunction {
private:
  double norms[100];
  static ScalingFunction* _instance;
  ScalingFunction() {
    for (int i = 0; i < 100; i++) {
      norms[i] = std::sqrt(2*i+1);
    } 
  }
  ~ScalingFunction() {}
public:
  static ScalingFunction* instance() {
    if (!_instance) _instance = new ScalingFunction();
    return _instance;
  }
  // Evaluate the shifted normalized Legendre polynomials up to the
  // given order at x defined on [0,1].
  // These are also our scaling functions, phi_i(x) , i=0..k-1
  // In addition to forming an orthonormal basis on [0,1] we have
  // phi_j(1/2-x) = (-1)^j phi_j(1/2+x)
  // (the wavelets are similar with phase (-1)^(j+k)).
  std::vector<double> phi(double x,int k) {
    auto order = k-1;
    auto p = legendre(2.*x-1, order);
    for (auto n = 0; n < k; n++)
      p[n] = p[n]*norms[n];
    return p;
  }
};
ScalingFunction* ScalingFunction::_instance = 0;


