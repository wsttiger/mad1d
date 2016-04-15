#include "Matrix.h"
#include "gauss_legendre.h"
#include "twoscale.h"

int main(int argc, char** argv) {
  real_matrix hg = TwoScaleCoeffs::instance()->hg(2);
  print(hg);
  return 0;
}
