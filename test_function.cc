#include "function1d.h"

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825;

double func(double x) {
  auto c = 1.0;
  auto mu = 2.0;
  return c*sin(2*PI*x)*std::exp(-mu*x);
}

int main(int argc, char** argv) {

  Function1D f(func, 8, 1e-10, 30, 2);
  f.print_tree();
  auto x = 0.23111;
  printf("x: %15.8e f: %15.8e func: %15.8e error: %15.8e\n", x, f(x), func(x), std::abs(f(x)-func(x)));
  auto g = compress(f);
  g.print_tree();
  auto h = reconstruct(g);
  h.print_tree();
  return 0;
}
