#include "function1d.h"

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825;
const int k = 4;
const double thresh = 1e-8;

double func1(double x) {
  auto c = 1.0;
  auto mu = 2.0;
  return c*sin(2*PI*x)*std::exp(-mu*x);
}

double func2(double x) {
  auto c = 1.0;
  auto mu = 2.0;
  return c*cos(2*PI*x)*std::exp(-mu*x);
}

double sum_f12(double x) {
  return func1(x) + func2(x);
}

void test_add() {
  auto fr = Function1D(func1, k, thresh, 30, 2);
  auto fc = compress(fr);
  printf("\nfunction fc:\n");
  fc.print_tree();
  auto gr = Function1D(func2, k, thresh, 30, 2);
  auto gc = compress(gr);
  printf("\nfunction gc:\n");
  gc.print_tree();
  auto fgc = fc + gc;
  printf("\nfunction f+g(c):\n");
  fgc.print_tree();
  auto f12r = Function1D(sum_f12, k, thresh, 30, 2);
  auto f12c = compress(f12r);
  printf("\nfunction f12c:\n");
  f12c.print_tree();
}

int main(int argc, char** argv) {

//   Function1D f(func, 8, 1e-10, 30, 2);
//   f.print_tree();
//   auto x = 0.23111;
//   printf("x: %15.8e f: %15.8e func: %15.8e error: %15.8e\n", x, f(x), func(x), std::abs(f(x)-func(x)));
//   auto g = compress(f);
//   auto g2 = g + g;
//   g.print_tree();
//   g2.print_tree();
//   auto h = 2.0*reconstruct(g);
//   h.print_tree();

  test_add();
  return 0;
}
