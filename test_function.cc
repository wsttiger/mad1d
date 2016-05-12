#include "function1d.h"

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825;
const int k = 8;
const double thresh = 1e-4;
const int initiallevel = 1;

double linear(double x) {
  return x;
}

double quadratic(double x) {
  return x*x;
}

double func1(double x) {
  auto c = 1.0;
  auto mu = 2.0;
  return c*sin(2*PI*x)*std::exp(-mu*x);
}

double dfunc1(double x) {
  auto c = 1.0;
  auto mu = 2.0;
  return c*cos(2*PI*x)*std::exp(-mu*x)-mu*c*sin(2*PI*x)*std::exp(-mu*x);
}

double func2(double x) {
  auto c = 1.0;
  auto mu = 2.0;
  return c*cos(2*PI*x)*std::exp(-mu*x);
}

double sum_f12(double x) {
  return func1(x) + func2(x);
}

double mul_f12(double x) {
  return func1(x) * func2(x);
}

void test_function_point() {
   Function1D f(func1, k, thresh, 30, initiallevel);
   f.print_tree();
   auto x = 0.23111;
   printf("x: %15.8e f: %15.8e func: %15.8e error: %15.8e\n", x, f(x), func1(x), std::abs(f(x)-func1(x)));
}

void test_diff_function_point() {
   Function1D f(func1, k, thresh, 30, initiallevel);
   f.print_tree();
   Function1D df = f.diff();
   df.print_tree();
   auto x = 0.23111;
   printf("x: %15.8e f: %15.8e func: %15.8e error: %15.8e\n", x, df(x), dfunc1(x), std::abs(df(x)-dfunc1(x)));
}
void test_function_compress() {
   Function1D f(func1, k, thresh, 30, initiallevel);
   printf("f : \n");
   f.print_tree();
   Function1D g = reconstruct(compress(f));
   //Function1D g = compress(f);
   printf("\ng : \n");
   g.print_tree();
}

void test_add() {
  auto fr = Function1D(func1, k, thresh, 30, initiallevel);
  auto fc = compress(fr);
  printf("\nfunction fc:\n");
  fr.print_tree();
  auto gr = Function1D(func2, k, thresh, 30, initiallevel);
  auto gc = compress(gr);
  printf("\nfunction gc:\n");
  gr.print_tree();
  auto fgc = fc + gc;
  printf("\nfunction f+g(c):\n");
  fgc.print_tree();
  auto f12r = Function1D(sum_f12, k, thresh, 30, initiallevel);
  auto f12c = compress(f12r);
  printf("\nfunction f12c:\n");
  f12c.print_tree();
}

void test_mul() {
  auto fr = Function1D(func1, k, thresh, 30, initiallevel);
  printf("\nfunction fr:\n");
  fr.print_tree();
  auto gr = Function1D(func2, k, thresh, 30, initiallevel);
  printf("\nfunction gr:\n");
  gr.print_tree();
  auto fgr = fr * gr;
  printf("\nfunction f+g(r):\n");
  fgr.print_tree();
  auto f12r = Function1D(mul_f12, k, thresh, 30, initiallevel+1);
  printf("\nfunction f12r:\n");
  f12r.print_tree();

  auto x = 0.23111;
  printf("x: %15.8e f: %15.8e func: %15.8e error: %15.8e\n", x, fgr(x), f12r(x), std::abs(fgr(x)-f12r(x)));
}

void test_simple_mul() {
  auto fr = Function1D(linear, k, thresh, 30, 0);
  printf("\nfunction fr:\n");
  fr.print_tree(true);
  auto gr = Function1D(linear, k, thresh, 30, 0);
  printf("\nfunction gr:\n");
  gr.print_tree(true);
  auto fgr = fr * gr;
  printf("\nfunction f+g(r):\n");
  fgr.print_tree(true);
  auto f12r = Function1D(quadratic, k, thresh, 30, 1);
  printf("\nfunction f12r:\n");
  f12r.print_tree(true);

  auto x = 0.23111;
  printf("x: %15.8e f: %15.8e func: %15.8e error: %15.8e\n", x, fgr(x), f12r(x), std::abs(fgr(x)-f12r(x)));
}

int main(int argc, char** argv) {
  //test_function_compress();
  //test_add();
  //test_mul();
  test_diff_function_point();
  return 0;
}
