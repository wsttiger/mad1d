#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <cassert>
#include <memory>
#include <complex>

using std::vector;
using std::pair;

struct Slice {
  Slice(int i0, int i1, int j0, int j1) 
   : i0(i0), i1(i1), j0(j0), j1(j1), ndim(2) {}
  
  Slice(int i0, int i1) 
   : i0(i0), i1(i1), j0(0), j1(0), ndim(1) {}

  // size of first dimension
  int n0() const {
    return i1-i0+1;
  }

  // size of second dimension
  int n1() const {
    return j1-j0+1;
  }

  int i0, i1, j0, j1, ndim;
};

void print(const Slice& sl) {
  if (sl.ndim == 1) printf("Slice: i0 = %d  i1 = %d\n", sl.i0, sl.i1);
  if (sl.ndim == 2) printf("Slice: i0 = %d  i1 = %d  j0 = %d  j1 = %d\n", sl.i0, sl.i1, sl.j0, sl.j1);
}

template <typename T>
class MatrixT {
private:
  // the dimensions
  unsigned int _dim0, _dim1;
  // the data pointer
  std::shared_ptr<T> _sp;
  T* _p;
  // is the matrix currently allocated?
  bool _allocated;

  friend MatrixT copy(const MatrixT& t, bool empty = false) {
    MatrixT r;
    if (t._allocated) {
      r._dim0 = t._dim0;
      r._dim1 = t._dim1;


      int sz = r.size();
      r._p = new T[sz];
      r._sp = std::shared_ptr<T>(r._p, [](T *p) {delete[] p;});
      if (empty) {
        for (int i = 0; i < sz; i++) {
          r._p[i] = 0.0;
        }
      }
      else {
        for (int i = 0; i < sz; i++) {
          r._p[i] = t._p[i];
        }
      }
      r._allocated = t._allocated;
    }
    return r;
  }

public:
  MatrixT() 
    :  _dim0(0), _dim1(0), _p(0),_allocated(false) {}
  
  virtual ~MatrixT() {
    _allocated = false;
    _dim0 = 0;
    _dim1 = 0;
  }
  
  MatrixT& operator=(const MatrixT& t) {
    if (this != &t) {
      _dim0 = t._dim0;
      _dim1 = t._dim1;
      _p = t._p;
      _sp = t._sp;
      _allocated = t._allocated;
    }
    return *this;
  }

  MatrixT(const MatrixT& t) {
    (*this) = t;
  }

  void create(int d0, int d1) {
    // dims
    _dim0 = d0; _dim1 = d1;
    // allocation
    _p = new T[d0*d1];
    _sp = std::shared_ptr<T>(_p, [](T *p) {delete[] p;});
    _allocated = true;
  }

  void create(std::function<T (double, double)> f, const vector<double>& x, const vector<double>& y, 
              int d0, int d1) {
    create(d0, d1);
    for (int i = 0; i < d0; i++)
      for (int j = 0; j < d1; j++)
        _p[i*d1+j] = f(x[i], y[j]);
  }

  int size() const {
    return _dim0*_dim1;
  }

  int nrows() const {
    return _dim0; 
  }

  int ncols() const {
    return _dim1; 
  }

  // set all elements to zero
  void empty() {
    int sz = this->size();
    for (int i = 0; i < sz; i++) {
      _p[i] = T(0);
    }
  }

  // set all elements to value
  void value(T val) {
    int sz = this->size();
    for (int i = 0; i < sz; i++) {
      _p[i] = val;
    }
  }

  // pointer at a base position
  T* ptr() {
    return _p;
  }

  // pointer at a base position
  const T* ptr() const {
    return _p;
  }

  // pointer at a given (i,j) index
  T* ptr(int i0, int i1) {
    return &_p[i0*_dim1+i1];
  }

  // pointer at a given (i,j) index
  const T* ptr(int i0, int i1) const {
    return &_p[i0*_dim1+i0];
  }

  // index operator
  T& operator[](int i) {
    return _p[i];
  }

  // index operator
  T& operator[](int i) const {
    return _p[i];
  }

  // index operator
  T& operator()(int i) {
    return _p[i];
  }

  // index operator
  T& operator()(int i) const {
    return _p[i];
  }

  // index operator
  T& operator()(int i0, int i1) {
    return _p[i0*_dim1+i1];
  }

  // index operator
  T& operator()(int i0, int i1) const {
    return _p[i0*_dim1+i1];
  }

  // perform inplace a*(*this) + b*t
  void gaxpy(const T& a, const MatrixT& t, const T& b) {
    int sz1 = this->size();
    int sz2 = t.size();
    assert(sz1 == sz2);
    for (int i = 0; i < sz1; i++) _p[i] = a*_p[i]+b*t._p[i];
  }

  // perform out-of-place a*(*this) + b*t
  MatrixT gaxpy_oop(const T& a, const MatrixT& t, const T& b) const {
    int sz1 = this->size();
    int sz2 = t.size();
    assert(sz1 == sz2);
    MatrixT r = copy(*this, true);
    for (int i = 0; i < sz1; i++) r._p[i] = a*_p[i]+b*t._p[i];
    return r;
  }

  // addition
  MatrixT operator+(const MatrixT& t) const {
    return gaxpy_oop(1.0, t, 1.0);
  }

  // subtraction
  MatrixT operator-(const MatrixT& t) const {
    return gaxpy_oop(1.0, t, -1.0);
  }

  // inplace scaling by a constant 
  void scale(T a) {
    int sz = this->size();
    for (int i = 0; i < sz; i++) _p[i] *= a;
  }

  // inner like a vector 
  T inner(const MatrixT& t) const {
    int sz1 = this->size();
    int sz2 = t.size();
    assert(sz1 == sz2);
    T rval = 0.0;
    for (int i = 0; i < sz1; i++) rval += _p[i]*t._p[i]; 
    return rval;
  }

  // two norm of the vector
  T norm2() const {
    int sz = this->size();
    T rval = 0.0;
    for (int i = 0; i < sz; i++) rval += _p[i]*_p[i]; 
    return std::sqrt(rval);
  }

  // normalize the values according to the two-norm of the vector
  void normalize() {
    T s = this->norm2();
    this->scale(1./s); 
  }

  // simple slicing
  MatrixT<T> operator()(const Slice& sl) const {
    MatrixT<T> R;
    R.create(sl.n0(),sl.n1());
    const T* p = ptr(sl.i0, sl.j0);
    for (int i = 0; i < sl.n0(); i++) {
      for (int j = 0; j < sl.n1(); j++) {
        R._p[i*R._dim1+j] = p[i*_dim1+j];
      }
    }
    return R;
  }

  // conversion from complex to real
  operator MatrixT<std::complex<T> > () const {
    MatrixT<std::complex<T> > Ac;
    Ac.create(nrows(), ncols());
    int sz = size();
    for (int i = 0; i < sz; i++) Ac(i) = std::complex<T>(_p[i],T(0));
    return Ac;
  }

  // extract a column vector from matrix
  MatrixT<T> col(int j) const {
    MatrixT<T> R;
    int nr = nrows();
    R.create(nr,1);
    for (int i = 0; i < nr; i++) {
      R(i) = _p[i*_dim0+j];
    }
    return R;
  }

  // extract a contiguous group of column vectors from matrix
  std::vector<MatrixT<T> > cols(const Slice& sl) const {
    assert(sl.ndim == 1);
    int nc = sl.n0();
    std::vector<MatrixT<T> > R(nc);
    for (int j = 0; j < nc; j++) {
      R[j] = col(j+sl.i0);
    }
    return R;
  }

  // extract all column vectors from matrix
  std::vector<MatrixT<T> > cols() const {
    int nc = ncols();
    std::vector<MatrixT<T> > R(nc);
    for (int j = 0; j < nc; j++) {
      R[j] = col(j);
    }
    return R;
  }
};

// typedefs
typedef MatrixT<double> real_matrix;
typedef MatrixT<std::complex<double> > complex_matrix;

template <typename Q>
Q sum(const MatrixT<Q>& A) {
  Q s = Q(0);
  for (int i = 0; i < A.size(); i++) s += A[i];
  return s;
}

template <typename Q>
Q product(const MatrixT<Q>& A) {
  Q s = Q(1);
  for (int i = 0; i < A.size(); i++) s *= A[i];
  return s;
}

template <typename Q>
MatrixT<Q> gaxpy(const Q& a, const MatrixT<Q>& T1, const Q& b, const MatrixT<Q>& T2) {
  return T1.gaxpy_oop(a, T2, b);
}

template <typename Q>
double inner(const MatrixT<Q>& t1, const MatrixT<Q>& t2) {
  return t1.inner(t2);
}

double norm2(const real_matrix& t) {
  return t.norm2();
}

double norm2(const complex_matrix& t) {
  return std::abs(t.norm2());
}

void print(const real_matrix& A) {
  int nr = A.nrows();
  int nc = A.ncols();
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      printf("%15.8e  ", A(i,j));
      //printf("%20.25e  ", A(i,j));
    }
    printf("\n");
  }
  printf("\n");
}

void print(const complex_matrix& A) {
  int nr = A.nrows();
  int nc = A.ncols();
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      printf("(%15.8e, %15.8e)  ", real(A(i,j)), imag(A(i,j)));
    }
    printf("\n");
  }
  printf("\n");
}

template <typename Q>
MatrixT<Q> operator*(const Q& s, const MatrixT<Q>& A) {
  MatrixT<Q> r = copy(A);
  r.scale(s);
  return r;
}

template <typename Q>
MatrixT<std::complex<Q> > operator*(const Q& s, const MatrixT<std::complex<Q> >& A) {
  MatrixT<std::complex<Q> > r = copy(A);
  r.scale(std::complex<Q>(s));
  return r;
}

// create zeros matrix
template <typename Q>
MatrixT<Q> zeros(int d0, int d1) {
  MatrixT<Q> r;
  r.create(d0, d1);
  r.empty();  
  return r;
}

// create an constant matrix
template <typename Q>
MatrixT<Q> constant(int d0, int d1, double val) {
  MatrixT<Q> r;
  r.create(d0, d1);
  r.value(val);  
  return r;
}

// create an constant matrix
template <typename Q>
MatrixT<Q> ones(int d0, int d1) {
  return constant(d0, d1, T(1.0));
}

template <typename Q>
MatrixT<Q> from_vector(int d0, int d1, const std::vector<Q>& v) {
  MatrixT<Q> r;
  r.create(d0, d1);
  for (int i = 0; i < d0; i++)
    for (int j = 0; j < d1; j++)
      r(i,j) = v[i*d1+j];
  return r;
}

template <typename Q> 
MatrixT<Q> real(const MatrixT<std::complex<Q> > Ac) {
  MatrixT<Q> Ar;
  Ar.create(Ac.nrows(), Ac.ncols()); 
  int sz = Ar.size();
  for (int i = 0; i < sz; i++) Ar(i) = std::real(Ac(i));
  return Ar;
}

template <typename Q> 
MatrixT<Q> imag(const MatrixT<std::complex<Q> > Ac) {
  MatrixT<Q> Ar;
  Ar.create(Ac.nrows(), Ac.ncols()); 
  int sz = Ar.size();
  for (int i = 0; i < sz; i++) Ar(i) = std::imag(Ac(i));
  return Ar;
}

template <typename Q> 
MatrixT<Q> abs(const MatrixT<std::complex<Q> > Ac) {
  MatrixT<Q> Ar;
  Ar.create(Ac.nrows(), Ac.ncols()); 
  int sz = Ar.size();
  for (int i = 0; i < sz; i++) Ar(i) = std::abs(Ac(i));
  return Ar;
}

real_matrix real(const real_matrix& A) {
  real_matrix Ar;
  Ar.create(A.nrows(), A.ncols()); 
  int sz = Ar.size();
  for (int i = 0; i < sz; i++) Ar(i) = A(i);
  return Ar;
}

real_matrix imag(const real_matrix& A) {
  real_matrix Ar;
  Ar.create(A.nrows(), A.ncols()); 
  int sz = Ar.size();
  for (int i = 0; i < sz; i++) Ar(i) = 0.0;
  return Ar;
}

real_matrix abs(const real_matrix& A) {
  real_matrix Ar;
  Ar.create(A.nrows(), A.ncols()); 
  int sz = Ar.size();
  for (int i = 0; i < sz; i++) Ar(i) = std::abs(A(i));
  return Ar;
}

template <typename Q>
MatrixT<Q> transpose(const MatrixT<Q> A) {
  int nr = A.nrows();
  int nc = A.ncols();
  MatrixT<Q> r = zeros<Q>(nc, nr);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      r(j, i) = A(i, j);
    }
  }
  return r;
}

real_matrix ctranspose(const real_matrix A) {
  return transpose(A);
}

complex_matrix ctranspose(const complex_matrix A) {
  int nr = A.nrows();
  int nc = A.ncols();
  complex_matrix r = zeros<std::complex<double> >(nc, nr);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      r(j, i) = std::conj(A(i, j));
    }
  }
  return r;
}

template<typename Q>
MatrixT<std::complex<Q> > make_complex(const MatrixT<Q>& rp, const MatrixT<Q>& ip) {
  int nr1 = rp.nrows();
  int nc1 = rp.ncols();
  int nr2 = ip.nrows();
  int nc2 = ip.ncols();
  assert(nr1 == nr2);
  assert(nc1 == nc2);
  MatrixT<std::complex<Q> > r = zeros<std::complex<Q> >(nr1, nc2);
  for (int i = 0; i < nr1; i++) {
    for (int j = 0; j < nc1; j++) {
      r(i,j) = std::complex<Q>(rp(i,j),ip(i,j)); 
    }
  }
  return r;
}

#endif
