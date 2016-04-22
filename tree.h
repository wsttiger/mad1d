#include <map>
#include <cmath>
#include "Matrix.h"

struct Key {
  int n;
  int l;
  Key() : n(0), l(0) {}
  Key(int n, int l) : n(n), l(l) {}

  bool operator< (const Key &k) const {
    printf("n: %d  l: %d    k.n: %d  k.l: %d\n", n, l, k.n, k.l);
    return ((1<<n)+l < (k.n<<n)+k.l);
  }

  // bool operator< (const Key& key) const {
  //   if ((n < key.n) || (n == key.n && l < key.l)) return true;
  //   return false;
  // }
};

struct MyHashCompare {
    static size_t hash( const Key& x ) {
        size_t prime = 31;
        size_t result = 1;
        result = prime * result + x.l;
        result = prime * result + x.n;
        return result;
    }
    //! True if strings are equal
    static bool equal( const Key& x, const Key& y ) {
        return ((x.n == y.n) && (x.l == y.l));
    }
};

// class Tree {
// private:
//   std::map<Key,real_matrix,MyHashCompare> _map;
// public:
//   Tree() {}
//   ~Tree() {}
// };
