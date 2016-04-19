#include <map>
#include "Matrix.h"

struct Key {
  int n;
  int l;
  Key(int n, int l) : n(n), l(l) {}
  bool operator< (const Key& key) const {
    if ((n < key.n) || (n == key.n && l < key.l)) return true;
    return false;
  }
};

class Tree {
private:
  std::map<Key,real_matrix> _map;
public:
  Tree() {}
  ~Tree() {}
};
