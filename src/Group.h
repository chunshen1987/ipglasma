#ifndef Group_h
#define Group_h

#include <vector>

#include "Matrix.h"

class Group {
  private:
    std::vector<Matrix> t;   // generators of the group
    std::vector<Matrix> tA;  // adjoint representation of generators of the group
    int Nc;                  // number of colors

  public:
    Group(int N);

    Matrix &getT(int i) const { return const_cast<Matrix &>(t[i]); };
    Matrix &getTA(int i) const { return const_cast<Matrix &>(tA[i]); };
};
#endif
