#ifndef SRC_JIMWLK_H_
#define SRC_JIMWLK_H_

#include <memory>

#include "FFT.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "Random.h"

class JIMWLK {
  private:
    Parameters &param_;
    std::shared_ptr<FFT> fft_ptr_;
    int nn_[2];

    int Nc2m1_;
    int Ncells_;

    Group *group_ptr_;
    Random *random_ptr_;
    Lattice *lat_ptr_;

  public:
    JIMWLK() = delete;
    JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random);
    ~JIMWLK();
};

#endif  // SRC_JIMWLK_H_
