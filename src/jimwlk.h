#ifndef SRC_JIMWLK_H_
#define SRC_JIMWLK_H_

#include <complex>
#include <memory>
#include <vector>

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

    const double fmgev = 5.068;

    int Ngrid_;
    int Nc2m1_;
    int Ncells_;

    Group *group_ptr_;
    Random *random_ptr_;
    Lattice *lat_ptr_;

    std::vector<std::complex<double> > **K_;
    std::vector<std::complex<double> > **S_;

  public:
    JIMWLK() = delete;
    JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random);
    ~JIMWLK();

    void initializeKandS();
    double getMassRegulator(const double x, const double y) const;
};

#endif  // SRC_JIMWLK_H_
