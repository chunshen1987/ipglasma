#include "jimwlk.h"

JIMWLK::JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random)
    : param_(param) {
    nn_[0] = param_.getSize();
    nn_[1] = param_.getSize();

    fft_ptr_ = std::make_shared<FFT>(nn_);

    Nc2m1_ = param_.getNc() * param_.getNc() - 1;
    Ncells_ = param_.getSize() * param_.getSize();

    group_ptr_ = group;
    random_ptr_ = random;
    lat_ptr_ = lat;
}

JIMWLK::~JIMWLK() {}
