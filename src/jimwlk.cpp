#include "jimwlk.h"

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

JIMWLK::JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random)
    : param_(param) {
    nn_[0] = param_.getSize();
    nn_[1] = param_.getSize();
    Ngrid_ = param_.getSize();

    fft_ptr_ = std::make_shared<FFT>(nn_);

    Nc2m1_ = param_.getNc() * param_.getNc() - 1;
    Ncells_ = param_.getSize() * param_.getSize();

    group_ptr_ = group;
    random_ptr_ = random;
    lat_ptr_ = lat;

    initializeKandS();
}

JIMWLK::~JIMWLK() {
    for (int i = 0; i < Ncells_; i++) {
        delete K_[i];
        delete S_[i];
    }
    delete[] K_;
    delete[] S_;
}

void JIMWLK::initializeKandS() {
    K_ = new std::vector<std::complex<double> > *[Ncells_];
    S_ = new std::vector<std::complex<double> > *[Ncells_];
    for (int i = 0; i < Ncells_; i++) {
        K_[i] = new std::vector<std::complex<double> >;
        S_[i] = new std::vector<std::complex<double> >;
    }

    for (int pos = 0; pos < Ncells_; pos++) {
        double x = pos / Ngrid_ - static_cast<double>(Ngrid_) / 2.;
        double y = pos % Ngrid_ - static_cast<double>(Ngrid_) / 2.;
        double r2 = x * x + y * y;
        if (r2 < 1e-16) {
            K_[pos]->push_back(0.);
            K_[pos]->push_back(0.);
            if (param_.getSimpleLangevin() == false) {
                S_[pos]->push_back(0.);
            }
            continue;
        }
        double mass_regulator = getMassRegulator(x, y);
    }
}

double JIMWLK::getMassRegulator(const double x, const double y) const {
    // if m suppresses long distance tails,
    // K is multiplied by this, which is m*r*K_1(m*r)
    double mass_regulator = 1.0;
    double m = param_.getm();
    if (m < 1e-16) {
        return mass_regulator;
    }
    const double fmgev = 5.068;
    double length = param_.getL();

    // Lattice units
    // Here x is [-N/2, N/2]
    double lat_x = sin(M_PI * x / Ngrid_) / (M_PI);
    double lat_y = sin(M_PI * y / Ngrid_) / (M_PI);
    // lat_x and lat_y are now in [-1/2,1/2] as x/nn[0] is in [-N/2, N/2]
    double lat_r = sqrt(lat_x * lat_x + lat_y * lat_y) * Ngrid_;
    // lat_r now tells how many lattice units the distance is
    double a = length / Ngrid_;
    double lat_m = m * a * fmgev;
    double bessel_argument = lat_m * lat_r;
    double bes = std::cyl_bessel_k(1, bessel_argument);
    mass_regulator = bessel_argument * bes;
    return mass_regulator;
}
