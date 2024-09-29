#include "jimwlk.h"

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

JIMWLK::JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random)
    : param_(param),
      Nc_(param.getNc()),
      Nc2m1_(param.getNc() * param.getNc() - 1),
      Ngrid_(param.getSize()),
      Ncells_(param.getSize() * param.getSize()) {
    nn_[0] = param_.getSize();
    nn_[1] = param_.getSize();

    fft_ptr_ = std::make_shared<FFT>(nn_);

    group_ptr_ = group;
    random_ptr_ = random;
    lat_ptr_ = lat;

    initializedKandS_ = initializeKandS();
}

JIMWLK::~JIMWLK() {
    if (initializedKandS_) {
        for (int i = 0; i < Ncells_; i++) {
            delete K_[i];
            delete S_[i];
        }
        delete[] K_;
        delete[] S_;
    }
}

bool JIMWLK::initializeKandS() {
    K_ = new std::vector<std::complex<double> > *[Ncells_];
    S_ = new std::vector<std::complex<double> > *[Ncells_];
    for (int i = 0; i < Ncells_; i++) {
        K_[i] = new std::vector<std::complex<double> >;
        S_[i] = new std::vector<std::complex<double> >;
    }

    double mu0 = param_.getMu0();
    double Lambda2 = param_.getLambdaQCD() * param_.getLambdaQCD();

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
        x /= Ngrid_;
        y /= Ngrid_;
        if (param_.getRunningCoupling() == 0) {
            // discretization without singularities
            double tmpk1 = cos(M_PI * y) * (sin(2. * M_PI * x) / (2. * M_PI))
                           / ((pow(sin(M_PI * x) / M_PI, 2.)
                               + pow(sin(M_PI * y) / M_PI, 2.)))
                           / Ngrid_;
            double tmpk2 = cos(M_PI * x) * (sin(2. * M_PI * y) / (2. * M_PI))
                           / ((pow(sin(M_PI * x) / M_PI, 2.)
                               + pow(sin(M_PI * y) / M_PI, 2.)))
                           / Ngrid_;

            // Regulate long distance tails, does nothing if m=0
            tmpk1 *= mass_regulator;
            tmpk2 *= mass_regulator;

            K_[pos]->push_back(tmpk1);
            K_[pos]->push_back(tmpk2);
            S_[pos]->push_back(
                (pow(cos(M_PI * y), 2.)
                     * pow(sin(2. * M_PI * x) / (2. * M_PI), 2.)
                 + pow(cos(M_PI * x), 2.)
                       * pow(sin(2. * M_PI * y) / (2. * M_PI), 2.))
                / pow(
                    (pow(sin(M_PI * x) / M_PI, 2.)
                     + pow(sin(M_PI * y) / M_PI, 2.)),
                    2.)
                / Ngrid_ / Ngrid_ * mass_regulator * mass_regulator);
        } else {
            double c = 0.2;
            double length = param_.getL();
            double phys_x = x / Ngrid_ * length;  // in fm
            double phys_y = y / Ngrid_ * length;
            double phys_r2 = phys_x * phys_x + phys_y * phys_y;
            int Nf = 3;

            // Alphas in physical units! Lambda2 is lambda_QCD^2 in GeV
            double alphas =
                4. * M_PI
                / ((11.0 * param_.getNc() - 2.0 * Nf) / 3.
                   * log(pow(
                       (pow(mu0 * mu0 / Lambda2, 1. / c)
                        + pow(
                            4. / (phys_r2 * Lambda2 * fmgev * fmgev), 1. / c)),
                       c)));

            // discretization without singularities
            K_[pos]->push_back(
                sqrt(alphas)
                * (cos(M_PI * y) * (sin(2. * M_PI * x) / (2. * M_PI))
                   / ((pow(sin(M_PI * x) / M_PI, 2.)
                       + pow(sin(M_PI * y) / M_PI, 2.))))
                / Ngrid_ * mass_regulator);
            K_[pos]->push_back(
                sqrt(alphas)
                * (cos(M_PI * x) * (sin(2. * M_PI * y) / (2. * M_PI))
                   / ((pow(sin(M_PI * x) / M_PI, 2.)
                       + pow(sin(M_PI * y) / M_PI, 2.))))
                / Ngrid_ * mass_regulator);
            S_[pos]->push_back(
                alphas
                * (pow(cos(M_PI * y), 2.)
                       * pow(sin(2. * M_PI * x) / (2. * M_PI), 2.)
                   + pow(cos(M_PI * x), 4.)
                         * pow(sin(2. * M_PI * y) / (2. * M_PI), 2.))
                / pow(
                    (pow(sin(M_PI * x) / M_PI, 2.)
                     + pow(sin(M_PI * y) / M_PI, 2.)),
                    2.)
                / Ngrid_ / Ngrid_ * mass_regulator);
        }
    }
    fft_ptr_->fftnVector(K_, K_, nn_, 1);
    if (param_.getSimpleLangevin()) {
        fft_ptr_->fftnVector(S_, S_, nn_, 1);
    }
    return true;
}

double JIMWLK::getMassRegulator(const double x, const double y) const {
    // if m suppresses long distance tails,
    // K is multiplied by this, which is m*r*K_1(m*r)
    double mass_regulator = 1.0;
    double m = param_.getm();
    if (m < 1e-16) {
        return mass_regulator;
    }
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
