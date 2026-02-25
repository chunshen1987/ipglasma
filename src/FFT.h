// FFT.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef FFT_H
#define FFT_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fftw3.h>

#include <algorithm>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
    // assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(
        a.begin(), a.end(), b.begin(), std::back_inserter(result),
        std::plus<T>());
    return result;
}

template <typename T>
std::vector<T> operator*(
    const std::vector<T> &a, const std::complex<double> &b) {
    std::vector<T> result;
    result = a;
    int size = a.size();
    for (int i = 0; i < size; i++) {
        result.at(i) = b * a.at(i);
    }
    return result;
}

template <typename T>
std::vector<T> operator/(const std::vector<T> &a, const double b) {
    std::vector<T> result;
    result = a;
    int size = a.size();
    for (int i = 0; i < size; i++) {
        result.at(i) = a.at(i) / b;
    }
    return result;
}

class FFT {
  private:
    fftw_complex *input, *output;
    fftw_complex *inputMany, *outputMany;
    fftw_plan p, pback, pmany, pmanyback;

    // Batched plans for array FFTs — created lazily on first use
    int nn_stored_[2];
    int matDim_;  // Nc*Nc: number of matrix components
    bool arrayPlansInitialized_ = false;
    fftw_plan pArrayNc2m1_, pArrayNc2m1Back_;    // Nc^2-1 transforms
    fftw_plan pArray2Nc2m1_, pArray2Nc2m1Back_;  // 2*(Nc^2-1) transforms

    void initArrayPlans() {
        if (arrayPlansInitialized_) return;
        int Nc2m1 = matDim_ - 1;
        int maxBatch = std::max(matDim_, 2 * Nc2m1);
        int ntot = nn_stored_[0] * nn_stored_[1];

        // Expand buffers from matDim_ to maxBatch to fit the array plans.
        // Rebuild pmany/pmanyback since their buffer pointer becomes invalid.
        fftw_free(inputMany);
        fftw_free(outputMany);
        inputMany  = 
            (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ntot * maxBatch);
        outputMany = 
            (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ntot * maxBatch);

        fftw_destroy_plan(pmany);
        fftw_destroy_plan(pmanyback);
        pmany = fftw_plan_many_dft(
            2, nn_stored_, matDim_, inputMany, nn_stored_, 1, ntot,
            outputMany, nn_stored_, 1, ntot,
            FFTW_FORWARD, FFTW_MEASURE);
        pmanyback = fftw_plan_many_dft(
            2, nn_stored_, matDim_, inputMany, nn_stored_, 1, ntot,
            outputMany, nn_stored_, 1, ntot,
            FFTW_BACKWARD, FFTW_MEASURE);

        // Array plans for Nc^2-1 simultaneous FFTs
        pArrayNc2m1_ = fftw_plan_many_dft(
            2, nn_stored_, Nc2m1, inputMany, nn_stored_, 1, ntot,
            outputMany, nn_stored_, 1, ntot,
            FFTW_FORWARD, FFTW_MEASURE);
        pArrayNc2m1Back_ = fftw_plan_many_dft(
            2, nn_stored_, Nc2m1, inputMany, nn_stored_, 1, ntot,
            outputMany, nn_stored_, 1, ntot,
            FFTW_BACKWARD, FFTW_MEASURE);

        // Array plans for 2*(Nc^2-1) simultaneous FFTs
        pArray2Nc2m1_ = fftw_plan_many_dft(
            2, nn_stored_, 2 * Nc2m1, inputMany, nn_stored_, 1, ntot,
            outputMany, nn_stored_, 1, ntot,
            FFTW_FORWARD, FFTW_MEASURE);
        pArray2Nc2m1Back_ = fftw_plan_many_dft(
            2, nn_stored_, 2 * Nc2m1, inputMany, nn_stored_, 1, ntot,
            outputMany, nn_stored_, 1, ntot,
            FFTW_BACKWARD, FFTW_MEASURE);

        arrayPlansInitialized_ = true;
    }

  public:
    // Constructor. Nc is stored for lazy creation of batched plans.
    FFT(const int nn[], int Nc = 3) {
        nn_stored_[0] = nn[0];
        nn_stored_[1] = nn[1];
        matDim_ = Nc * Nc;

        input =
            (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nn[0] * nn[1]);
        output =
            (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nn[0] * nn[1]);
        p = fftw_plan_dft_2d(
            nn[0], nn[1], input, output, FFTW_FORWARD, FFTW_MEASURE);
        pback = fftw_plan_dft_2d(
            nn[0], nn[1], input, output, FFTW_BACKWARD, FFTW_MEASURE);

        // Allocate many-buffers for matrix FFTs (matDim_ = Nc*Nc = 9 for SU(3)).
        // initArrayPlans() will reallocate to maxBatch when first needed by JIMWLK.
        inputMany = (fftw_complex *)fftw_malloc(
            sizeof(fftw_complex) * nn[0] * nn[1] * matDim_);
        outputMany = (fftw_complex *)fftw_malloc(
            sizeof(fftw_complex) * nn[0] * nn[1] * matDim_);

        // Matrix many plans
        pmany = fftw_plan_many_dft(
            2, nn, matDim_, inputMany, nn, 1, nn[0] * nn[1],
            outputMany, nn, 1, nn[0] * nn[1],
            FFTW_FORWARD, FFTW_MEASURE);
        pmanyback = fftw_plan_many_dft(
            2, nn, matDim_, inputMany, nn, 1, nn[0] * nn[1],
            outputMany, nn, 1, nn[0] * nn[1],
            FFTW_BACKWARD, FFTW_MEASURE);
    };
    // Destructor
    ~FFT() {
        if (arrayPlansInitialized_) {
            fftw_destroy_plan(pArray2Nc2m1Back_);
            fftw_destroy_plan(pArray2Nc2m1_);
            fftw_destroy_plan(pArrayNc2m1Back_);
            fftw_destroy_plan(pArrayNc2m1_);
        }
        fftw_destroy_plan(pmany);
        fftw_destroy_plan(pmanyback);
        fftw_destroy_plan(p);
        fftw_destroy_plan(pback);
        fftw_free(input);
        fftw_free(output);
        fftw_free(inputMany);
        fftw_free(outputMany);
    };
    void fftnVector(
        vector<complex<double>> **data, vector<complex<double>> **outdata,
        const int nn[], const int isign);
    void fftnArray(
        complex<double> **data, complex<double> **outdata, const int nn[],
        const int isign, const int mDim);
    void fftnArrayMany(
        complex<double> **data, complex<double> **outdata, const int nn[],
        const int isign, const int mDim);

    template <class T>
    void fftn(T **data, T **outdata, const int nn[], const int isign);
    template <class T>
    void fftnMany(T **data, T **outdata, const int nn[], const int isign);

    void fftnComplex(
        complex<double> *data, complex<double> *outdata, const int nn[],
        const int isign);
};

#endif  // FFT_H
