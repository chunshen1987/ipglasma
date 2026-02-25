#ifndef Cell_h
#define Cell_h

#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <memory>

#include "Matrix.h"

// Hydro tensor data — bundled into one allocation, only present in mode=1 cells
struct HydroData {
    std::array<double, 10> Tmunu = {};
    std::array<double, 10> pimunu = {};
    std::array<double, 4> umu = {};
};

class Cell {
  private:
    int mode_;
    double epsilon;  // energy density after collision

    // nucleus A
    double g2mu2A;  // color charge density of nucleus A
    double TpA;     // sum over the proton T(b) in this cell for nucleus A
    Matrix *U;  // U is in the fundamental rep. (Nc*Nc matrix) // duobles as x
                // component of electric field

    // nucleus B
    double g2mu2B;  // color charge density of nucleus B
    double TpB;     // sum over the proton T(b) in this cell for nucleus B
    Matrix *U2;  // Ui is the initial U in the fundamental rep. (Nc*Nc matrix)
                 // // doubles as y component of electric field

    Matrix *Ux;  // U is in the fundamental rep. (Nc*Nc matrix)
    Matrix *Uy;  // U is in the fundamental rep. (Nc*Nc matrix)

    Matrix *Ux1;  // U is in the fundamental rep. (Nc*Nc matrix) nucleus 1 (also
                  // room to save g, the gauge fixing matrix)
    Matrix *Uy1;  // U is in the fundamental rep. (Nc*Nc matrix) nucleus 1 (also
                  // room to save Uplaq, the plaquette)

    Matrix *Ux2;  // U is in the fundamental rep. (Nc*Nc matrix) nucleus 2
                  // (doubles as longitudinal electric field pi)
    Matrix *Uy2;  // U is in the fundamental rep. (Nc*Nc matrix) nucleus 2
                  // (doubles as scalar field (longitudinal) )

    //  bool parity; // Parity of the cell (needed for Gauge fixing)

    std::unique_ptr<HydroData> hydro_;  // only allocated for mode=1 cells

  public:
    Cell(const int Nc, const int mode);
    ~Cell();

    //  void setParity(bool in) { parity = in; };
    //  bool getParity() { return parity; };

    void setg2mu2A(double in) { g2mu2A = in; };
    void setg2mu2B(double in) { g2mu2B = in; };

    double getg2mu2A() { return g2mu2A; };
    double getg2mu2B() { return g2mu2B; };

    void setTpA(double in) { TpA = in; };
    void setTpB(double in) { TpB = in; };

    double getTpA() { return TpA; };
    double getTpB() { return TpB; };

    void setU(const Matrix &x) { *U = x; };
    void setU2(const Matrix &x) { *U2 = x; };
    void setUplaq(const Matrix &x) {
        *Uy1 = x;
    };  // using unused Uy1 to store Uplaq

    void setUx(const Matrix &x) { *Ux = x; };
    void setUy(const Matrix &x) { *Uy = x; };
    void setUx1(const Matrix &x) { *Ux1 = x; };
    void setUy1(const Matrix &x) { *Uy1 = x; };
    void setUx2(const Matrix &x) { *Ux2 = x; };
    void setUy2(const Matrix &x) { *Uy2 = x; };

    void setEpsilon(const double in) { epsilon = in; };
    double getEpsilon() { return epsilon; };

    void setTtautau(const double in) { hydro_->Tmunu[0] = in; };
    double getTtautau() { return hydro_->Tmunu[0]; };
    void setTxx(const double in) { hydro_->Tmunu[4] = in; };
    double getTxx() { return hydro_->Tmunu[4]; };
    void setTyy(const double in) { hydro_->Tmunu[7] = in; };
    double getTyy() { return hydro_->Tmunu[7]; };
    void setTxy(const double in) { hydro_->Tmunu[5] = in; };
    double getTxy() { return hydro_->Tmunu[5]; };
    void setTetaeta(const double in) { hydro_->Tmunu[9] = in; };
    double getTetaeta() { return hydro_->Tmunu[9]; };
    void setTtaux(const double in) { hydro_->Tmunu[1] = in; };
    double getTtaux() { return hydro_->Tmunu[1]; };
    void setTtauy(const double in) { hydro_->Tmunu[2] = in; };
    double getTtauy() { return hydro_->Tmunu[2]; };
    void setTtaueta(const double in) { hydro_->Tmunu[3] = in; };
    double getTtaueta() { return hydro_->Tmunu[3]; };
    void setTxeta(const double in) { hydro_->Tmunu[6] = in; };
    double getTxeta() { return hydro_->Tmunu[6]; };
    void setTyeta(const double in) { hydro_->Tmunu[8] = in; };
    double getTyeta() { return hydro_->Tmunu[8]; };

    void setpitautau(const double in) { hydro_->pimunu[0] = in; };
    double getpitautau() { return hydro_->pimunu[0]; };
    void setpixx(const double in) { hydro_->pimunu[4] = in; };
    double getpixx() { return hydro_->pimunu[4]; };
    void setpiyy(const double in) { hydro_->pimunu[7] = in; };
    double getpiyy() { return hydro_->pimunu[7]; };
    void setpixy(const double in) { hydro_->pimunu[5] = in; };
    double getpixy() { return hydro_->pimunu[5]; };
    void setpietaeta(const double in) { hydro_->pimunu[9] = in; };
    double getpietaeta() { return hydro_->pimunu[9]; };
    void setpitaux(const double in) { hydro_->pimunu[1] = in; };
    double getpitaux() { return hydro_->pimunu[1]; };
    void setpitauy(const double in) { hydro_->pimunu[2] = in; };
    double getpitauy() { return hydro_->pimunu[2]; };
    void setpitaueta(const double in) { hydro_->pimunu[3] = in; };
    double getpitaueta() { return hydro_->pimunu[3]; };
    void setpixeta(const double in) { hydro_->pimunu[6] = in; };
    double getpixeta() { return hydro_->pimunu[6]; };
    void setpiyeta(const double in) { hydro_->pimunu[8] = in; };
    double getpiyeta() { return hydro_->pimunu[8]; };

    void setutau(const double in) { hydro_->umu[0] = in; };
    double getutau() { return hydro_->umu[0]; };
    void setux(const double in) { hydro_->umu[1] = in; };
    double getux() { return hydro_->umu[1]; };
    void setuy(const double in) { hydro_->umu[2] = in; };
    double getuy() { return hydro_->umu[2]; };
    void setueta(const double in) { hydro_->umu[3] = in; };
    double getueta() { return hydro_->umu[3]; };

    Matrix &getg() const { return *Ux1; };  // use unused Ux1 to store g
    Matrix &getU() const { return *U; };
    Matrix &getUx() const { return *Ux; };
    Matrix &getUy() const { return *Uy; };
    Matrix &getU2() const { return *U2; };
    Matrix &getUx1() const { return *Ux1; };
    Matrix &getUy1() const { return *Uy1; };
    Matrix &getUx2() const { return *Ux2; };
    Matrix &getUy2() const { return *Uy2; };
    Matrix &getUplaq() const { return *Uy1; };  // use unused Uy1 to store Uplaq

    void setE1(const Matrix &x) { *U = x; };  // use unused U to store E1
    Matrix &getE1() const { return *U; };
    void setE2(const Matrix &x) { *U2 = x; };  // use unused U2 to store E2
    Matrix &getE2() const { return *U2; };
    void setphi(const Matrix &x) { *Uy2 = x; };  // use unused Uy2 to store phi
    Matrix &getphi() const { return *Uy2; };
    void setpi(const Matrix &x) { *Ux2 = x; };  // use unused Ux2 to store pi
    Matrix &getpi() const { return *Ux2; };
    void setg(const Matrix &x) { *Ux1 = x; };  // using unused Ux1 to store g

    //  void computeAdjointU();
};

class SmallCell {
  private:
    Matrix *buffer1;
    Matrix *buffer2;

  public:
    SmallCell(const int Nc);
    ~SmallCell();

    Matrix &getbuffer1() const { return *buffer1; };
    void setbuffer1(const Matrix &x) { *buffer1 = x; };
    Matrix &getbuffer2() const { return *buffer2; };
    void setbuffer2(const Matrix &x) { *buffer2 = x; };
};

#endif
