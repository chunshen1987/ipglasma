#ifndef Cell_h
#define Cell_h

#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "Matrix.h"

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

    std::vector<double> Tmunu;

    std::vector<double> pimunu;

    std::vector<double> umu;

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

    void setTtautau(const double in) { Tmunu[0] = in; };
    double getTtautau() { return Tmunu[0]; };
    void setTxx(const double in) { Tmunu[4] = in; };
    double getTxx() { return Tmunu[4]; };
    void setTyy(const double in) { Tmunu[7] = in; };
    double getTyy() { return Tmunu[7]; };
    void setTxy(const double in) { Tmunu[5] = in; };
    double getTxy() { return Tmunu[5]; };
    void setTetaeta(const double in) { Tmunu[9] = in; };
    double getTetaeta() { return Tmunu[9]; };
    void setTtaux(const double in) { Tmunu[1] = in; };
    double getTtaux() { return Tmunu[1]; };
    void setTtauy(const double in) { Tmunu[2] = in; };
    double getTtauy() { return Tmunu[2]; };
    void setTtaueta(const double in) { Tmunu[3] = in; };
    double getTtaueta() { return Tmunu[3]; };
    void setTxeta(const double in) { Tmunu[6] = in; };
    double getTxeta() { return Tmunu[6]; };
    void setTyeta(const double in) { Tmunu[8] = in; };
    double getTyeta() { return Tmunu[8]; };

    void setpitautau(const double in) { pimunu[0] = in; };
    double getpitautau() { return pimunu[0]; };
    void setpixx(const double in) { pimunu[4] = in; };
    double getpixx() { return pimunu[4]; };
    void setpiyy(const double in) { pimunu[7] = in; };
    double getpiyy() { return pimunu[7]; };
    void setpixy(const double in) { pimunu[5] = in; };
    double getpixy() { return pimunu[5]; };
    void setpietaeta(const double in) { pimunu[9] = in; };
    double getpietaeta() { return pimunu[9]; };
    void setpitaux(const double in) { pimunu[1] = in; };
    double getpitaux() { return pimunu[1]; };
    void setpitauy(const double in) { pimunu[2] = in; };
    double getpitauy() { return pimunu[2]; };
    void setpitaueta(const double in) { pimunu[3] = in; };
    double getpitaueta() { return pimunu[3]; };
    void setpixeta(const double in) { pimunu[6] = in; };
    double getpixeta() { return pimunu[6]; };
    void setpiyeta(const double in) { pimunu[8] = in; };
    double getpiyeta() { return pimunu[8]; };

    void setutau(const double in) { umu[0] = in; };
    double getutau() { return umu[0]; };
    void setux(const double in) { umu[1] = in; };
    double getux() { return umu[1]; };
    void setuy(const double in) { umu[2] = in; };
    double getuy() { return umu[2]; };
    void setueta(const double in) { umu[3] = in; };
    double getueta() { return umu[3]; };

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
