#ifndef glauber_h  // avoids multiple inclusions of the header file
#define glauber_h

#include <string>

#include "Random.h"

#define TOL (1.0e-6)
#define tiny (1.0e-10)
#define limit 10000

struct ReturnValue {
    double x;
    double y;
    double z;
    double phi;
    int collided;
    bool proton;
    // int acceptances;
};

typedef struct nucleus {
    std::string name;
    double A;
    double Z;
    int AnumFunc;
    int AnumFuncIntegrand;
    int DensityFunc;
    double w_WS;
    double a_WS;
    double R_WS;
    double rho_WS;
    double beta2;
    double beta3;
    double beta4;
    double gamma;
    bool forceDminFlag;
    double d_min;
    double dR_np;
    double da_np;
} Nucleus;

typedef struct data {
    double SigmaNN;
    Nucleus Target;
    Nucleus Projectile;
    double SCutOff;
    int InterMax;
    /* trap door */
} Data;

class Glauber {
  private:
  public:
    typedef double (*ptr_func)(double);

    double AnumR, NuInS_S;
    Nucleus *Nuc_WS;
    Data GlauberData;
    ptr_func tempFunc;
    double b;  // impact parameter
    double currentTAB;
    double currentA1;
    double currentA2;
    double currentZ1;
    double currentZ2;

    Glauber() {};
    ~Glauber() { remove("tmp.dat"); }

    int nucleusA1() const { return static_cast<int>(currentA1); }
    int nucleusA2() const { return static_cast<int>(currentA2); }
    int nucleusZ1() const { return static_cast<int>(currentZ1); }
    int nucleusZ2() const { return static_cast<int>(currentZ2); }
    int IsFile(char *file_name);
    void FindNucleusData(
        Nucleus *nucleus, std::string target, std::string file_name, int rank);
    void FindNucleusData2(
        Nucleus *nucleus, std::string name, bool setWSDeformParams, double R_WS,
        double a_WS, double beta2, double beta3, double beta4, double gamma,
        bool force_dmin, double d_min, double dR_np, double da_np);
    void PrintGlauberData();
    void PrintNucleusData(Nucleus *nucleus);
    int LinearFindXorg(double x, double *Vx, int ymax);
    double FourPtInterpolate(
        double x, double *Vx, double *Vy, double h, int x_org);
    void MakeCoeff(
        double *a, double *b, double *c, double *d, double *Vy, double h,
        int x_org);
    double VInterpolate(double x, double *Vx, double *Vy, int ymax);
    int FindXorg(double x, double *Vx, int ymax);
    double *MakeVx(double down, double up, int maxi_num);
    double *MakeVy(double *vx, int maxi_num);
    double *ReadInVx(char *, int maxi_num, int quiet);
    double *ReadInVy(char *, int maxi_num, int quiet);

    double InterNuPInSP(double s);
    double InterNuTInST(double s);
    void CalcRho(Nucleus *nucleus);
    double NuInS(double s);

    double Anum3Fermi(double R_WS);
    double Anum3FermiInt(double xi);
    double NuInt3Fermi(double xi);
    double Anum3Gauss(double R_WS);
    double Anum3GaussInt(double xi);
    double NuInt3Gauss(double xi);
    double Anum2HO();
    double Anum2HOInt(double xi);
    double NuInt2HO(double xi);
    double AnumHulthen();
    double AnumHulthenInt();
    double NuIntHulthen(double xi);

    double integral(int id, double down, double up, double tol, int *count);
    double qnc7(
        int id, double tol, double down, double dx, double *f_of,
        double pre_sum, double area, int *count);
    double OLSIntegrand(double s);
    double TAB();
    double PAB(double x, double y);
    void initGlauber(
        double SigmaNN, std::string Target, std::string Projectile, double inb,
        bool setWSDeformParams, double R_WS, double a_WS, double beta2,
        double beta3, double beta4, double gamma, bool force_dmin, double d_min,
        double dR_np, double da_np, int imax);
    double areaTA(double x, double A);
    ReturnValue SampleTARejection(Random *random, int PorT);
};
#endif
