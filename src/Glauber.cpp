#include "Glauber.h"
#include "Util.h"

using namespace std;

Glauber::Glauber()
{
  setup = new Setup;
}

// destructor
Glauber::~Glauber() {
    delete setup;
    remove("tmp.dat");
}

void Glauber::FindNucleusData2(Nucleus *nucleus, string name, int rank)
{
  string densityFunction;
  if(name.compare("Au") == 0)
    {
      nucleus->A = 197;
      nucleus->Z =  79;
      densityFunction = "3Fermi";
      nucleus->R_WS = 6.37;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.535;
      nucleus->beta2 = -0.13;
      nucleus->beta4 = -0.03;
    }
  else if(name.compare("Pb") == 0)
    {
      nucleus->A = 208.;
      nucleus->Z =  82.;
      densityFunction = "3Fermi";
      nucleus->R_WS = 6.62;
      nucleus->w_WS = 0.;
      nucleus->a_WS = 0.546;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("p") == 0)
    {
      nucleus->A = 1.;
      nucleus->Z =  1.;
      densityFunction = "3Fermi";
      nucleus->R_WS = 1.;
      nucleus->w_WS = 0.;
      nucleus->a_WS = 1.;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("He3") == 0)
    {
      nucleus->A = 3;
      nucleus->Z =  2;
      densityFunction = "readFromFile";
      nucleus->R_WS = 0;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("d") == 0)
    {
      nucleus->A = 2;
      nucleus->Z =  1;
      densityFunction = "Hulthen";
      nucleus->R_WS = 1.0;
      nucleus->w_WS = 1.18;
      nucleus->a_WS = 0.228;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("C") == 0)
    {
      nucleus->A = 12;
      nucleus->Z =  6;
      densityFunction = "2HO";
      nucleus->R_WS = 2.44;
      nucleus->w_WS = 1.403;
      nucleus->a_WS = 1.635;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("O") == 0)
    {
      nucleus->A = 16;
      nucleus->Z =  8;
      densityFunction = "3Fermi";
      nucleus->R_WS = 2.608;
      nucleus->w_WS = -0.051;
      nucleus->a_WS = 0.513;
      nucleus->beta2 = -0.01;                 // from arXiv:1508.06294
      nucleus->beta4 = -0.122;                // from arXiv:1508.06294
    }
  else if(name.compare("S") == 0)
    {
      nucleus->A = 32;
      nucleus->Z =  16;
      densityFunction = "3Gauss";
      nucleus->R_WS = 2.54;
      nucleus->w_WS = 0.16;
      nucleus->a_WS = 2.191;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("W") == 0)
    {
      nucleus->A = 184;
      nucleus->Z =  74;
      densityFunction = "3Fermi";
      nucleus->R_WS = 6.51;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.535;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("Al") == 0)
    {
      nucleus->A = 27;
      nucleus->Z =  13;
      densityFunction = "3Fermi";
      nucleus->R_WS = 3.07;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.519;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("Ca") == 0)
    {
      nucleus->A = 40;
      nucleus->Z =  20;
      densityFunction = "3Fermi";
      nucleus->R_WS = 3.766;
      nucleus->w_WS = -0.161;
      nucleus->a_WS = 0.586;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("Cu") == 0)
    {
      nucleus->A = 63;
      nucleus->Z =  29;
      densityFunction = "3Fermi";
      nucleus->R_WS = 4.163;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.606;
      nucleus->beta2 = 0.162;
      nucleus->beta4 = 0.006;
    }
  else if(name.compare("Fe") == 0)
    {
      nucleus->A = 56;
      nucleus->Z =  26;
      densityFunction = "3Fermi";
      nucleus->R_WS = 4.106;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.519;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("Pt") == 0)
    {
      nucleus->A = 195;
      nucleus->Z =  78;
      densityFunction = "3Fermi";
      nucleus->R_WS = 6.78;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.54;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("U") == 0)
    {
      nucleus->A = 238;
      nucleus->Z =  92;
      densityFunction = "3Fermi";
      nucleus->R_WS = 6.81;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.55;
      nucleus->beta2 = 0.28;
      nucleus->beta4 = 0.093;
    }
  else if(name.compare("Ru") == 0)
    {
      nucleus->A = 96;
      nucleus->Z = 44;
      densityFunction = "3Fermi";
      nucleus->R_WS = 5.085;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.46;
      nucleus->beta2 = 0.158;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("Zr") == 0)
    {
      nucleus->A = 96;
      nucleus->Z = 40;
      densityFunction = "3Fermi";
      nucleus->R_WS = 5.02;
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.46;
      nucleus->beta2 = 0.0;
      nucleus->beta4 = 0.0;
    }
  else if(name.compare("Xe") == 0)
    {
      nucleus->A = 129;
      nucleus->Z = 54;
      densityFunction = "3Fermi";
      nucleus->R_WS = 5.42; // 5.36       // new values from arXiv:1508.06294
      nucleus->w_WS = 0;
      nucleus->a_WS = 0.57; // 0.590;     // new values from arXiv:1508.06294
      nucleus->beta2 = 0.162;             // from arXiv:1508.06294
      nucleus->beta4 = -0.003;            // from arXiv:1508.06294
    }

  nucleus->rho_WS = nucleus->R_WS;
  
  if(densityFunction.compare("2HO")==0) 
    {
      nucleus->AnumFunc = 1; //Anum2HO;
      nucleus->AnumFuncIntegrand = 1; //Anum2HOInt;
      nucleus->DensityFunc = 1; //NuInt2HO;
    }
  else if(densityFunction.compare("3Gauss")==0) 
    {
      nucleus->AnumFunc = 2; //Anum3Gauss;
      nucleus->AnumFuncIntegrand = 2; //Anum3GaussInt;
      nucleus->DensityFunc = 2; //NuInt3Gauss;
    }
  else if(densityFunction.compare("3Fermi")==0) 
    {
      nucleus->AnumFunc = 3; //Anum3Fermi;
      nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
      nucleus->DensityFunc = 3; //NuInt3Fermi;
    }
  else if(densityFunction.compare("Hulthen")==0) 
    {
      nucleus->AnumFunc = 8;//AnumHulthen;
      nucleus->AnumFuncIntegrand = 8; //AnumHulthenInt;
      nucleus->DensityFunc = 8; //NuIntHulthen;  
    }
  else if(densityFunction.compare("readFromFile")==0) 
    {
      nucleus->AnumFunc = 1;
      nucleus->AnumFuncIntegrand = 1;
      nucleus->DensityFunc = 1;  
    }
}/* FindNucleusData2 */



void Glauber::PrintGlauberData()
{
 fprintf(stderr, "GlauberData.SigmaNN = %e\n",  GlauberData.SigmaNN);
 fprintf(stderr, "GlauberData.InterMax = %d\n", GlauberData.InterMax);
 fprintf(stderr, "GlauberData.SCutOff = %f\n", GlauberData.SCutOff);
}/* PrintGlauberData */


void Glauber::PrintNucleusData(Nucleus *nucleus)
{
  cout << "Nucleus Name: " << nucleus->name << endl;
  cout << " Nucleus.A = " << nucleus->A << endl;
  cout << " Nucleus.Z = " << nucleus->Z << endl;
  cout << " Nucleus.w_WS = " << nucleus->w_WS << endl;
  cout << " Nucleus.a_WS = " << nucleus->a_WS << endl;
  cout << " Nucleus.R_WS = " << nucleus->R_WS << endl;

}


int Glauber::LinearFindXorg(double x, double *Vx, int ymax)
{
/* finds the first of the 4 points, x is between the second and the third */
 
 int x_org;
 double nx;

 nx = ymax*(x - Vx[0])/(Vx[ymax] - Vx[0]);
 
 x_org = (int) nx;
 x_org -= 1;

 if( x_org <= 0 ) return 0;
 else if(x_org >= ymax - 3) return ymax - 3;
 else return x_org;

}/* Linear Find Xorg */


double Glauber::FourPtInterpolate(double x, double *Vx, double *Vy, double h, int x_org, int ymax)
{
 /* interpolating points are x_org, x_org+1, x_org+2, x_org+3 */
 /* cubic polynomial approximation */

 double a, b, c, d, f;

 MakeCoeff(&a, &b, &c, &d,  Vy, Vx, h, x_org);
 
 f = a*pow(x - Vx[x_org], 3.);
 f += b*pow(x - Vx[x_org], 2.);
 f += c*(x - Vx[x_org]);
 f += d;
 
 return f;
}/* FourPtInterpolate */


void Glauber::MakeCoeff(double *a, double *b, double *c, double *d, 
			double *Vy, double *Vx, double h, int x_org)
{
 double f0, f1, f2, f3;

 f0 = Vy[x_org];
 f1 = Vy[x_org+1];
 f2 = Vy[x_org+2];
 f3 = Vy[x_org+3];

 *a =  (-f0 + 3.0*f1 - 3.0*f2 + f3)/(6.0*h*h*h);
 
 *b =  (2.0*f0 - 5.0*f1 + 4.0*f2 - f3)/(2.0*h*h);

 *c =  (-11.0*f0 + 18.0*f1 - 9.0*f2 + 2.0*f3)/(6.0*h);

 *d = f0;

}/* MakeCoeff */

int Glauber::FindXorg(double x, double *Vx, int ymax)
{
 int i, x_org;

 i = 0;
 while(Vx[i] < x) i++;
 
 x_org = i - 2;
 
 if( x_org <= 1 ) return 1;
 else if(x_org >= ymax - 3) return ymax - 3;
 else return x_org;

}/* Find Xorg */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::VInterpolate(double x, double *Vx, double *Vy, int ymax)
{
 int x_org;
 double h;

 if( (x < Vx[0])||(x > Vx[ymax]) )
  {
   fprintf(stderr, 
           "VInterpolate: x = %le is outside the range (%le, %le).\n", 
	    x,Vx[0],Vx[ymax]);
   fprintf(stderr, "This can't happen.  Exiting...\n");
   exit(0);
  }

/* we only deal with evenly spaced Vx */
/* x_org is the first of the 4 points */

 x_org = LinearFindXorg(x, Vx, ymax);

 h = (Vx[ymax] - Vx[0])/ymax;

 return FourPtInterpolate(x, Vx, Vy, h, x_org, ymax);

}/* VInterpolate */


double *Glauber::MakeVx(double down, double up, int maxi_num)
{
 static double dx, *vx;
 int i;

 vx = Util::vector_malloc(maxi_num + 1);
 dx = (up - down)/maxi_num;
 for(i=0; i<=maxi_num; i++)
  {
   vx[i] = dx*i;
  }
 return vx;
 

}/* MakeVx */


double *Glauber::MakeVy(double *vx, int maxi_num)
{
 int i;
 static double *vy;

 // if(maxi_num > 200) di = 100;
 // if(maxi_num <= 200) di = 20;
 
 vy = Util::vector_malloc(maxi_num + 1);

 // ofstream data_file(st.c_str());
 
 // data_file << "EndOfData" << endl;
 
 for(i=0; i<=maxi_num; i++)
  {
   vy[i] = NuInS(vx[i]);
 //   if(i % di == 0)
//     {
//       cerr << st << "[" << i << "] = " << vy[i] << endl; 
//     }
   //data_file << vx[i] << " " << vy[i] << endl;
  }

 //data_file.close();

 return vy;
}/* MakeVy */


double *Glauber::ReadInVx(char *file_name, int maxi_num, int quiet)
{
 static double x, *vx;
 int i;
 FILE *input;
 static char *s, *sx;
 //int bytes_read;
 s = Util::char_malloc(120);
 sx = Util::char_malloc(120);
 
 vx = Util::vector_malloc(maxi_num + 1);

 if(quiet == 1)
  {
   fprintf(stderr, "Reading in Vx from %s ...\n", file_name);
   }
 
 input = fopen(file_name, "r");
 fscanf(input, "%s", s);
 while(strcmp(s, "EndOfData") != 0)
  {
   fscanf(input, "%s", sx);
   fscanf(input, "%s", s);
  }

 for(i=0; i<=maxi_num; i++)
  {
   fscanf(input, "%lf", &x);
   vx[i] = x;
   fscanf(input, "%lf", &x);
  }
 fclose(input);

 Util::char_free(sx);
 Util::char_free(s);
 return vx;

}/* ReadInVx */


double *Glauber::ReadInVy(char *file_name, int maxi_num, int quiet) {
    static double y, *vy;
    int i;
    FILE *input;
    static char *s, *sy;
    // int bytes_read;
    s = Util::char_malloc(120);
    sy = Util::char_malloc(120);

    vy = Util::vector_malloc(maxi_num + 1);

    if (quiet == 1) {
        fprintf(stderr, "Reading in Vy from %s ...\n", file_name);
    }

    input = fopen(file_name, "r");
    fscanf(input, "%s", s);
    while (strcmp(s, "EndOfData") != 0) {
        fscanf(input, "%s", sy);
        fscanf(input, "%s", s);
    }

    for (i=0; i<=maxi_num; i++) {
        fscanf(input, "%lf", &y);
        fscanf(input, "%lf", &y);
        vy[i] = y;
    }
    fclose(input);

    Util::char_free(s);
    Util::char_free(sy);
    return vy;
}  /* ReadInVy */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::InterNuPInSP(double s)
{
 double y;
 static int ind = 0;
 static double up, down; 
 static int maxi_num; 
 static double *vx, *vy;
 ind++;

 if(GlauberData.Projectile.A == 1.0) return 0.0;

 if(ind == 1)
  {
    CalcRho(&(GlauberData.Projectile));
    up = 2.0*GlauberData.SCutOff;
    down = 0.0; 
    maxi_num = GlauberData.InterMax; 
    vx = MakeVx(down, up, maxi_num);
    vy = MakeVy(vx, maxi_num);
  }/* if ind */
 
 if(s > up) return 0.0;
 else{
   y = VInterpolate(s, vx, vy, maxi_num);
   if( y < 0.0 ) return 0.0; 
   else 
     {
       return y;
     }
  }
}/* InterNuPInSP */

double Glauber::InterNuTInST(double s)
{
  double y;
  static int ind = 0;
  static double up, down; 
  static int maxi_num; 
  static double *vx, *vy;
  
  ind++;
  if(GlauberData.Target.A == 1.0) return 0.0;
  
  if(ind == 1) 
    {
      CalcRho(&(GlauberData.Target));
      
      up = 2.0*GlauberData.SCutOff;
      down = 0.0; 
      maxi_num = GlauberData.InterMax; 
      
      vx = MakeVx(down, up, maxi_num);
      vy = MakeVy(vx, maxi_num);
    }/* if ind */
  
  // cout << *vx << " " << *vy << endl;
  
  if(s > up) return 0.0;
  else{
    y = VInterpolate(s, vx, vy, maxi_num);
   if( y < 0.0 ) return 0.0; 
   else return y;
  }
}/* InterNuTInST */

void Glauber::CalcRho(Nucleus *nucleus)
{
 double f, R_WS;
/* to pass to AnumIntegrand */ 
 
 Nuc_WS = nucleus;
 

 R_WS = nucleus->R_WS;   
 
 if (nucleus->AnumFunc==1)
   f = Anum2HO(R_WS)/(nucleus->rho_WS);
 else if (nucleus->AnumFunc==2)
   f = Anum3Gauss(R_WS)/(nucleus->rho_WS);
 else if (nucleus->AnumFunc==3)
   f = Anum3Fermi(R_WS)/(nucleus->rho_WS);
 else if (nucleus->AnumFunc==8)
   f = AnumHulthen(R_WS)/(nucleus->rho_WS);
 else
   f = Anum3Fermi(R_WS)/(nucleus->rho_WS);
 
 nucleus->rho_WS = (nucleus->A)/f;
 //cout << " nucleus->rho_WS=" << nucleus->rho_WS << endl;
}/* CalcRho */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::NuInS(double s)
{
 double y;
 int count;
 int id;

/* to pass to the DensityFunc's */
 NuInS_S = s;

 id = Nuc_WS->DensityFunc;

 count = 0;
 // cout << "calling integral" << endl;
 y = integral(id, 0.0, 1.0, TOL, &count); 
 
 return y;
}

double Glauber::Anum3Fermi(double R_WS)
{
 int count=0;
 double up, down, a_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 // w_WS = Nuc_WS->w_WS;
 rho = Nuc_WS->rho_WS;

/* to pass to Anumintegrand */
 AnumR = R_WS/a_WS;
 
 down = 0.0;
 up = 1.0;
 
 f = integral(4, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum3Fermi */


double Glauber::Anum3FermiInt(double xi) 
{
 double f;
 double r;
 double R_WS, w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0-tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 R_WS = AnumR;
 r = -log(xi);

 f = r*r;
 f *= ( 1.0 + w_WS*pow(r/R_WS, 2.) );
 f /= (xi + exp(-R_WS));
  
 return f;
}/* Anum3FermiInt */


double Glauber::NuInt3Fermi(double xi)
{
 double f;
 double c;
 double z, r, s;
 double w_WS, R_WS, a_WS, rho;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 
   
/* devide by a_WS, make life simpler */
   
   s = NuInS_S/a_WS;
   z = -log(xi);
   r = sqrt(s*s + z*z);
   R_WS /= a_WS;

   c = exp(-R_WS);
   
   f = 2.0*a_WS*rho*(GlauberData.SigmaNN);
   f *= 1.0 + w_WS*pow(r/R_WS, 2.);
   f /= xi + c*exp(s*s/(r + z)); 
 
 return f;
}/* NuInt3Fermi */


/* %%%%%%% 3 Parameter Gauss %%%%%%%%%%%% */


double Glauber::Anum3Gauss(double R_WS)
{
 int count=0;
 double up, down, a_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 // w_WS = Nuc_WS->w_WS;
 rho = Nuc_WS->rho_WS;

/* to pass to Anumintegrand */
 AnumR = R_WS/a_WS;
 
 down = 0.0;
 up = 1.0;
 
 f = integral(5, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum3Gauss */


double Glauber::Anum3GaussInt(double xi) 
{
 double y;
 double r_sqr;
 double R_WS, w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0 - tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 
 R_WS = AnumR;
 
 r_sqr = -log(xi);

 y = sqrt(r_sqr);
 
 y *= 1.0 + w_WS*r_sqr/pow(R_WS, 2.);

/* 2 comes from dr^2 = 2 rdr */

 y /= 2.0*(xi + exp(-R_WS*R_WS));
 
 return y;
}/* Anum3GaussInt */


double Glauber::NuInt3Gauss(double xi)
{
 double f;
 double c;
 double z_sqr, r_sqr, s;
 double w_WS, R_WS, a_WS, rho;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z*z/a/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 

/* devide by a_WS, make life simpler */

   s = NuInS_S/a_WS;
   z_sqr = -log(xi);
   r_sqr = s*s + z_sqr;
   R_WS /= a_WS;

   c = exp(-R_WS*R_WS);

   f = a_WS*rho*(GlauberData.SigmaNN);
   f *= 1.0 + w_WS*r_sqr/pow(R_WS, 2.);
   f /= sqrt(z_sqr)*(xi + c*exp(s*s)); 
 
 return f;
}/* NuInt3Gauss */



/* %%%%%%% 2 Parameter HO %%%%%%%%%%%% */


double Glauber::Anum2HO(double R_WS)
{
 int count=0;
 double up, down, a_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 // w_WS = Nuc_WS->w_WS; /* take this to be alpha */
 rho = Nuc_WS->rho_WS;

 down = 0.0;
 up = 1.0;
 
 f = integral(6, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum2HO */


double Glauber::Anum2HOInt(double xi) 
{
 double y;
 double r_sqr, r;
 double w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0 - tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 
 r_sqr = -log(xi);
 
 r = sqrt(r_sqr);
 
/* 2 comes from dr^2 = 2 rdr */
 
 y = r + w_WS*r*r_sqr;
 y /= 2.0;
 
 return y;
}/* Anum2HOInt */


double Glauber::NuInt2HO(double xi)
{
 double f;
 double z_sqr, r_sqr, s;
 double w_WS, a_WS, rho;
 
 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 // R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z*z/a/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 

/* devide by a_WS, make life simpler */
   
   s = NuInS_S/a_WS;
   z_sqr = -log(xi);
   r_sqr = s*s + z_sqr;

/* no need to divide by 2 here because -infty < z < infty and 
   we integrate only over positive z */
   
   if(z_sqr < 0.0) z_sqr = tiny;
   f = a_WS*rho*(GlauberData.SigmaNN);
   f *= (1.0 + w_WS*r_sqr)*exp(-s*s)/sqrt(z_sqr);
 
 return f;
}/* NuInt2HO */


/* %%%%%%% Hulthen %%%%%%%%%%%% */


double Glauber::AnumHulthen(double R_WS)
{
 double a_WS, b_WS, rho, f;

 /* R_WS is dummy */

 a_WS = Nuc_WS->a_WS;
 b_WS = Nuc_WS->w_WS; /* take this to be b */
 rho = Nuc_WS->rho_WS;

 /* density function has the form */
 /* rho = (1/r^2) (exp(-a r) - exp(-b r))^2 */
 /* int d^3r rho is given by */
 
 f =  2.0*(a_WS - b_WS)*(a_WS - b_WS)*M_PI/a_WS/b_WS/(a_WS + b_WS);
 
 f *= rho;

 /* so we return f*rho  CalcRho will do f/rho and rho = A/f */
 return f;

}/* AnumHulthen */

double Glauber::AnumHulthenInt(double xi)
{
 /* this is dummy */
 
 return 0.0;

}/* AnumHulthenInt */

double Glauber::NuIntHulthen(double xi)
{
 double f, g;
 double z, r, s;
 double b_WS, a_WS, rho;

 a_WS = Nuc_WS->a_WS;
 b_WS = Nuc_WS->w_WS;
 // R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 
   
/* multiply by a_WS (in fm^-1), make life simpler */
   
   s = NuInS_S*a_WS;
   z = -log(xi);
   r = sqrt(s*s + z*z);

   /* mult by 2 because the integral is originally over -infty to infty */

   f = 2.0*a_WS*rho*(GlauberData.SigmaNN);
   g = (1.0/r)*(exp(-r) - exp(-(b_WS/a_WS)*r));
   f *= g*g;

   //   fprintf(stderr, "%e\n", f);
 return f;
}/* NuIntHulthen */


double Glauber::integral (int id, double down, double up, double tol, int *count)
{
  double dx, y, g1[7];
  int i;
  
  if (down == up) y = 0.0;
  else
    {
      dx = (up-down)/6.0;
      for( i=0; i<7; i++) 
	{
	  if (id==1) g1[i] = NuInt2HO(down + i*dx);
	  else if (id==2) g1[i] = NuInt3Gauss(down + i*dx);
	  else if (id==3) g1[i] = NuInt3Fermi(down + i*dx);
	  else if (id==4) g1[i] = Anum3FermiInt(down + i*dx);
	  else if (id==5) g1[i] = Anum3GaussInt(down + i*dx);
	  else if (id==6) g1[i] = Anum2HOInt(down + i*dx);
	  else if (id==7) g1[i] = OLSIntegrand(down + i*dx);
	  else if (id==8) g1[i] = NuIntHulthen(down + i*dx);
	  //
	  //cout << *count << " " << id << " " << down << " " << up << ", g1[" << i << "]=" << g1[i] << endl;
	}
      *count = 7;
      y = qnc7(id, tol, down, dx, g1, 0.0, 0.0, count);
    }
  return y;
} /* end of integral */

double Glauber::qnc7(int id, double tol, double down, double dx, double *f_of, 
		     double pre_sum, double area, int *count)
{
  int i;
  double left_sum, right_sum, ans;
  static double w[] = 
    {41.0/140.0, 54.0/35.0, 27.0/140.0, 68.0/35.0, 27.0/140, 54.0/35.0,
       41.0/140.0};
  double fl[7];
  double fr[7];
  /*
    qnc7 calculates integral over left and right half of the given interval
    and branches
    to do so, first halve dx
    */
  
  dx /= 2.0;


  /*
    first calculate the left estimate
    f_of[] contains the evaluated values at down+i, 0< i <7
    store half distanced values for the left sum in fl[]
    */

  
  if (id==1)
    {
      fl[1] = NuInt2HO(down + dx);
      fl[3] = NuInt2HO(down + 3.0*dx);
      fl[5] = NuInt2HO(down + 5.0*dx);
    }
  else if (id==2)
    {
      fl[1] = NuInt3Gauss(down + dx);
      fl[3] = NuInt3Gauss(down + 3.0*dx);
      fl[5] = NuInt3Gauss(down + 5.0*dx);
    }
  else if (id==3)
    {
      fl[1] = NuInt3Fermi(down + dx);
      fl[3] = NuInt3Fermi(down + 3.0*dx);
      fl[5] = NuInt3Fermi(down + 5.0*dx);
    }
  else if (id==4)
    {
      fl[1] = Anum3FermiInt(down + dx);
      fl[3] = Anum3FermiInt(down + 3.0*dx);
      fl[5] = Anum3FermiInt(down + 5.0*dx);
    }  
  else if (id==5)
    {
      fl[1] = Anum3GaussInt(down + dx);
      fl[3] = Anum3GaussInt(down + 3.0*dx);
      fl[5] = Anum3GaussInt(down + 5.0*dx);
    }  
  else if (id==6)
    {
      fl[1] = Anum2HOInt(down + dx);
      fl[3] = Anum2HOInt(down + 3.0*dx);
      fl[5] = Anum2HOInt(down + 5.0*dx);
    }  
  else if (id==7)
    {
      fl[1] = OLSIntegrand(down + dx);
      fl[3] = OLSIntegrand(down + 3.0*dx);
      fl[5] = OLSIntegrand(down + 5.0*dx);
    }  
  else if (id==8)
    {
      fl[1] = NuIntHulthen(down + dx);
      fl[3] = NuIntHulthen(down + 3.0*dx);
      fl[5] = NuIntHulthen(down + 5.0*dx);
    }  
  
  fl[0] = f_of[0];
  fl[2] = f_of[1];
  fl[4] = f_of[2];
  fl[6] = f_of[3];
  
  
  *count += 3;


  left_sum = 0.0;
  for(i=0; i<7; i++) left_sum += w[i]*fl[i];
  left_sum *= dx;

/*printf("leftsum is %le\n", left_sum);*/


  /*
    like wise, the right sum is in fr[]
    */

  if (id==1)
    {
      fr[1] = NuInt2HO(down + 7.0*dx);
      fr[3] = NuInt2HO(down + 9.0*dx);
      fr[5] = NuInt2HO(down + 11.0*dx);
    }
  else if (id==2)
    {
      fr[1] = NuInt3Gauss(down + 7.0*dx);
      fr[3] = NuInt3Gauss(down + 9.0*dx);
      fr[5] = NuInt3Gauss(down + 11.0*dx);
    }
  else if (id==3)
    {
      fr[1] = NuInt3Fermi(down + 7.0*dx);
      fr[3] = NuInt3Fermi(down + 9.0*dx);
      fr[5] = NuInt3Fermi(down + 11.0*dx);
    }
  else if (id==4)
    {
      fr[1] = Anum3FermiInt(down + 7.0*dx);
      fr[3] = Anum3FermiInt(down + 9.0*dx);
      fr[5] = Anum3FermiInt(down + 11.0*dx);
    }
  else if (id==5)
    {
      fr[1] = Anum3GaussInt(down + 7.0*dx);
      fr[3] = Anum3GaussInt(down + 9.0*dx);
      fr[5] = Anum3GaussInt(down + 11.0*dx);
    }
  else if (id==6)
    {
      fr[1] = Anum2HOInt(down + 7.0*dx);
      fr[3] = Anum2HOInt(down + 9.0*dx);
      fr[5] = Anum2HOInt(down + 11.0*dx);
    }
  else if (id==7)
    {
      fr[1] = OLSIntegrand(down + 7.0*dx);
      fr[3] = OLSIntegrand(down + 9.0*dx);
      fr[5] = OLSIntegrand(down + 11.0*dx);
    }
  else if (id==8)
    {
      fr[1] = NuIntHulthen(down + 7.0*dx);
      fr[3] = NuIntHulthen(down + 9.0*dx);
      fr[5] = NuIntHulthen(down + 11.0*dx);
    }  
  
  fr[0] = f_of[3];
  fr[2] = f_of[4];
  fr[4] = f_of[5];
  fr[6] = f_of[6];

  
  *count += 3;

  right_sum = 0.0;
  
  for(i=0; i<7; i++) 
    {
      right_sum += w[i]*fr[i];
    }
  right_sum *= dx;
  

  ans = left_sum + right_sum;

  area += -fabs(pre_sum) + fabs(left_sum) + fabs(right_sum);


  if( fabs(ans - pre_sum) > tol*fabs(area) && (*count < limit))
    {
      /*
	branch by calling the function itself
	by calling the qnc7 twice, we are branching
	since left hand side is being calculated first, until the condition
	is satisfied, the left branch keeps branching
	when finally the condition is met by one left-most interval, 
	qnc7 returns the right hand side of one up level,
	and the same process resumes 
	until the criterion is met by all the branched
	intervals,
	then qnc7 returns to the original right branch and resumes halving 
	until the condition is met by all intervals
	(funk, down, dx, f_of[7], pre_ans, ans)
	*/

      tol /= 1.414213562;
      left_sum = qnc7(id, tol, down, dx, fl, left_sum, area, count);
      right_sum = qnc7(id, tol, down+dx*6., dx, fr, right_sum, area, count);

      ans = left_sum + right_sum;

      /* printf("ans is %le\n", ans);*/

    }/* belongs to if*/

  return ans;
} /* end of qnc */


double Glauber::OLSIntegrand(double s)
{
  double sum, arg, x, r;
  int k, m;
  m = 20;
  sum = 0.0;
  for(k=1; k<=m; k++)
    {
      arg = M_PI*(2.0*k - 1.0)/(2.0*m);
      x = cos(arg);
      r = sqrt(s*s + b*b + 2.0*s*b*x);
      sum += InterNuTInST(r);
    }/* k */
  
  return s*sum*M_PI/(1.0*m)*InterNuPInSP(s);
  
}/* OLSIntegrand */


double Glauber::TAB()
{
  double f;
  int count = 0;
  f = integral(7, 0.0, GlauberData.SCutOff, TOL, &count); // integrate OLSIntegrand(s)
  f *= 2.0/(GlauberData.SigmaNN); //here TAB is the number of binary collisions, dimensionless 
                                //(1/fm^4 integrated over dr_T^2 (gets rid of 1/fm^2), divided by sigma (gets rid of the other))
  return f;
}/* TAB */


double Glauber::PAB(double x, double y)
{
  double s1=sqrt(pow(x+b/2.,2.)+y*y);
  double s2=sqrt(pow(x-b/2.,2.)+y*y);
  return InterNuPInSP(s1)*InterNuTInST(s2)/(currentTAB*GlauberData.SigmaNN);
}/* PAB */


void Glauber::initGlauber(double SigmaNN, string Target, string Projectile, double inb, int imax, int rank)
{
  string Target_Name;
  Target_Name = Target;

  string Projectile_Name;
  Projectile_Name = Projectile;
 
  string p_name;
  stringstream sp_name;

  const char* EOSPATH = "HYDROPROGRAMPATH";
  char* envPath = getenv(EOSPATH);
  if (envPath != 0 && *envPath != '\0') 
    {
      sp_name << envPath;
      sp_name << "/known_nuclei.in";
      p_name = sp_name.str();
    }  
  else
    {
      sp_name << "./known_nuclei.in";
      p_name = sp_name.str();
    }

  string paf;
  paf = p_name;

  FindNucleusData2(&(GlauberData.Target), Target_Name, rank);
  FindNucleusData2(&(GlauberData.Projectile), Projectile_Name, rank);
 
  GlauberData.SigmaNN = 0.1*SigmaNN; // sigma in fm^2 
  currentA1 = GlauberData.Projectile.A;
  currentA2 = GlauberData.Target.A;
  currentZ1 = GlauberData.Projectile.Z;
  currentZ2 = GlauberData.Target.Z;  
  
  GlauberData.InterMax = imax;
  GlauberData.SCutOff = 12.;
  
  b=inb; 
}/* init */

double Glauber::areaTA(double x, double A)
{
  double f;
  f = A*220.*(1.-exp(-0.025*x*x)); 
  return f;
}

ReturnValue Glauber::SampleTARejection(Random *random, int PorT)
{
  ReturnValue returnVec;
  
  double r, x, y, tmp;
  double phi;
  double A=1.2*GlauberData.SigmaNN/4.21325504715; // increase the envelope for larger sigma_inel (larger root_s) (was originally written
  // for root(s)=200 GeV, hence the cross section of 4.21325504715 fm^2 (=42.13 mb)
  cout.precision(10);
  if(PorT == 1)
    {
      do
	{
	  phi = 2.*M_PI*random->genrand64_real1();
	  r = 6.32456*sqrt(-log((-0.00454545*(-220.*A+areaTA(15.,A)*random->genrand64_real1()))/A));
	  // here random->genrand64_real1()*areaTA(GlauberData.SCutOff) 
	  // is a uniform random number on [0, area under f(x)]
	  tmp = random->genrand64_real1();

	  // x is uniform on [0,1]
	  if ( r*InterNuPInSP(r) > A*r*11.*exp(-r*r/40.) ) 
	    cout << "WARNING: TA>envelope: " << "TA=" << r*InterNuPInSP(r) 
		 << ", f=" << A*r*11.*exp(-r*r/40.) << endl;
	} while( tmp > r*InterNuPInSP(r)/(A*r*11.*exp(-r*r/40.))); 
    }
  else
    { 
      do
	{
	  phi = 2.*M_PI*random->genrand64_real1();
	  r = 6.32456*sqrt(-log((-0.00454545*(-220.*A+areaTA(15.,A)*random->genrand64_real1()))/A));
	  // here random->genrand64_real1()*areaTA(GlauberData.SCutOff) 
	  // is a uniform random number on [0, area under f(x)]
	  tmp = random->genrand64_real1();
	  // x is uniform on [0,1]
	  if ( r*InterNuTInST(r) > A*r*11.*exp(-r*r/40.) ) 
	    cout << "WARNING: TA>envelope: " << "TA=" << r*InterNuTInST(r) 
		 << ", f=" << A*r*11.*exp(-r*r/40.) << endl;
	} while( tmp > r*InterNuTInST(r)/(A*r*11.*exp(-r*r/40.)));
    } 
  // reject if tmp is larger than the ratio p(y)/f(y), f(y)=A*r*11.*exp(-r*r/40.))
  x=r*cos(phi);
  y=r*sin(phi);
  returnVec.x=x;
  returnVec.y=y;
  returnVec.collided=0;
  return returnVec; 
}

