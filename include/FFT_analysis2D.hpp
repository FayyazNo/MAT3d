#include <iostream>
#include <eigen/Eigen/Dense>
#include "fftw3.h"
#include <CImg.h>
#include <fstream>

#ifndef FFT_ANALYSIS2D_HPP
#define FFT_ANALYSIS2D_HPP
using namespace std;
using namespace Eigen;
using namespace cimg_library;


typedef Matrix<complex<double>, Dynamic, Dynamic> MatrixCdd;
typedef Matrix<double, Dynamic, Dynamic> MatrixDdd;
typedef Matrix<int, Dynamic, Dynamic> MatrixIdd;

inline double nrm2 (double *zet) {return (zet[0]*zet[0]+zet[1]*zet[1]);}
inline double dlt  (int i, int j) {if (i==j) return 1.; else return 0.;}

class FFT_Analysis2D
{
  private:
      fftw_complex *t1, *t2; 
      fftw_plan p3, p4;
      MatrixDdd RVE;
      double *GF;
      int itNo=0;
      double convrg=0;

  public:
      MatrixDdd Ci;
      MatrixDdd Cm;
      MatrixDdd C0;
      MatrixCdd epsilon;
      MatrixCdd sigma;
      MatrixCdd tau;
      MatrixCdd sigmaOld;
      MatrixCdd E;
      MatrixDdd epHomo;
      MatrixDdd CHomo;
      MatrixDdd SHomo;
      MatrixDdd Si;
      MatrixDdd Sm;
      MatrixDdd sigHomo;
      int Nx, Ny;
      int xy;
      double la0;
      double mu0;
      double Em;
      double vm;
      double Ei;
      double vi;
      double mum;
      double lam;
      double mui;
      double lai;
      
  private:
      // Creat Green operator and stor it in 1D array for all voxells
      // each voxels has just 6 independent elements 
      void setGF(); 
      // perform sigma(x)=C(x)*esilon(x) where x can be each voxel 
      void multipCx();
      
      // FFT ot tau where tau= (C(x)-C0)*epsilon(x)
      //First stors tau in a 1D array t1 and after FFT tranform
      //stores FFT(t1) in t2
      void FFT_tau();
      
      // Perform t1=-GF*t2
      void multiplyGF();
      // Perform t2=FFT_backward (t1)
      // copy t2 in epsilon
      void FFTBCKWRD_epsilon();
      
      // Check equilibrium state and diffrence of sigma in two successive iteration
      void convergenceCheck();
      
  public:
      // Constructor, Nx1 and Ny1 are dimension of RVE
      FFT_Analysis2D(int Nx1, int Ny1);
      // Destructor
      ~FFT_Analysis2D();
      //Assigne material properties, properties[0] and properties [1]
      // are E and v matrix and properties [3] and properties [4] are 
      // E and v inclusion. RVE1 is characterstic function where contain 0 and 1
      // 1 for matrix phase and 0 for inclusion phase.
      void material (double *properties, MatrixDdd &RVE1);
  
      // E1 is strain vector [E11, E22, E12] applied to all voxels
      void setLoad(MatrixCdd &E1);
      //Start analysis ant tray to catch equilibrated state after applying E1 to all
      // voxels
      void runAnalysis();
      //Run 3 Analysis with 6 diffrent E1 and gives homogenized elastic tensor
      void homogenization();
      // Export Results    
      void printer();
};
#endif
