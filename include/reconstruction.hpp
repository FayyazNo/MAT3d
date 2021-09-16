#ifndef RECCLASS_HPP
#define RECCLASS_HPP

#include <iostream>
#include <eigen/Eigen/Dense>
#include "fftw3.h"
#include "CImg.h"
#include <fstream>
#include "omp.h"


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkRenderWindow.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
using namespace std;
using namespace Eigen;
using namespace cimg_library;

typedef Matrix<complex<double>, Dynamic, Dynamic> MatrixCdd;
typedef Matrix<double, Dynamic, Dynamic> MatrixDdd;
typedef Matrix<int, Dynamic, Dynamic> MatrixIdd;

namespace recon
{
  void  fftshift(double *, double *, int , int , int );
  void ifftshift(double *, double *, int , int , int );
  
  inline double funcNormC (fftw_complex *MTilda, double *normC, int xyz) 
  {
     double err=0;
     for (int k=0; k<xyz; ++k)
         err+=abs(sqrt(MTilda[k][0]*MTilda[k][0]+MTilda[k][1]*MTilda[k][1])-normC[k]);
     err=err/xyz;
     return err;
  }

  class recPlan {
	public:
	   /* Constructors */
	   /* Using c1 and c2 images where c1 and c2 are 2 perpenicular cut section of TPCF*/
	   recPlan(CImg<unsigned char> c1, CImg<unsigned char> c2);
	   /* Using TPCF values and [Nx1, Ny1, Nz1] where  TPCF is full set of
		two point correlation function the RVE (TPCF=(1/Nx/Ny/Nz)*(norm(FFT3(RVE)^)2)
		[Nx1, Ny1, Nz1] are dimension of RVE
		*/
	   recPlan(double *TPCF, int Nx1, int Ny1, int Nz1);
	   // destructor
	   ~recPlan();

	   /* Helpper to estimate norm of the FFT of RVE */
	   void estimateNormFFT_RVE();
	   /* Runs reconstruction iterations where n, b are number of iteration and beta 
	   respectively. beta is a parameter used in phase recovery algorithm and
		0<beta<1
	   
	   */
	   void runRec(int n, double b);
	   /* Saves ouput cross sections to the locFldr folder */
	   void saveRecImages(string locFldr);
	   //vtk format writer
	   void vtiWriter(string fileName, double*);
	
    public:
	   MatrixDdd trackrealerrorA;
	   fftw_complex *mTilda;
	   void  fftShift(double*, double*, int, int, int);
	   void ifftShift(double*, double*, int, int, int);

	private:
	   int xyz;
	   fftw_complex  *MHat, *mHat, *MTilda, *CC1, *CC2, *CC3, *CC4;
	   double *normC;
	   MatrixDdd RVE;
	   MatrixCdd temp;
	   MatrixDdd ang ;
	   MatrixDdd ang1;
	   MatrixDdd ang2;
	   fftw_plan p1, p2, p3, p4, p5, p6, p7, p8;
	   int Nx, Ny, Nz, Nzz;
	   int numloops=0;
	   double beta;
	   int cd;
    };
}
#endif
