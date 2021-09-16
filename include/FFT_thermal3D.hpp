#ifndef FFT_THERMAL3D_HPP
#define FFT_THERMAL3D_HPP

#include <iostream>
#include <eigen/Eigen/Dense>
#include "fftw3.h"
#include <CImg.h>
#include <fstream>
#include <ctime> 
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



class FFT_analysis3d_thermal
{
    private:
// two temperary matrixfor fft and ifft calculations
        fftw_complex *t3, *t2; 
		double *t1;
// two plan for fft and ifft respectively
        fftw_plan p3, p4;
// charactestic function of the RVE which just hase 0 and 1 elements
// for inclusion (0) and for matrix(1)
        MatrixIdd RVE;
// Green operator matrix
        double *GF;
// iteration number
        int itNo;
		fftw_iodim64 *dims, *dims_o;
		fftw_iodim64 *howmany_dims;

// sum of the sigma in all voxels in current and previous iteration and is used in calculation 
        double sum_field_old;
		double sum_field;
        double diff_sum_field;

// Equlibrum parameter in current and previous iteration

		double avr_eqVar;
		double aver_eqVarOld;
		double diff_avr_eqVar;
		
/* There are two methods to solve the problem
   1-Heat flux based: which is suitable when Inclusion conductivity is greater than matrix conductivity (ki>km)
   2- Temperatur gradient Base which is well suited for inclusion with lower thermal conductivity (km>ki)
   Based on the inouts the program choose the aproprate method
*/
        bool heatFluxBase;
        bool tempGradientBase;
  
    public:
// thermal conductivity tensor of inclusion
        MatrixDdd Ki;
// thermal conductivity tensor of matrix
        MatrixDdd Km;
// thermal conductivity tensor of reference phase K0=alf*(Km+Ki)/2
        MatrixDdd K0;
// teperature gradient
        MatrixDdd tempGradient;
//
		MatrixDdd field_c;
//
	 	MatrixDdd K1;
		MatrixDdd K2;

	 	MatrixDdd R1;
		MatrixDdd R2;
// heat flux
        MatrixDdd heatFlux;

//polarization fiel, tau= heatFlux-K0*tempGradient
        MatrixDdd tau;
//predefined load
        MatrixDdd E;
        MatrixDdd tempGradHomo;
//homogenized thermal conductivity tensor
        MatrixDdd KHomo;
//homogenized thermal Resistivity tensor
        MatrixDdd RHomo;
// resistivity tensor of inclusion
        MatrixDdd Ri;
// resistivity tensor of matrix
        MatrixDdd Rm;
// resistivity tensor of reference material
        MatrixDdd R0;
        MatrixDdd heatFluxHomo;
//Number of voxels in each direction
        int Nx, Ny, Nz;
//xyz=Nx*Ny*Nz
        int xyz;

	double K11, K22, K33;

    private:
// Calculate and assigne Green operator to GF array
        void setGF();
//Multiply Km or Ki in tempGradient heatFlux=K(x)*tempGradient
        void multipKx();
// Calculating FFT
        void FFT();
// Multiplying Green operator in FFT(tau)
        void multiplyGF();
//FFT backward
        void FFTBCKWRD();
        void convergenceCheck();
        void strainUpdate();

    public:
// constructor of the class which take pixel in each direction
        FFT_analysis3d_thermal (int Nx1, int Ny1, int Nz1);
//Destracture of the Class
        ~FFT_analysis3d_thermal();
// Assigne material properties, km thermal conductivity of matix, ki thermal conductivity inclusion
// RVE characterstic function of RVE
        void material (double km, double ki, MatrixIdd &RVE1);
//Assine predefined load
        void setLoad(MatrixDdd &E1);

// start a single Analysis, arguments are error tolerances it is recommended to use 1e-6 and 1e-5 respectively
// tol_field means average dimensionless error of the stress in two successive iterations
// tol_eq means difference of the equlibrium paarameter for two successive iterations
        void runAnalysis(double tol_field, double tol_eq );
// start a homogenizatio analysis (include 3 single analysis), first two arguments are 
//error tolerances it is recommended to use 1e-9 and 1e-6 respectively, third argument is output file name
        void homogenization(double tol_field, double tol_eq, string filename );
//gives outputs
        void printer (string);
// vti file writer, the argument is file name;		
		void vtiWriter(std::string );

};

#endif

