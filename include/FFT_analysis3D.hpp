#ifndef FFT_ANALYSIS3D_HPP
#define FFT_ANALYSIS3D_HPP

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


typedef Matrix<double, Dynamic, Dynamic> MatrixDdd;
typedef Matrix<int, Dynamic, Dynamic> MatrixIdd;



class FFT_Analysis3D 
{
    private:
// two temperary 1D array for fft and ifft
        fftw_complex  *t2, *t3, *t4; 
		double *t1;
// two plan for fft and ifft
        fftw_plan p_t1_t2_frw,  p_t1_t4_frw,  p_t3_t1_bkw;
// characteristic function of the RVE which should has 0(for inclusion)
//and 1 (for matrix)
        MatrixIdd RVE;

// 1D array of Green operator
        double *GF;
// iteration counter in each analysis
        int itNo;
// convergence parameter
        double eqVar;
// the problem can be solved based on the two methods:
// 1- Strain based method which is approprate when inclusion is soft (Ei<Em)
// 2- Stress Base method which is approprate when inclusion is rigid (Ei>Em)
// after assigning the material properties the program choose the approprate method
        bool strainBase;
        bool stressBase;

// sum of the sigma in all voxels in current and previous iteration and is used in calculation 
        double sum_field_old;
		double sum_field;
        double diff_sum_field;

// Equlibrum parameter in current and previous iteration

		double avr_eqVar;
		double aver_eqVarOld;
		double diff_avr_eqVar;
// fftw3 parameters		
		fftw_iodim64 *dims, *dims_o, *howmany_dims;

// parameter for E_update function
		MatrixDdd sigma_iter_old;
		MatrixDdd epsilon_iter_old;
//Number of increments
		double nInc; 
//Identity tensor in vector form (i.e <1,1,1,0,0,0>)
        MatrixDdd I; 
//      material behavior e:elastic, ep:elastoplastic, ed:elastodamage
		string mater_behav;
      
		double incNo;
    public:

// stiffness tensors of the inclusion, matrix and refrence materials 
        MatrixDdd *C;
        MatrixDdd C0;

		MatrixDdd *C_acc;
		MatrixDdd *S_acc;
//compliance tensors of the inclusion, matrix and refrence materials 
        MatrixDdd *S;
        MatrixDdd S0;

// strain. stress and polarization tensors
        MatrixDdd epsilon;
        MatrixDdd sigma;

// predefined load, uniform strain in stran based methed and uniform stress
// in stress based methos
        MatrixDdd E;
		MatrixDdd dE;
//predefined stress direction look at the Mulinec and Suquet paper
        MatrixDdd X0;
// number of the voxels in each direction
        int Nx, Ny, Nz;
// xyz=NX*Ny*Nz
        int xyz;
// Lamme constants of the reference material,
        double la0;
        double mu0;
        double *la;
        double *mu;
        double *El;
        double *v;
		double *K;

//  sigma= S_y+H(epsilon_plastic)^n    Plastic Constitutive
		double *H;
		double *n;
		double *Sy;

//  sigma=(1-d)epsilon, d=exp(P*Y-Y0), Y=epsilon:C:epsilon  Damage Constitutive
		double *P;
		double *Y0;
		int noPhases;
		
		double E11, E22, E33;

// Homogenized strain= <epsilon>
        MatrixDdd epHomo;
// Homogenized stess= <sigma>
        MatrixDdd sigHomo;
// Homogenized compliance matrix
        MatrixDdd SHomo;
// Homogenized stiffness matrix
        MatrixDdd CHomo;  
// volume fractions;
		double *vol_fractions;

		MatrixDdd epsilon_old;
		MatrixDdd sigma_old;
       // MatrixDdd epsilon_p;
        MatrixDdd epsilon_p_eff;
        MatrixDdd sigma_eff;
        MatrixDdd sigma_d;
		MatrixDdd phi;
		MatrixDdd d_ep_p;
//---------------
		MatrixDdd d;
		MatrixDdd d_old;	
		MatrixDdd Y;
		
		MatrixDdd prop_homogenized;	

    private:
// Assigne Green operator to GF array
        void setGF();
// Multiplying Cm or Ci in epsilon:  sigma=C(x)*epsilon in strainBased method
// Multiplying Sm or Si in sigma :   epsilon=S(x)*sigma in stressBased method 
        void multipCx();
// Fourier transform of the polrization field
        void FFT_sigma();
        void FFT_epsilon();
// Multiplying Green operator in FFT(tau)
        void multiplyGF();
// IFFT
        void FFTBCKWRD_sigma();
        void FFTBCKWRD_epsilon();
//Checking convegence
        void convergenceCheck(); 

//  update fields in accelerated scheme
        void update(); 
		void elastoDamage();
		void elastoPlastic();

//update E
		void updateE();
		void old_iter_fields_update();
		void export_homogen_var(int _i);
		void update_old_Variable();
		void constitutive();
		void sigma_eff_calculator(int _i);

		double nrm2 (double *zet);
		double dlt  (int i, int j); 
		void max_min_Array(const double  *array, int length, double out[4]);
		// inner product of a and b
		double prod(MatrixDdd &a, MatrixDdd &b);

    public:
// constructor, inputs are number of voxels in 3 directions
        FFT_Analysis3D(int Nx1, int Ny1, int Nz1, string mater_behav);
//Destructur
        ~FFT_Analysis3D();
// Assigning Material properties, porperties=[Em, vm, Ei, vi], RVE: characterstic function of the RVR
        void material (MatrixDdd &mat, MatrixIdd &RVE1, int noPhases);
// Assigne pre defined Load E1=[E11, E22, E33, E31, E32, E12]
        void setLoad(MatrixDdd &E1, MatrixDdd &X0);
// start a single Analysis, arguments are error tolerances it is recommended to use 1e-9 and 1e-6 respectively
// tol_field means average dimensionless error of the stress in two successive iterations
// tol_eq means difference of the equlibrium paarameter for two successive iterations
        void runAnalysis(double tol_field, double tol_eq, int nInc1);
// start a homogenizatio analysis (include 6 single analysis), first two arguments are 
//error tolerances it is recommended to use 1e-9 and 1e-6 respectively, third argument is output file name
        double homogenization(double, double, string);
// gives the output in the folder results
        void printer (std::string outFldr="results");
// gives vti file, input is file name, example: "Myfile.vti"
        void vtiWriter(std::string );
    
};
#endif
