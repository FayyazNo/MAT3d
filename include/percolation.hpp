#include <iostream>
#include <eigen/Eigen/Dense>
#include "fftw3.h"
#include <CImg.h>
#include <fstream>
#include <algorithm> 
#include <vector>

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

#ifndef PERCOLATION_HPP
#define PERCOLATION_HPP


typedef Matrix<complex<double>, Dynamic, Dynamic> MatrixCdd;
typedef Matrix<double, Dynamic, Dynamic> MatrixDdd;
typedef Matrix<int, Dynamic, Dynamic> MatrixIdd;

//---------------------------------------------------
// TODO : Add comments for the methods
class percolation 
{
  public:
// finding cluster of the phase intrested phase
	percolation (MatrixDdd & RVE, int Nx1, int Ny1, int Nz1, int noOfPhases, int intrestPhase);
	~percolation();
// increasing total_iter until output become constant (clusters data), 8 is recommended for rve with concave shapes
// for seperated convex shape 1 is enough
	void runPercolationAnalysis(int total_iter);
	int commonSurface;
//calculate common surface between phase1 and phase2
	void SurfaceOfTwoPhase (MatrixDdd &RVE, int phase1, int phase2);
//
// fold is name if the folder to store results, clustresSizedata is string name od the txt file 
// includes clusters size data, largestCluster is name of the RVE which include the largest cluster
// of the intrested phase
	void printer(std::string fold, std::string clustresSizedata, std::string largestCluster);
	int NX;
	int NY;
	int NZ;
	int noOfPhases;
	int intrestPhase;
	void vtiWriter(std::string fileName);
	vector<int>	all_clusters;
   	vector<double>	clusters_vol_frac; 
	double vol_intrstd_phs;  
  private:
	void setMargin(int);
	void neighborCheck(int i, int j, int k);
	void removeMargin();
	int Nx;
	int Ny;
	int Nz;
	int *C, *C1;
	int neighbers[3];
	int phasePlus;
	MatrixIdd intrestedNeigbor;
};
#endif
