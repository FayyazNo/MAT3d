#include "FFT_analysis2D.hpp"
#include<fstream>


// Main Program
int main() 
{
    CImg<unsigned char> phaseImg("C2.tiff");        // Read a picture
    int Nx=phaseImg.height();
    int Ny =phaseImg.width ();
    cout<<"Nx="<<Nx<<"   NY="<<Ny<<endl;
    int Max=0;
    int xy=Nx*Ny;
    MatrixDdd RVE(1,xy);

    for (int i=0; i<Nx; ++i)                    // Transorm the picture to binary (0 or 1) array
      for(int j=0; j<Ny; ++j)
         Max=max(Max, (int)(phaseImg.atXY(i,j)));


    for (int i=0;  i<Nx; ++i)
      for(int j=0; j<Ny; ++j)
         RVE( j+i*Ny)=(int)(phaseImg.atXY(i,j))/Max;

    MatrixCdd E(3, 1);                       //Predefined strain vector
    E<<0, 0, 0.005;

    double q [4];                           // Material properties of the sample.
    q[0]=68.9e9;                            //Elastic modules of matrix
    q[1]=.35;                               // Poisson's ratio of matrix
    q[2]=400e9;                             // Elastic modules of inclusion
    q[3]=.23;                               // poisson's ratio of inclusion

    FFT_Analysis2D test(Nx, Ny);           
    test.material(q, RVE);
    test.setLoad(E);
   
    test.runAnalysis();
    test.printer();

    return 0;
}

