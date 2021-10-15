#include "FFT_thermal3D.hpp"
#include <fstream>
#include <stdlib.h>  

// Main Program
int main(int argc, char *argv[]) 
{

	ifstream inFile;
	string input;
	string file=argv[1];

	double in[5];
	string in2[6];

	inFile.open(file);
	if(inFile.fail())
	{
		cerr << "Error while opening the file, the "<<in2[0]<< " does not exist "  <<endl;	
		exit(EXIT_FAILURE);
	}
		
	int s=-1, c=-1;
	while(inFile>>input){
		inFile>> input;
		if(atof(input.c_str()) && s<5){
			s+=1;
			if (s==5) {s-=1; break;}
			in[s]=atof(input.c_str());

		cout<<"s="<<in[s]<<endl;
			
		}
		else{
			c+=1;
			in2[c]=input;

		cout<<"c="<<in2[c]<<endl;
		}
	}
	/*cout<<"ccc="<<c<<endl;
	cout<<"sss="<<s<<endl;*/
	double TG1[3];
	for (auto i=0; i<3; ++i){

		TG1[i] = atof(input.c_str());
		int y=input.find(",");
		input.erase(input.begin(), input.begin()+y+1);

	}
	inFile.close();

	if(in2[1]!="Y" and in2[1]!="N")
	{cerr<<"Just 'Y' and 'N' are allowed character for vtiOutputYorN!"<<endl; exit(EXIT_FAILURE);}

	if(in2[3]!="S" and in2[3]!="H")
	{cerr<<"Just 'S' and 'H' are allowed character for SingleAnalysisOrHomogenization!"<<endl; exit(EXIT_FAILURE);}

	//if(c!=5 and s!=6)
	//{cerr<<"The input data is incomplete or wrong!"<<endl; exit(EXIT_FAILURE);}
//---------------------------------------------------------------------------------
//------------------------------------------------------------------------------
    int Nx=in[0], Ny=in[0], Nz=in[0];
    int xyz=Nx*Ny*Nz;
    MatrixIdd RVE(1, xyz);  
    ifstream x (in2[0]);                      // Read the charcteristic function of the RVE as single coloumn .txt file
   	
    if(x.fail())
    {
        cerr << "Error while opening the file, the "<<in2[0]<< "does not exist "  <<endl;	
        exit(EXIT_FAILURE);
    }
    
    double a;
    int id=0;
    while (x>>a)
    {
        RVE(id)=a;
        id+=1;
    }
    if (id!=xyz)
	{
		cerr << "Dimension is not correct "  <<endl;	
		exit(EXIT_FAILURE);
	}

    MatrixDdd TG(3, 1);              //Predefined temperature gradient which will be applied to all voxels
    TG << 1, 0, 0;



    //----------------------------------------------------
	FFT_analysis3d_thermal test(Nx, Ny,Nz);      //Class constructor which take dimension of the RVE
    test.setLoad(TG);                            // Set uniform temperature gradient to all voxels
    test.material(in[1],in[2] , RVE);           // The first and second argument are km and ki, where 
                                                 // Ki=ki*I and Km=km*I and KI, Km and I are thermal conductivity 
                                                 //  tensors of inclusion, matrix and identity tensor respectively

	if (in2[3]=="S")
		test.runAnalysis(in[3], in[4]);                   // This function run a single analysis and try to find equilibrium state.
	else                                                  // The argument is number of the iterations.
		test.homogenization(in[3], in[4], in2[4]);       // This function run 3 analysis and gives effective thermal conductivity tensor
		                                                 // The argument is number of the iterations per analysis

	 if (in2[1]=="Y")
		 test.vtiWriter(in2[2]);                     // The argument is folder name.

    return 0;
}
