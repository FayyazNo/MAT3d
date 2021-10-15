#include "FFT_analysis3D.hpp"
#include<fstream>


// Main Program
int main(int argc, char *argv[]) 
{
	ifstream inFile;
	string input;
	string file=argv[1];

	double in[7];
	string in2[6];

	inFile.open(file.c_str());
	if(inFile.fail())
	{
		cerr << "Error while opening the file, the "<<in2[0]<< " does not exist "  <<endl;	
		exit(EXIT_FAILURE);
	}
		
	int s=-1, c=-1;
	while(inFile>>input){
		inFile>> input;
		if(atof(input.c_str())){

			s+=1;
			if (s==7) {s-=1; break;}
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
	double E1[6];
	for (auto i=0; i<6; ++i){

		E1[i] = atof(input.c_str());
		int y=input.find(",");
		input.erase(input.begin(), input.begin()+y+1);

	}
	cout<<E1[0]<<"  "<<E1[1]<<"  "<<E1[2]<<"  "<<E1[3]<<"  "<<E1[4]<<"  "<<E1[5]<<"  "<<std::endl;
	//---------------------------------------------------------------------------------
	inFile.close();

	if(in2[1]!="Y" and in2[1]!="N")
	{cerr<<"Just 'Y' and 'N' are allowed character for vtiOutputYorN!"<<endl; exit(EXIT_FAILURE);}

	if(in2[3]!="S" and in2[3]!="H")
	{cerr<<"Just 'S' and 'H' are allowed character for SingleAnalysisOrHomogenization!"<<endl; exit(EXIT_FAILURE);}

	if(c!=5 and s!=6)
	{cerr<<"The input data is incomplete or wrong!"<<endl; exit(EXIT_FAILURE);}
	//--------------------------------------------

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

    MatrixDdd E(6, 1);                            //Predefined strain vector which will be applied to all voxels
    E<<E1[0], E1[1], E1[2],E1[3], E1[4], E1[5];
/*
    double q [4];                                  // Mechanical properties of the sample
    q[0]=in[1];                                      // Elastic modules of white phase (or 1)
    q[1]=in[2];                                       // Poisson's ratio of white phase (or 1)
    q[2]=in[3];                                     // Elastic modules of Black phase (or 0)
    q[3]=in[4];                                      // Poisson ratio of   Black phase (or 0)

*/



int noPhases=2;
	MatrixDdd mat(noPhases, 7);
	//  E, v                                                    Elastic Constants;
	//  sigma= s_y+H(epsilon_plastic)^n                         Plastic Constitutive
	//  sigma=(1-d)epsilon, d=exp(P*Y-Y0), Y=epsilon:C:epsilon  Damage  Constitutive

	//   E           v      s_y          H       n      P       Y0     
	mat<<  in[1],   in[2],   40e6,     100e6,     1,   8e-2,   8e2,    //material 0
	       in[3],   in[4],  250e6,      8006,     .2,   8e-4,   8e3;    //material 1

	MatrixDdd X(6,1);
	X<<0,1,0,0,0,0;
//******************************************************************************
	string mater_behav="ep";
	FFT_Analysis3D test(Nx, Ny,Nz, mater_behav);                 // Construct an analysis plan with assignment of the RVE size


    test.material(mat, RVE, 2);       // Set material properties and characteristic function, no of phases, and material  //
                                                    //behavior

    test.setLoad(E, X);                             // Set strain vector to all voxels

	if (in2[3]=="S")
		test.runAnalysis(in[5], in[6], 15);               // Run a single analysis and try to catch equlibrium state, two argument are error  
	else                                                //  tolerances.
		test.homogenization(in[5], in[6], in2[4]);

	cout<<in2[1]<<endl;

	 if (in2[1]=="Y")
		test.vtiWriter(in2[2]);

	//test.printer();
    return 0;
}
