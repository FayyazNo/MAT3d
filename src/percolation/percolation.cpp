#include"percolation.hpp"
 
// Constructor Which take characterstic function as RVE, and its dimension as
// Nx, Ny and Nz. Nomber of phase (noPhases)   and the phase that
//percolation is calculated for (intrestPhase) are needed too.
     
percolation::percolation (MatrixDdd &RVE, int Nx1, int Ny1, int Nz1, int noOfPhases1, int intrestPhase1) 
{
    NX=Nx1;
    NY=Ny1;
    NZ=Nz1;
    Nx=NX+2;
    Ny=NY+2;
    Nz=NZ+2;

    noOfPhases=noOfPhases1;
    intrestPhase=intrestPhase1;
	vol_intrstd_phs=0;  

    neighbers [0]=-1;
    neighbers [1]=0;
    neighbers [2]=1;
    intrestedNeigbor.resize(4,9);
    phasePlus=noOfPhases+1;
    C1 = new int [NX*NY*NZ];
    C=   new int [Nx*Ny*Nz];
    for (int i=0; i<NX*NY*NZ; ++i)
        C1[i]=RVE(i);

}


//set -1 to all boundaries voxels and enumerate intresed phase voxels
// from (number_of_phases+1)
//--------------------------------------------------------
void percolation::setMargin(int cd)
{
    int p=phasePlus;
    for(int i=0; i<Nx; ++i)
        for (int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
            {
                if ((i==(Nx-1)) || (j==(Ny-1)) || (k==(Nz-1))|| (i==0)|| (j==0)|| (k==0))
                    C[k+j*Nz+i*Ny*Nz]=-1;
                else
                    {
                        if (C1[(k-1)+(j-1)*NZ+(i-1)*NY*NZ]==intrestPhase && cd==0)
                        {
                            C[k+j*Nz+i*Ny*Nz]=p;
                            p+=1;
							vol_intrstd_phs+=1;
                        }
                        else
                            C[k+j*Nz+i*Ny*Nz]=C1[(k-1)+(j-1)*NZ+(i-1)*NY*NZ];
                    }
            }

}

// remove -1 marginal voxels
//---------------------------------------------------------------
void percolation::removeMargin()
{
    for(int i=0; i<NX; ++i)
        for (int j=0; j<NY; ++j)
            for(int k=0; k<NZ; ++k)
                C1[k+j*NZ+i*NY*NZ]=C[(k+1)+(j+1)*Nz+(i+1)*Ny*Nz];	
}


// For  voxel (i, j, k) check all 9 neighbor voxels and for the voxel with
// intrested phase gives the minimum number of them to all other voxels
//----------------------------------------------------------------

void percolation::neighborCheck(int i, int j, int k)
{   
    int counter=0;
    intrestedNeigbor.setZero();
    for (int t=0; t<2; ++t)
        for(int r=0; r<2; ++r)
            for(int w=0; w<2; ++w)
            {
                int ii=i+neighbers[t];
                int jj=j+neighbers[r];
                int kk=k+neighbers[w];
                if ( C[kk+jj*Nz+ii*Ny*Nz]>=phasePlus)
                {
                    intrestedNeigbor(0,counter)=ii;
                    intrestedNeigbor(1,counter)=jj;
                    intrestedNeigbor(2,counter)=kk;
                    intrestedNeigbor(3,counter)=C[kk+jj*Nz+ii*Ny*Nz];
                    counter+=1;
                }
            }
    
    
    double Min=1e99;
    for (int i=0; i<counter; ++i)
        Min=min(Min, intrestedNeigbor(3,i)*1.);
    
    for (int i=0; i<counter; ++i)
    {
        int vv1=intrestedNeigbor(0,i);
        int vv2=intrestedNeigbor(1,i);
        int vv3=intrestedNeigbor(2,i);
        C[vv3+vv2*Nz+vv1*Ny*Nz]=Min;
    }

}


 // Start to check the neighbors for all voxels
// And perform the check in all possible direction that leads
// to very accurate and fast recogenition of all clusters.
//Just Works for Nx=Ny=Nz
//-----------------------------------------------------------------------------
void percolation::runPercolationAnalysis(int total_iter)
{
    setMargin(0);
for(int X=0; X<5; ++X)
{
    //1--------------------
	for(int t=1; t<7; ++t)
		for(int i=1; i<Nx-1; ++i)
		    for (int j=1; j<Ny-1; ++j)
		        for(int k=1; k<Nz-1; ++k)
		        {
					if (t==1  and C[k+j*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==2  and C[j+k*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==3  and C[k+i*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==4  and C[i+k*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==5  and C[i+j*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==6  and C[j+i*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);
		        }

//2--------------------
	for(int t=1; t<7; ++t)
		for(int i=1; i<Nx-1; ++i)
		    for (int j=1; j<Ny-1; ++j)
		        for(int k=Nz-2; k>0; --k)
		        {
					if (t==1  and C[k+j*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==2  and C[j+k*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==3  and C[k+i*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==4  and C[i+k*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==5  and C[i+j*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==6  and C[j+i*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);
		        }
	//3--------------------
	for(int t=1; t<7; ++t)
		for(int i=1; i<Nx-1; ++i)
		    for (int j=Ny-2; j>0; --j)
		        for(int k=1; k<Nz-1; ++k)
		        {
					if (t==1  and C[k+j*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==2  and C[j+k*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==3  and C[k+i*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==4  and C[i+k*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==5  and C[i+j*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==6  and C[j+i*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);
		        }
//4--------------------
	for(int t=1; t<7; ++t)
		for(int i=1; i<Nx-1; ++i)
		    for (int j=Ny-2; j>0; --j)
		        for(int k=Nz-2; k>0; --k)
		        {
					if (t==1  and C[k+j*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==2  and C[j+k*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==3  and C[k+i*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==4  and C[i+k*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==5  and C[i+j*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==6  and C[j+i*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);
		        }
//5--------------------
	for(int t=1; t<7; ++t)
		for(int i=Nx-2; i>0; --i)
		    for (int j=1; j<Ny-1; ++j)
		        for(int k=1; k<Nz-1; ++k)
		        {
					if (t==1  and C[k+j*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==2  and C[j+k*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==3  and C[k+i*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==4  and C[i+k*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==5  and C[i+j*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==6  and C[j+i*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);
		        }
//6--------------------
	for(int t=1; t<7; ++t)
		for(int i=Nx-2; i>0; --i)
		    for (int j=1; j<Ny-1; ++j)
		        for(int k=Nz-2; k>0; --k)
		        {
					if (t==1  and C[k+j*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==2  and C[j+k*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==3  and C[k+i*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==4  and C[i+k*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==5  and C[i+j*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==6  and C[j+i*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);
		        }
//7--------------------
	for(int t=1; t<7; ++t)
		for(int i=Nx-2; i>0; --i)
		    for (int j=Ny-2; j>0; --j)
		        for(int k=1; k<Nz-1; ++k)
		        {
					if (t==1  and C[k+j*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==2  and C[j+k*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==3  and C[k+i*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==4  and C[i+k*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==5  and C[i+j*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==6  and C[j+i*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);
		        }
//8--------------------
	for(int t=1; t<7; ++t)
		for(int i=Nx-2; i>0; --i)
		    for (int j=Ny-2; j>0; --j)
		        for(int k=Nz-2; k>0; --k)
		        {
					if (t==1  and C[k+j*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==2  and C[j+k*Nz+i*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==3  and C[k+i*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==4  and C[i+k*Nz+j*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==5  and C[i+j*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);

					else if (t==6  and C[j+i*Nz+k*Ny*Nz]>=phasePlus)
		                neighborCheck( i, j,  k);
		        }
}
    removeMargin();
	
//--------------------------------------------------------
//find all clusters
	vector<int>::iterator it;

	for (int i=0; i<NX*NY*NZ; ++i)
		if (C1[i]>=phasePlus)
		{
			if (all_clusters.size()<1)

				all_clusters.push_back(C1[i]);

			else 
			{
	 			it= std::find(all_clusters.begin(), all_clusters.end(), C1[i]);
				if (it==all_clusters.end())
					all_clusters.push_back(C1[i]);
			}
			
		}	

	//std::sort(all_clusters.begin(), all_clusters.end(), greater<int>());

//-----------------------------------------------------
	for (auto & elem : all_clusters)
	{	
		double s=0;
		for (int i=0; i<NX*NY*NZ; ++i)
			if (C1[i]==elem)
				s+=1;

   		clusters_vol_frac.push_back(s/vol_intrstd_phs); 
	}	

	vector<double> clusters_vol_frac_unsorted=clusters_vol_frac;
	std::sort(clusters_vol_frac.begin(), clusters_vol_frac.end(),  greater<double>());	
//-------------------------------------------------------------------------------
	int index1, index2;	
	vector<int>::iterator it1; 
	vector<double>::iterator it2; 
	
	for (int i=0; i<NX*NY*NZ; ++i)
		if (C1[i]>=phasePlus)
		{
 			it1= std::find(all_clusters.begin(), all_clusters.end(), C1[i]);
			index1 = std::distance(all_clusters.begin(), it1);


 			it2= std::find(clusters_vol_frac.begin(), clusters_vol_frac.end(), clusters_vol_frac_unsorted[index1]);
			index2 = std::distance(clusters_vol_frac.begin(), it2);
//to export largest cluster
/*			if (index2==0)
				C1[i]=1;

			else
				C1[i]=0;*/

//to export all clusters
			C1[i]=phasePlus+index2;

		}

	std::sort(clusters_vol_frac.begin(), clusters_vol_frac.end(),  greater<double>());	
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void percolation::SurfaceOfTwoPhase (MatrixDdd &RVE,int phase1, int phase2)
{
	int phas1=phase1, phas2=phase2;
	
	for (int i=0; i<NX*NY*NZ; ++i)
        C1[i]=RVE(i);

    setMargin(1);

	commonSurface=0;

	int ii, jj, kk, diff;
	for (int i=0; i<Nx; ++i)
		for (int j=0; j<Ny; ++j)
			for (int k=0; k<Nz; ++k)

				if ( i!=(Nx-1) && j!=(Ny-1) && k!=(Nz-1) && i*j*k !=0 &&  C[k+j*Nz+i*Ny*Nz]==phas1)
				{	
					for (int t=0; t<2; ++t)
						for(int r=0; r<2; ++r)
							for(int w=0; w<2; ++w)
							{
								ii=i+neighbers[t];
								jj=j+neighbers[r];
								kk=k+neighbers[w];
								diff=abs(ii-i)+abs(jj-j)+abs(kk-k);
								if ( diff==1 && C[kk+jj*Nz+ii*Ny*Nz]==phas2)
									commonSurface+=1;
							}
				}				

}
//--------------------------------------------------------------
void percolation::printer(std::string fold, std::string clustresSizedata, std::string largestCluster)
{
    std::cout << fold << "\n";
    int ss = mkdir(fold.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

    if (ss)
        std::cerr << "Failed to create directory "<< fold
            << " with errno " << errno << "\n";

    ofstream op0(fold+"/"+"clustersFile.txt");

    for (int i=0; i<NX*NY*NZ; ++i)
        op0<<C1[i]<<endl;

    op0.close();

    ofstream op01(fold+"/"+"clustersdata.txt");
        op01<<"           intrested phase :  "<<intrestPhase<<"     Number of the phases: "<<noOfPhases<<endl;
        op01<<     "-------------------------------------------------------------------------"<<endl;
        op01<<"connected cluster name"<<"        "<<"volume fraction custer/volume fraction intrested phase"<<endl;

	int c=0;
	for (auto & elem : clusters_vol_frac)
        {op01<<"        "<<c+phasePlus<<"                                   "<<clusters_vol_frac[c]<<endl; c++;}
}
//------------------------------------------------------------------------------

void percolation::vtiWriter(std::string fileName)
{	
	vtkSmartPointer<vtkImageData> imageData =
	vtkSmartPointer<vtkImageData>::New();
	imageData->SetDimensions(NX, NY, NZ);//
	#if VTK_MAJOR_VERSION <= 5
	imageData->SetNumberOfScalarComponents(1);
	imageData->SetScalarTypeToDouble();
	#else
	imageData->AllocateScalars(VTK_DOUBLE, 1);
	#endif
	int* dims = imageData->GetDimensions();

	// Fill every entry of the image data with "2.0"
	for (int z = 0; z < dims[2]; z++)
		for (int y = 0; y < dims[1]; y++)
			for (int x = 0; x < dims[0]; x++)
			{
				double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
				int len=x+y*dims[0]+z*dims[0]*dims[1];
				pixel[0] =  C1[len];
			}

	vtkSmartPointer<vtkXMLImageDataWriter> writer =
	vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(fileName.c_str());
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(imageData->GetProducerPort());
	#else
	writer->SetInputData(imageData);
	#endif
	writer->Write();

}

//TODO find cluster size and statistically informations



// Destructor
//-----------------------------------------------------------------------------
percolation::~percolation()
{
    delete [] C;
    delete [] C1;

}





