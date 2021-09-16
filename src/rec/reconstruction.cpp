#include "reconstruction.hpp"
using namespace recon;


//fftshift****************************************************************************
/* This function send zero ferequency to the center in a 3d matrix

examples of 2d fftshift are shown, 2D is choosed just to 
portray logic of fftshift and the function fftshift works just for 3D matrix.

1-Where Nx and Ny are odd

a=[  1   2   3  |4  5
     6   7   8  |9  10
     11  12  13 |14 15
    ------------------
     16  17  18 |19 20 
     21  22  23 |24 25]

fftshift (a)=

  [19  20| 18  17  18
   24  25| 21  22  23
   -------------------
   4   5 | 1   2   3
   9   10| 6   7   8
   14  15| 11  12  13]

2-Where Nx and Ny are even
  a=[1   2  | 3   4  
     5   6  | 7   8
    ------------------
     9   10 | 11  12  
    13   14 | 15  16]

fftshift (a)=

  a=[11   12 | 9   10 
     15   16 | 13  14
    ------------------
     3    4  | 1   2  
     7    8  | 5   6]

*/

//******************************************************************************
void recon::fftshift(double *in, double *out, int Nx, int Ny, int Nz)
{    


	int i_mid= Nx/2,   j_mid= Ny/2,   k_mid= Nz/2;
	cout<<k_mid<<"   "<<j_mid<<"  "<<i_mid<<endl;

	//swaping
	int i_shift, j_shift, k_shift;
	for (int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{

				if      (~Nx%2 && i<i_mid)  {i_shift= Nx/2;}

				else if (~Nx%2 && i>=i_mid) {i_shift= -Nx/2;}

				else if ( Nx%2 && i<i_mid)  {i_shift= (Nx/2);}

				else if ( Nx%2 && i>=i_mid) {i_shift= -(Nx/2)-1;}
				//-------------------------------------------     
				if      (~Ny%2 && j<j_mid)  {j_shift= Ny/2;}

				else if (~Ny%2 && j>=j_mid) {j_shift= -Ny/2;}

				else if ( Ny%2 && j<j_mid)  {j_shift= (Ny/2);}

				else if ( Ny%2 && j>=j_mid) {j_shift= -(Ny/2)-1;}
				//---------------------------------------------
				if      (~Nz%2 && k<k_mid)  {k_shift= Nz/2;}

				else if (~Nz%2 && k>=k_mid) {k_shift= -Nz/2;}

				else if ( Nz%2 && k<k_mid)  {k_shift= (Nz/2);}

				else if ( Nz%2 && k>=k_mid) {k_shift= -(Nz/2)-1;}
				//---------------------------------------------

				out[(k+k_shift)+(Nz)*(j+j_shift)+(Nz*Ny)*(i+i_shift)]=in[k+(Nz)*j+(Nz*Ny)*i];

			}
}

//ifftshift**********************************************************************
/*
Just look at the fftshift comments and pay attention:

a=ifftshift(fftshift(a))

*/

void recon::ifftshift(double *in, double *out, int Nx, int Ny, int Nz)
{  

	int i_mid= Nx/2,   j_mid= Ny/2,   k_mid= Nz/2;
	cout<<k_mid<<"   "<<j_mid<<"  "<<i_mid<<endl;

	//swaping
    int i_shift, j_shift, k_shift;
	for (int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{

				if      (~Nx%2 && i<i_mid)  {i_shift= Nx/2;}

				else if (~Nx%2 && i>=i_mid) {i_shift= -Nx/2;}

				else if ( Nx%2 && i<i_mid)  {i_shift= (Nx/2)+1;}

				else if ( Nx%2 && i>=i_mid) {i_shift= -(Nx/2);}
				//-------------------------------------------     
				if      (~Ny%2 && j<j_mid)  {j_shift= Ny/2;}

				else if (~Ny%2 && j>=j_mid) {j_shift= -Ny/2;}

				else if ( Ny%2 && j<j_mid)  {j_shift= (Ny/2)+1;}

				else if ( Ny%2 && j>=j_mid) {j_shift= -(Ny/2);}
				//---------------------------------------------
				if      (~Nz%2 && k<k_mid)  {k_shift= Nz/2;}

				else if (~Nz%2 && k>=k_mid) {k_shift= -Nz/2;}

				else if ( Nz%2 && k<k_mid)  {k_shift= (Nz/2)+1;}

				else if ( Nz%2 && k>=k_mid) {k_shift= -(Nz/2);}
				//---------------------------------------------

				out[(k+k_shift)+(Nz)*(j+j_shift)+(Nz*Ny)*(i+i_shift)]=in[k+(Nz)*j+(Nz*Ny)*i];

			}
}


//reconstruction class
//*****************************************************************************
recPlan::recPlan(CImg<unsigned char> c1, CImg<unsigned char> c2)
{ 
	cd=1;
	CImg<unsigned char> C1=c1;
	CImg<unsigned char> C2=c2;
	Nx = C1.width();
	Ny = C1.height();
	Nz = C2.height();
	Nzz= C2.width();

	fftw_init_threads();
	fftw_plan_with_nthreads(3);

	xyz=Nx*Ny*Nz;
	cout<<"Nx=  "<<Nx<<"  Ny= "<<Ny<<"  Nz= "<<Nz<<"  xyz= "<<xyz<<endl;
	normC  = new double [Nx*Ny*Nz];
	mHat   = new fftw_complex [Nx*Ny*Nz];
	MHat   = new fftw_complex [Nx*Ny*Nz];
	mTilda = new fftw_complex [Nx*Ny*Nz];
	MTilda = new fftw_complex [Nx*Ny*Nz];
	CC1    = new fftw_complex [Nx*Ny];
	CC3    = new fftw_complex [Nx*Ny];
	CC2    = new fftw_complex [Nz*Nzz];
	CC4    = new fftw_complex [Nz*Nzz];

	MatrixDdd RVE (1,Nx*Ny*Nz);
	temp.resize(1,Nx*Ny*Nz);
	ang.resize (1,Nx*Ny*Nz);
	ang1.resize(1,Nx*Ny*Nz);
	ang2.resize(1,Nx*Ny*Nz);

	p1=fftw_plan_dft_3d(Nx, Ny, Nz, mTilda, MTilda, FFTW_FORWARD,  FFTW_ESTIMATE); //fft3 plan 
	p2=fftw_plan_dft_3d(Nx, Ny, Nz, MTilda, mTilda, FFTW_BACKWARD, FFTW_ESTIMATE); // ifft3 plan
	p3=fftw_plan_dft_3d(Nx, Ny, Nz, mHat, MHat, FFTW_FORWARD,  FFTW_ESTIMATE); //fft3 plan 
	p4=fftw_plan_dft_3d(Nx, Ny, Nz, MHat, mHat, FFTW_BACKWARD, FFTW_ESTIMATE); // ifft3 plan

	p5=fftw_plan_dft_2d(Nx, Ny, CC1, CC3, FFTW_FORWARD,  FFTW_ESTIMATE); //fft3 plan 
	p6=fftw_plan_dft_2d(Nx, Ny, CC1, CC3, FFTW_BACKWARD, FFTW_ESTIMATE); // ifft3 plan

	p7=fftw_plan_dft_2d(Nz, Nzz, CC2, CC4, FFTW_FORWARD,  FFTW_ESTIMATE); //fft3 plan 
	p8=fftw_plan_dft_2d(Nz, Nzz, CC2, CC4, FFTW_BACKWARD, FFTW_ESTIMATE); // ifft3 plan

	int Max=0;
	for (int i=0; i<Nx; ++i)
	  for(int j=0; j<Ny; ++j)
		  Max=max(Max, (int)(C1.atXY(i,j)));
	for (int j=0; j<Nx; ++j)
	  for (int k=0; k<Ny; ++k)
	  {
		  CC1[k+j*Ny][0]=(int)(C1.atXY(j,k))/Max;
		  CC1[k+j*Ny][1]=0;
	  }
	Max=0;
	for (int i=0; i<Nx; ++i)
	  for(int j=0; j<Ny; ++j)
		  Max=max(Max, (int)(C2.atXY(i,j)));
	for (int j=0; j<Nx; ++j)
	  for (int k=0; k<Ny; ++k)
	  {
		  CC2[k+j*Ny][0]=(int)(C2.atXY(j,k))/Max;
		  CC2[k+j*Ny][1]=0;
	  }
	}


//*****************************************************************************
recPlan::recPlan(double *TPCF, int Nx1, int Ny1, int Nz1)
{ 
	cd=2;
	Nx=Nx1;
	Ny=Ny1;
	Nz=Nz1;
	xyz=Nx*Ny*Nz;

	fftw_init_threads();
	fftw_plan_with_nthreads(3);

	cout<<"Nx=  "<<Nx<<"  Ny= "<<Ny<<"  Nz= "<<Nz<<"  xyz= "<<xyz<<endl;
	normC  = new double [Nx*Ny*Nz];
	mHat   = new fftw_complex [Nx*Ny*Nz];
	MHat   = new fftw_complex [Nx*Ny*Nz];
	mTilda = new fftw_complex [Nx*Ny*Nz];
	MTilda = new fftw_complex [Nx*Ny*Nz];

	MatrixDdd RVE (1,Nx*Ny*Nz);
	temp.resize(1,Nx*Ny*Nz);
	ang.resize (1,Nx*Ny*Nz);
	ang1.resize(1,Nx*Ny*Nz);
	ang2.resize(1,Nx*Ny*Nz);

	p1 = fftw_plan_dft_3d(Nx, Ny, Nz, mTilda, MTilda, FFTW_FORWARD,  FFTW_ESTIMATE); //fft3 plan 
	p2 = fftw_plan_dft_3d(Nx, Ny, Nz, MTilda, mTilda, FFTW_BACKWARD, FFTW_ESTIMATE); // ifft3 plan
	p3 = fftw_plan_dft_3d(Nx, Ny, Nz, mHat, MHat, FFTW_FORWARD,  FFTW_ESTIMATE); //fft3 plan 
	p4 = fftw_plan_dft_3d(Nx, Ny, Nz, MHat, mHat, FFTW_BACKWARD, FFTW_ESTIMATE); // ifft3 plan

	for (int i=0; i<Nx*Ny*Nz; ++i)
	   normC[i]=sqrt(abs(TPCF[i])*xyz);
}

//Estimation of the TPCF, it needed when first constructor is called*************
void recPlan::estimateNormFFT_RVE()
{
    fftw_execute (p5);
    for (int i=0; i<Nx*Ny; ++i)
    {
        CC1[i][0]=(CC3[i][0]*CC3[i][0]+CC3[i][1]*CC3[i][1])/(Nx*Ny*Nx*Ny);
        CC1[i][1]=0;
    }

    fftw_execute (p6); //put Norm FFT CC1 (i.e. ifft(fft(CC1)*conjugate(fft(CC1))/Nx/Ny))----> CC3

    fftw_execute (p7);
    for (int i=0; i<Nz*Nzz; ++i)
    {
        CC2[i][0]=(CC4[i][0]*CC4[i][0]+CC4[i][1]*CC4[i][1])/(Nz*Nzz*Nz*Nzz);
        CC2[i][1]=0;
    }
    fftw_execute (p8); //put Norm FFT CC2 (i.e. ifft(fft(CC2)*conjugate(fft(CC2))/Nx/Ny))----> CC4
    

//TPCF Estimation***************************************************************
    double f1=CC3[0][0];
    cout<<endl<<f1;
    for (int i=0; i<Nx; ++i)
        for(int j=0; j<Ny; ++j) 
            for(int k=0; k<Nz; ++k)
            {
                mHat[k+j*Nz+i*Nz*Ny][0]=(CC3[i*Ny+j][0]*CC4[k][0]/f1+(f1-CC3[i*Ny+j][0])*(f1-CC4[k][0])/(1-f1))*xyz;
                mHat[k+j*Nz+i*Nz*Ny][1]=0;
            }

    fftw_execute (p3);
    for (int i=0; i<Nx*Ny*Nz; ++i)
        normC[i]=sqrt(sqrt((MHat[i][0]*MHat[i][0]+MHat[i][1]*MHat[i][1])));
    ofstream op4("normc.txt");
	for(int i=0; i<Nx*Ny*Nz; i++)
	    op4<<normC[i]<<endl;
}

//Running reconstruction procedure**********************************************
void recPlan::runRec(int n, double b)
{      
    // cout<<"Reconstruction is running. Please wait.";
    numloops=n;
    beta=b;
    trackrealerrorA.resize(1,n);
	ang=(ang.setRandom(1,xyz)+MatrixDdd::Constant(1, xyz, 1.))*0.5;

	for (int s=0; s<xyz; ++s) { mTilda[s][0]=ang(s);  mTilda[s][1]=0; }       
	int j=0;
    int i=0;
	 
    while (true)
	{ 
        j+=1; 
	    fftw_execute(p1);
            
        if (i>=n)
        {
           cout<<"The Iteration Number Exceeded"<<endl;
           break; 
        }

	    trackrealerrorA(i)=funcNormC (MTilda, normC, xyz) ;  // Norm (MTilda)-->mHat
             
        if (trackrealerrorA(i)<0.1)
        { 
           cout<< "Phase Recovery Converged"<<endl;
           break;
        } 
             
	    //cout<<"Itteration="<<i<<"     Error= "<<trackrealerrorA(i)<<endl;
	    memcpy(temp.data(), *MTilda, 2*xyz*sizeof MTilda);
	    ang=temp.array().arg();

	    if (remainder(i+1, 150)==0)
	    {
	       ang1=(ang.setRandom(1,xyz)+MatrixDdd::Constant(1, xyz, 1.))*0.5;
	              
	       for (int s=0; s<xyz; ++s){ mHat[s][0]=ang(s);  mHat[s][1]=0; }
	       
	       fftw_execute(p3);   
	       
	       memcpy(temp.data(), *MHat, 2*xyz*sizeof MHat);
	       ang2=temp.array().arg()*.3;
	       		
	       memcpy(temp.data(), *MTilda, 2*xyz*sizeof MHat);
	       ang1=temp.array().arg()*.7;
	       ang=ang2+ang1;
	    }
   
	    for (int s=0; s<xyz; ++s)
        { 
            MHat[s][0]=normC[s]*cos(ang(s));  
            MHat[s][1]=normC[s]*sin(ang(s)); 
        } 
              
	    fftw_execute(p4);     
	    memcpy(temp.data(), *mHat, 2*xyz*sizeof MTilda);
	    temp=temp/xyz;
	    memcpy(*mHat, temp.data(), 2*xyz*sizeof MTilda);
	   
	    if (j<100)
	    {
	       for (int s=0; s<xyz; ++s)
		   {
               if (mHat[s][0]<0){mHat[s][0]=mTilda[s][0]-beta*(mHat[s][0]);     mHat[s][1]=0;}
		       if (mHat[s][0]>1){mHat[s][0]=mTilda[s][0]-beta*(mHat[s][0]-1);   mHat[s][1]=0;} 
		       if (mHat[s][0]>1){mHat[s][0]=1;                                  mHat[s][1]=0;} 
		   }
	       memcpy(*mTilda, *mHat, 2*xyz*sizeof MTilda);
	    } 

	    if(j>=100)
	    {
	       for (int s=0; s<xyz; ++s)
		   {
		      if (mHat[s][0]<0){mHat[s][0]=0;     mHat[s][1]=0;}
		      if (mHat[s][0]>1){mHat[s][0]=1;     mHat[s][1]=0;}  
		   }
	       if(j>200) j=0;
	       memcpy(*mTilda, *mHat, 2*xyz*sizeof MTilda);
	    } 
	    i+=1;
	  }




	int ctr=0; double x=0, vol_dif=0;
	while (true)
	{		
		for (int s=0; s<xyz; ++s)  
			ang(s)=mTilda[s][0];

		for (int s=0; s<xyz; ++s)
		{	
			if (ang(s)< (.5-x))
				ang(s)=0; 
			else 
				ang(s)=1;
		}

		vol_dif=(ang.sum()/(1.*xyz)-normC[0]/(1.*xyz));

cout<<vol_dif<<endl;

		if (vol_dif>0) 
			x-=.001; 
		else 
			x+=.001;       

		if (abs(vol_dif)<.0001 or ctr>200) 
			break;

		ctr+=1;
	}	 


		for (int s=0; s<xyz; ++s) mTilda[s][0]= ang(s);

}

//Save all cut sections of the reconstructed RVE in a folder********************
void recPlan::saveRecImages(string locFldr)
{   
    mkdir(locFldr.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    for (int i=0; i<Nx; ++i)
    {
        std::string tFName=locFldr+"/";
        tFName= tFName+to_string(i+1)+".png";
        CImg<unsigned char> C(Nx, Ny, 1);
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
                C(j, k, 0)=mTilda[k+j*Nz+i*Ny*Nz][0]*255;
        C.save((tFName).c_str());
    }
	ofstream op2(locFldr+"/"+"rve_rec.txt");
	double s=0;
	for(int i=0; i<xyz; i++)
	{
	    op2<<mTilda[i][0]<<endl;
		s+=	mTilda[i][0];

	}
	s=s/xyz;

	ofstream op1(locFldr+"/"+"err.txt");
	op1<<"voloum fraction= "<<s<<std::endl;

	for(int i=0; i<numloops-1; i++)
    	op1<<trackrealerrorA(i)<<endl;
}

// vti format file writer***************************************************
void recPlan::vtiWriter(string fileName, double *tpcf)
{
	vtkSmartPointer<vtkImageData> imageData =
	vtkSmartPointer<vtkImageData>::New();
	imageData->SetDimensions(Nx, Ny, Nz);//
#if VTK_MAJOR_VERSION <= 5
	imageData->SetNumberOfScalarComponents(2);
	imageData->SetScalarTypeToDouble();
#else
	imageData->AllocateScalars(VTK_DOUBLE, 2);
#endif
	int* dims = imageData->GetDimensions();

	// Fill every entry of the image data with "2.0"
	for (int z = 0; z < dims[2]; z++)
	{
		for (int y = 0; y < dims[1]; y++)
		{
			for (int x = 0; x < dims[0]; x++)
			{
				double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
				pixel[0] = mTilda[x+y*dims[0]+z*dims[0]*dims[1]][0];
                pixel[1] =   tpcf[x+y*dims[0]+z*dims[0]*dims[1]];

			}
		}
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


recPlan::~recPlan()
{
    if (normC) delete [] normC;
    if (mTilda) delete [] mTilda;
    if (MTilda) delete [] MTilda;
    if (mHat) delete [] mHat;
    if (MHat) delete [] MHat;
	fftw_cleanup_threads();
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_destroy_plan(p4);
    
    if (cd==1)
    { 
        if (CC1) delete [] CC1;
        if (CC2) delete [] CC2;
        if (CC3) delete [] CC3;
        if (CC4) delete [] CC4;
        fftw_destroy_plan(p5);
        fftw_destroy_plan(p6);
        fftw_destroy_plan(p7);
        fftw_destroy_plan(p8);
    }

    temp.resize(0,0);
    ang.resize (0,0);
    ang1.resize(0,0);
    ang2.resize(0,0);
}

