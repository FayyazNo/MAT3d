#include"FFT_thermal3D.hpp"

/*
This code calculate effective thermal conductivity of the hetergenous materials
The accelarated algorithm presented in "Computation of thermal properties via 3D 
homogenization of multiphase materials using FFT-based accelerated scheme, arXive, 2015"
is used in this code.
*/


//----------------------------------------------------------------------------------

//Constructor set up the FFT plans and matrices size
FFT_analysis3d_thermal::FFT_analysis3d_thermal (int Nx1, int Ny1, int Nz1)
{  
    Nx=Nx1;
    Ny=Ny1;
    Nz=Nz1;
    xyz=Nx*Ny*Nz;

    itNo=0;

	sum_field_old=0;
	sum_field=0;
	diff_sum_field=0;

	avr_eqVar=0;
	aver_eqVarOld=0;
	diff_avr_eqVar=0;

  	fftw_init_threads();
	fftw_plan_with_nthreads(3);
   
    t1  =   new double [3*xyz];
    t2  =   new fftw_complex [3*Nx*Ny*(Nz/2+1)];
	t3  =   new fftw_complex [3*Nx*Ny*(Nz/2+1)];
    
    GF= new double [6*Nx*Ny*(Nz/2+1)];
    RVE.resize(1, xyz);
    heatFlux.resize(3, xyz);
    tempGradient.resize(3, xyz);

	//field_c.resize(3, xyz);

    //tau.resize(3, xyz);
    E.resize(3,1);
    
  /*  for(int i=0; i<3*xyz; i++)
    {
        t2[i][1]=0;
    }*/
    
    int rank=3;

    dims= new fftw_iodim64[rank];
    dims_o= new fftw_iodim64[rank];

    if(dims==NULL || dims_o==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    dims[0].n=Nx;
    dims[0].is=Nz*Ny*3;
    dims[0].os=(Nz/2+1)*Ny*3;
    dims[1].n=Ny;
    dims[1].is=Nz*3;
    dims[1].os=(Nz/2+1)*3;
    dims[2].n=Nz;
    dims[2].is=3;
    dims[2].os=3;
    
    dims_o[0].n=Nx;
    dims_o[0].is=(Nz/2+1)*Ny*3;
    dims_o[0].os=Nz*Ny*3;
    dims_o[1].n=Ny;
    dims_o[1].is=(Nz/2+1)*3;
    dims_o[1].os=Nz*3;
    dims_o[2].n=Nz;
    dims_o[2].is=3;
    dims_o[2].os=3;
    
    int howmany_rank=1;
    
    howmany_dims=new fftw_iodim64 [howmany_rank];
    
    if(howmany_dims==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    howmany_dims[0].n=3;
    howmany_dims[0].is=1;
    howmany_dims[0].os=1;
    
    p3=fftw_plan_guru64_dft_r2c(rank, dims,  howmany_rank, howmany_dims,t1, t2, FFTW_ESTIMATE);
    p4=fftw_plan_guru64_dft_c2r(rank, dims_o,howmany_rank, howmany_dims,t3, t1, FFTW_ESTIMATE);

}



// Calculating Green operator and store it in 1D array
// Ir should be noted for each voxel only 6 independent 
// elemnt (out of nine) of the green operator is stored in GF array
void FFT_analysis3d_thermal::setGF()
{
    double mm;
    int x1, y1, z1;
    int cc;
    
    if (tempGradientBase)
    {
		 for (int x=0; x<Nx; x++)
		    for (int y=0; y<Ny; y++)
		        for (int z=0; z<(Nz/2+1); z++)
		        {
		            //cout<<x<<"    "<<y<<endl;
		            cc= (z+(Nz/2+1)*y+Ny*(Nz/2+1)*x)*6;
		            
		           if(x<=Nx/2) {x1=x*1.;}
		            else        {x1=-(Nx-x)*1.;}
		            								   
		            if(y<=Ny/2) {y1=y*1.;}
		            else        {y1=-(Ny-y)*1.;}
		            
		            if(z<=Nz/2) {z1=z*1.;}
		            else        {z1=-(Nz-z)*1.;}
		            
			
          
		            double zet [3]={z1*1.0, y1*1.0, x1*1.0};
		            
		            mm=zet[0]*(K0(0,0)*zet[0]+K0(0,1)*zet[1]+K0(0,2)*zet[2])+\
		               zet[1]*(K0(1,0)*zet[0]+K0(1,1)*zet[1]+K0(1,2)*zet[2])+\
		               zet[2]*(K0(2,0)*zet[0]+K0(2,1)*zet[1]+K0(2,2)*zet[2]);
		            
		            
		            GF[cc+0]=zet[0]*zet[0]/mm;  //GF11
		            GF[cc+1]=zet[0]*zet[1]/mm;  //GF12
		            GF[cc+2]=zet[0]*zet[2]/mm;  //GF13
		            GF[cc+3]=zet[1]*zet[1]/mm;  //GF22
		            GF[cc+4]=zet[1]*zet[2]/mm;  //GF23
		            GF[cc+5]=zet[2]*zet[2]/mm;  //GF33
		            
		        }

	}
    else
    {
		for (int x=0; x<Nx; x++)
		    for (int y=0; y<Ny; y++)
		        for (int z=0; z<(Nz/2+1); z++)
		        {
		            //cout<<x<<"    "<<y<<endl;
		            cc= (z+(Nz/2+1)*y+Ny*(Nz/2+1)*x)*6;
		            
		            if(x<=(Nx-1)/2) {x1=x*1.;}
		            else            {x1=-(Ny-x)*1.;}
		            								   
		            if(y<=(Ny-1)/2) {y1=y*1.;}
		            else            {y1=-(Ny-y)*1.;}
		            
		            if(z<=(Nz-1)/2) {z1=z*1.;}
		            else            {z1=-(Nz-z)*1.;}
		            
		            double zet [3]={z1*1.0, y1*1.0, x1*1.0};
		            
		            mm=(zet[0]*zet[0]+zet[1]*zet[1]+zet[2]*zet[2]);
		            
		            
		            GF[cc+0]=(1-zet[0]*zet[0]/mm)/R0(1,1);  //GF11
		            GF[cc+1]=( -zet[0]*zet[1]/mm)/R0(1,1);  //GF12
		            GF[cc+2]=( -zet[0]*zet[2]/mm)/R0(1,1);  //GF13
		            GF[cc+3]=(1-zet[1]*zet[1]/mm)/R0(1,1);  //GF22
		            GF[cc+4]=( -zet[1]*zet[2]/mm)/R0(1,1);  //GF23
		            GF[cc+5]=(1-zet[2]*zet[2]/mm)/R0(1,1);  //GF33
		           
		            
		        } 
	}

} 

//Multiply heat_flux(x)=-K(x)*tepmetature_gradient(x)
// In this case only 2 phases is assumed
// tempGradient and heatFlux dimension are (3,Nx*Ny*Nz)
void FFT_analysis3d_thermal::multipKx()
{

    if (tempGradientBase)
		for(int i=0; i<xyz; ++i)
		{
			if (RVE(i))
				heatFlux.col(i)=Km*tempGradient.col(i);
			else
				heatFlux.col(i)=Ki*tempGradient.col(i);
		}


	else

		for(int i=0; i<xyz; ++i)
		{
			if (RVE(i))
				tempGradient.col(i)=Rm*heatFlux.col(i);
			else
				tempGradient.col(i)=Ri*heatFlux.col(i);
		}
	
}


//Calculate FFT(tau) and store it in t2
// Guru interface is used to calculate FFT of the tau vector field
//It means that FFT of all 3 rows of the tau is done at same time in one transform
void FFT_analysis3d_thermal::FFT()
{   
	if (tempGradientBase)
    	memcpy(t1, heatFlux.data(), 3*xyz*sizeof t1);
	else
		memcpy(t1, tempGradient.data(), 3*xyz*sizeof t1);

    fftw_execute(p3);
}


//FFT backward of the GF*tau-----------------------------
void FFT_analysis3d_thermal::FFTBCKWRD()
{
    fftw_execute(p4);

	if (tempGradientBase)
	{
	    memcpy(heatFlux.data(), t1, 3*xyz*sizeof t1);
    	heatFlux=heatFlux/xyz;
	}

    else
	{
		memcpy(tempGradient.data(), t1,  3*xyz*sizeof t1);
		tempGradient=tempGradient/xyz;
	}

}


void FFT_analysis3d_thermal::strainUpdate()
{
	if (tempGradientBase)
	{
		K1=2*(Km+K0).inverse()*K0;
		K2=2*(Ki+K0).inverse()*K0;

		for(int i=0; i<xyz; ++i)
		{
			if (RVE(i))
				tempGradient.col(i)+=K1*(heatFlux.col(i)-tempGradient.col(i));
			else
				tempGradient.col(i)+=K2*(heatFlux.col(i)-tempGradient.col(i));
		}		
	}
	else
	{
		R1=2*(Rm+R0).inverse()*R0;
		R2=2*(Ri+R0).inverse()*R0;
 
		for(int i=0; i<xyz; ++i)
		{
			if (RVE(i))
				heatFlux.col(i)+=R1*(tempGradient.col(i)-heatFlux.col(i));
			else
				heatFlux.col(i)+=R2*(tempGradient.col(i)-heatFlux.col(i));
		}	
	}
}
//Multiply GF*FFT(tau)----------------------------
void FFT_analysis3d_thermal::multiplyGF()
{
    MatrixIdd  w(3,3);
    w<< 0,  1,   2,  
        1,  3,   4,  
        2,  4,   5; 

			for(int _k=0; _k<(Nz/2+1)*Nx*Nz; ++_k)
			{  
				int u=_k*3;
				int d=_k*6;
				for(int j=0; j<3; j++)
				{ 
				    t3[u+j][0]=(GF[d+w(j,0)]*t2[u+0][0]+GF[d+w(j,1)]*t2[u+1][0]+GF[d+w(j,2)]*t2[u+2][0]);
				    t3[u+j][1]=(GF[d+w(j,0)]*t2[u+0][1]+GF[d+w(j,1)]*t2[u+1][1]+GF[d+w(j,2)]*t2[u+2][1]);         
				}
			}
    
    t3[0][0]=E(0)*xyz;
    t3[0][1]=0;
    
    t3[1][0]=E(1)*xyz;
    t3[1][1]=0;
    
    t3[2][0]=E(2)*xyz;
    t3[2][1]=0;
    
}
 


//Checking Equilibrium in Fourier space and diffrence of the
// tempGradient in two successive iterations
void FFT_analysis3d_thermal::convergenceCheck()
{
    memcpy(t1, heatFlux.data(), 3*xyz*sizeof t1);
    fftw_execute(p3); 
    itNo+=1;
    double cc=0;
    double dd=sqrt((t2[0][0] *t2[0][0])+(t2[0][1] *t2[0][1])
                  +(t2[1][0] *t2[1][0])+(t2[1][1] *t2[1][1])
                  +(t2[2][0] *t2[2][0])+(t2[2][1] *t2[2][1]));

    double x1=0, y1=0, z1=0;
    
    for (int x=0; x<Nx; x++)
        for (int y=0; y<Ny; y++)
            for (int z=0; z<Nz/2+1; z++)
            {
                int i= (z+y*(Nz/2+1)+ x*Ny*(Nz/2+1))*3;
                
	            if(x<=Nx/2) {x1=x*1.;}
	            else        {x1=-(Nx-x)*1.;}
	            								   
	            if(y<=Ny/2) {y1=y*1.;}
	            else        {y1=-(Ny-y)*1.;}
	            
	            if(z<=Nz/2) {z1=z*1.;}
	            else        {z1=-(Nz-z)*1.;}
                
                double zet [3]={z1, y1, x1}; // frequency vector or each voxel
                
                cc+= sqrt ((t2[i+0][0] *t2[i+0][0]+t2[i+0][1] *t2[i+0][1])*zet[0]*zet[0]
                          +(t2[i+1][0] *t2[i+1][0]+t2[i+1][1] *t2[i+1][1])*zet[1]*zet[1]
                          +(t2[i+2][0] *t2[i+2][0]+t2[i+2][1] *t2[i+2][1])*zet[2]*zet[2]);
            }
 
    
    avr_eqVar=cc/dd/xyz;

    cout<<"Iter="<<itNo<<"   eqVar= "<<avr_eqVar<<"   ";


	if (tempGradientBase==true)
	{

        diff_sum_field=abs((heatFlux.sum()-sum_field_old)/(xyz*(0.5*Ki(0,0)+0.5*Km(0,0))*E.sum()));

		diff_avr_eqVar=abs(avr_eqVar-aver_eqVarOld);

		cout<<"diff_sum_field= "<<diff_sum_field<<"    "<<"diff_avr_eqVar= "<<diff_avr_eqVar<<"    "<<endl<<endl;
	}
    
	else if (heatFluxBase==true)
	{

        diff_sum_field=abs((tempGradient.sum()-sum_field_old)/(xyz*(0.5*Ri(0,0)+0.5*Rm(0,0))*E.sum()));

		diff_avr_eqVar=abs(avr_eqVar-aver_eqVarOld);

		cout<<"diff_sum_field= "<<diff_sum_field<<"    "<<"diff_avr_eqVar= "<<diff_avr_eqVar<<"    "<<endl<<endl;
	}
    
		aver_eqVarOld=avr_eqVar;

}




// Set up materil properties where ki and km are scaler thermal conductivity of the inclusion and matrix
// respectively and Ki=ki*I and Km=km*I where Ki, Km and I are thermal conductivity of the inclusion
// matrix and identity tensors respectively.
void FFT_analysis3d_thermal::material (double km, double ki, MatrixIdd &RVE1)
{ 
    heatFluxBase=false;
    tempGradientBase=false;


    Km.resize(3,3);
    Ki.resize(3,3);
    K0.resize(3,3);
    Rm.resize(3,3);
    Ri.resize(3,3);
    R0.resize(3,3);

    Km<<   km,   0,    0,    
           0,    km,   0,    
           0 ,   0,    km;    
    
    
    Ki<<   ki,   0,    0,    
           0,    ki,   0,    
           0 ,   0,    ki;    

    Rm=Km.inverse();
    Ri=Ki.inverse();

	RVE=RVE1;
    double vol_1= RVE.sum()/xyz;

   if ((vol_1>=0.5 and km>ki) or (vol_1<0.5 and km<ki) )
	{
		tempGradientBase=true;
        std::cout<< " temperature Gradient-Based Algorithm is Running"<<std::endl;
	}
   else
	{
		heatFluxBase=true; 
		std::cout<< " Heat Flux-Based  Algorithm is Running"<<std::endl;
	}
    


	double k0; 
	k0=sqrt(ki*km);

	K0<<   k0,   0,    0,    
		   0,    k0,   0,    
		   0 ,   0,    k0; 
	R0=K0.inverse();
    setGF();
}


//Impose uniform temperature gradient E(1,3) to all voxels
void FFT_analysis3d_thermal::setLoad(MatrixDdd &E1)
{ 
    E=E1;

    if (tempGradientBase)
	{
		for (int i=0; i<xyz; ++i)
        	tempGradient.col(i)=E1; 
	}
	else 
	{	
		for (int i=0; i<xyz; ++i)
        	heatFlux.col(i)=E1; 
	}

}


// Run a single analysis and try to find equilibrium state after impose uniform teperature gradient to all voxels.
// The argument is number of iterations. It give also the time of the each section in each iteration.
void FFT_analysis3d_thermal::runAnalysis(double tol_field, double tol_eq)
{ 
    multipKx(); 
	while (true)
	{
 		if (tempGradientBase)
		{
			sum_field_old=heatFlux.sum(); 
			heatFlux=K0*tempGradient-heatFlux;
		}
		else
		{
			sum_field_old=tempGradient.sum();
			tempGradient=R0*heatFlux-tempGradient;
		}

		FFT();

		multiplyGF();

		FFTBCKWRD();

		strainUpdate();

		multipKx(); 

		convergenceCheck();

		if (itNo>1 && diff_avr_eqVar<tol_eq  &&  diff_sum_field<tol_field)
				break;
		if (itNo>500)
		{
			std::cout<<" iteration excede from allowable number"<<std::endl;
			break;
		}

	}        


}

// Perform 3 analysis with different initial temperature gradient and gives homogenized thermal conductivity tensor of the RVE
// The argument is number of iteration in each analysis
void FFT_analysis3d_thermal::homogenization(double tol_field, double tol_eq, string filename)
{  
    ofstream homo(filename+".txt");
    tempGradHomo.resize(3,3);
    KHomo.resize(3,3);
    RHomo.resize(3,3);
    heatFluxHomo.resize(3,3);

	if (tempGradientBase)
	{ 

{
		for (int w=0; w<3; ++w)
		{
            itNo=0;
            sum_field_old=0;
		    E.setZero();
		    E(w)=1;
            setLoad(E); 
            cout<<E<<endl;
		    runAnalysis(tol_field, tol_eq);
		    
		    double e1=0, e2=0, e3=0;
		    double s1=0, s2=0, s3=0;
		    
		    for (int t=0; t<xyz; ++t)
		    {
		        e1+=tempGradient(0,t);
		        e2+=tempGradient(1,t);
		        e3+=tempGradient(2,t);
		        
		        s1+=heatFlux(0,t);
		        s2+=heatFlux(1,t);
		        s3+=heatFlux(2,t);
		    }
		    
		    tempGradHomo(w,0)=e1/xyz;
		    tempGradHomo(w,1)=e2/xyz;
		    tempGradHomo(w,2)=e3/xyz;
		    
		    heatFluxHomo(0,w)=s1/xyz;
		    heatFluxHomo(1,w)=s2/xyz;
		    heatFluxHomo(2,w)=s3/xyz;
		}
}		
		KHomo.col(0)=tempGradHomo.inverse()* heatFluxHomo.col(0);
		KHomo.col(1)=tempGradHomo.inverse()* heatFluxHomo.col(1);
    	KHomo.col(2)=tempGradHomo.inverse()* heatFluxHomo.col(2);

		RHomo=KHomo.inverse();

	}

	else    
	{ 
		for (int w=0; w<3; ++w)
		{    
            itNo=0;
		    E.setZero();
		    E(w)=1;
		    setLoad(E);
			cout<<E<<endl;
		    runAnalysis(tol_field, tol_eq);
		    
		    double e1=0, e2=0, e3=0;
		    double s1=0, s2=0, s3=0;
		    
		    for (int t=0; t<xyz; ++t)
		    {
		        e1+=heatFlux(0,t);
		        e2+=heatFlux(1,t);
		        e3+=heatFlux(2,t);
		        
		        s1+=tempGradient(0,t);
		        s2+=tempGradient(1,t);
		        s3+=tempGradient(2,t);
		    }
		    
		    heatFluxHomo(w,0)=e1/xyz;
		    heatFluxHomo(w,1)=e2/xyz;
		    heatFluxHomo(w,2)=e3/xyz;
		    
		    tempGradHomo(0,w)=s1/xyz;
		    tempGradHomo(1,w)=s2/xyz;
		    tempGradHomo(2,w)=s3/xyz;
		}
		
		RHomo.col(0)=heatFluxHomo.inverse()* tempGradHomo.col(0);
		RHomo.col(1)=heatFluxHomo.inverse()* tempGradHomo.col(1);
		RHomo.col(2)=heatFluxHomo.inverse()* tempGradHomo.col(2);
		
		KHomo=RHomo.inverse();
	}		
    //Mori-Tanaka effective thermal conductivity
	double k_homo_mori_tanaka;
	double km=Km(0,0);
	double ki=Ki(0,0);
	double fi_m=RVE.sum()/xyz;     // volume fraction of the white phase
	double fi_i=1-fi_m;            // volume fraction of the black phase
	k_homo_mori_tanaka=km+(fi_i* (ki-km))/(1+fi_m*(ki-km)/(3*km));
    //Lower and upper Hashin-Shtrikman bounds---------

	double k_up_hashin_shtrikman=  fi_i* ki+ fi_m*km-fi_i*fi_m*(km-ki)*(km-ki)/(3*ki-fi_i*(ki-km));
	double k_low_hashin_shtrikman= fi_i* ki+ fi_m*km-fi_i*fi_m*(km-ki)*(km-ki)/(3*km+fi_m*(ki-km));
//------------------
	K11=KHomo(0,0); K22=KHomo(1,1); K33=KHomo(2,2); 


    homo<<"Volume Fraction of the White Phase (or 1)="<<fi_m<<endl<<endl;
    homo<<"K phase matrix  ( or 1)="<<endl<<Km<<"     "<<endl;
    homo<<"K phase inclusion(or 0)="<<endl<<Ki<<"     "<<endl;
    homo<<"Nx="<<Nx<<"   "<<"Ny="<<Ny<<"   "<<Nz<<endl<<endl;
    homo<<endl<<"KHomogenized Mori-Tanaka="<<endl<<endl;
    homo<<k_homo_mori_tanaka<<endl;
    homo<<endl<<"Lower bound HS="<<endl<<endl;
    homo<<k_low_hashin_shtrikman<<endl;
    homo<<endl<<"Upper bound LS="<<endl<<endl;
    homo<<k_up_hashin_shtrikman<<endl;
    homo<<endl<<"KHomogenized="<<endl<<endl;
    homo<<KHomo<<endl;

    homo<<endl<<"RHomogenized(Homogenized Resistivity tensor)="<<endl<<endl;
    homo<<RHomo<<endl;
    homo<<endl<<"Rm="<<endl<<endl;
    homo<<Rm<<endl;
    homo<<endl<<"Ri="<<endl<<endl;
    homo<<Km<<endl;
    homo<<endl<<endl<<"--------------------------------------------------------------------"<<endl<<endl;
    homo<<endl<<"Homogenized Load Cases (Each row is one load case i.e (E11 E22 E33 ) and Eij is summation of temperatur gradient components in in  all voxels )"<<endl<<endl;
    homo<<tempGradHomo<<endl;
    homo<<endl<<"Homogenized stess (Each coulemn is (SS11 SS22 SS33) and SSij is summation of heat flux in in all voxels )"<<endl<<endl;
    homo<<heatFluxHomo<<endl;

}


// Gives all 3 component of heat flux and all 3 component of temperatuer gradient
// Argument is name of the folder which files store.
// Each file is single coloumn .txt file with Nx*Ny*Nz size
void FFT_analysis3d_thermal::printer (string folderName)
{

    int ss= mkdir(folderName.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    
    int x1=0;
    int x2=Nx*Ny*Nz;
    ofstream op0("heatFlux1.txt");
    for (int i=x1; i<x2; ++i)
        op0<<heatFlux(0,i)<<endl;
    
    ofstream op1("heatFlux2.txt");
    for (int i=x1; i<x2; ++i)
        op1<<heatFlux(1,i)<<endl;
    
    ofstream op2("heatFlux3.txt");
    for (int i=x1; i<x2; ++i)
        op2<<heatFlux(2,i)<<endl;
    
    
    ofstream op00("tempGradient1.txt");
    for (int i=x1; i<x2; ++i)
        op00<<tempGradient(0,i)<<endl;
    
    ofstream op01("tempGradient2.txt");
    for (int i=x1; i<x2; ++i)
        op01<<tempGradient(1,i)<<endl;
    
    ofstream op02("tempGradient3.txt");
    for (int i=x1; i<x2; ++i)
        op02<<tempGradient(2,i)<<endl;
    
    
    
    string s1= "/";
    string s2= ".txt";
    for(int i=0; i<3; ++i)
    {
        string s3= to_string(i+1);
        string s4="heatFlux";
        s4+=s3+s2;
        ifstream source(s4.c_str()  , ios::binary);
        string s5=folderName+s1+s4;
        ofstream dest(s5.c_str(), ios::binary);
        dest << source.rdbuf();
        remove (s4.c_str());
        s4.clear();
        s5.clear();
    }
    
    for(int i=0; i<3; ++i)
    {
        string s3= to_string(i+1);
        string s4="tempGradient";
        s4+=s3+s2;
        ifstream source(s4.c_str()  , ios::binary);
        string s5=folderName+s1+s4;
        ofstream dest(s5.c_str(), ios::binary);
        dest << source.rdbuf();
        remove (s4.c_str());
        s4.clear();
        s5.clear();
    }

}
//---

//vtiwriter---------------------------------------------------------------------------

void FFT_analysis3d_thermal::vtiWriter(std::string fileName)
{
	
	  vtkSmartPointer<vtkImageData> imageData =
	  vtkSmartPointer<vtkImageData>::New();
	  imageData->SetDimensions(Nx, Ny, Nz);//
#if VTK_MAJOR_VERSION <= 5
	  imageData->SetNumberOfScalarComponents(6);
	  imageData->SetScalarTypeToDouble();
#else
	  imageData->AllocateScalars(VTK_DOUBLE, 6);
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
            int len=x+y*dims[0]+z*dims[0]*dims[1];
		    pixel[0] =  tempGradient(0, len);
		    pixel[1] =  tempGradient(1, len);
		    pixel[2] =  tempGradient(2, len);

		    
            pixel[3] =    heatFlux(0, len);
		    pixel[4] =    heatFlux(1, len);
		    pixel[5] =    heatFlux(2, len);



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

// Destructure----------------------------------------
FFT_analysis3d_thermal::~FFT_analysis3d_thermal()
{ 
    if(t1) delete [] t1;
    if(t2) delete [] t2;
    if(GF) delete [] GF;
	if (t3)delete [] t3;

    delete [] dims;
	delete [] dims_o;
    delete [] howmany_dims;
    
    fftw_destroy_plan(p3);
    fftw_destroy_plan(p4);
    fftw_cleanup();
	fftw_cleanup_threads();

    RVE.resize(0,0);   
    heatFlux.resize(0,0);
    tempGradient.resize(0,0);
    tau.resize(0,0);

}


