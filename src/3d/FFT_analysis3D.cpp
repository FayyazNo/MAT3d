#include "FFT_analysis3D.hpp"

double FFT_Analysis3D::nrm2 (double *zet){return (zet[0]*zet[0]+zet[1]*zet[1]+zet[2]*zet[2]);}
double FFT_Analysis3D::dlt  (int i, int j){if (i==j) return 1.; else return 0.;} 
void FFT_Analysis3D::max_min_Array(const double  *array, int length, double out[4])
{
// out[0]=max(array), out[1]=index(max(array)), out[2]=min(array), out[3]=index(min[array]) 
	out[0]=array[0], out[1]=0; 	out[2]=array[0], out[3]=0;
	for (auto i=0; i<length-1; ++i)
	{
		if (array[i+1]>out[0])
		{
			out[0]=array[i+1];
			out[1]=i+1;
		}

		if (array[i+1]<out[2])
		{
			out[2]=array[i+1];
			out[3]=i+1;
		}
	}
}
//------------------------------------------------------------------------------
FFT_Analysis3D::FFT_Analysis3D(int Nx1, int Ny1, int Nz1, string mater_behav1)
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

    GF= new double [21*Nx*Ny*(Nz/2+1)];
    RVE.resize(1, xyz);

    epsilon.resize(6, xyz);
    sigma.resize(6, xyz);

    E.resize(6,1);
	X0.resize(6,1);

	I.resize(6,1);
	I<<1,1,1,0,0,0;

	mater_behav=mater_behav1;
//---------------------------------------------------------------
	if (mater_behav=="ep")
	{
		epsilon_old.resize(6, xyz);
	 	epsilon_old.setZero(6, xyz);

		sigma_old.resize(6, xyz);
		sigma_old.setZero(6, xyz);

		sigma_iter_old.resize(6,1);
		epsilon_iter_old.resize(6,1);

		epsilon_p_eff.resize(1, xyz);
		epsilon_p_eff.setOnes();
		epsilon_p_eff=epsilon_p_eff*.000000001236547;

		sigma_eff.resize(1, xyz);

		sigma_d.resize(6, 1);

		phi.resize(1, xyz);
		d_ep_p.resize(1, xyz);
		d_ep_p.setZero();
	// epsilon_p.resize(6, xyz);
	// epsilon_p.setZero();
	}

	if (mater_behav=="ed")
	{
cout<<"llllllllllllllllllll"<<endl;

		epsilon_old.resize(6, xyz);
	 	epsilon_old.setZero(6, xyz);

		d.resize(1,xyz);
		d.setZero(1, xyz);

		d_old.resize(1,xyz);
		d_old.setZero(1, xyz);

		// Starin energy
		Y.resize(1,xyz);
		Y.setZero(1, xyz);
	}

// fftw initialization------------------------------------------
 	fftw_init_threads();
	fftw_plan_with_nthreads(3);

    t1  =   new double [6*xyz];
    t2  =   new fftw_complex [6*Nx*Ny*(Nz/2+1)];
    t3  =   new fftw_complex [6*Nx*Ny*(Nz/2+1)];
    t4  =   new fftw_complex [6*Nx*Ny*(Nz/2+1)];


    int rank=3;
    dims= new fftw_iodim64[rank];
	dims_o= new fftw_iodim64[rank];
    
    if(dims==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    dims[0].n=Nx;
    dims[0].is=Nz*Ny*6;
    dims[0].os=(Nz/2+1)*Ny*6;
    dims[1].n=Ny;
    dims[1].is=Nz*6;
    dims[1].os=(Nz/2+1)*6;
    dims[2].n=Nz;
    dims[2].is=6;
    dims[2].os=6;


    dims_o[0].n=Nx;
    dims_o[0].is=(Nz/2+1)*Ny*6;
    dims_o[0].os=Nz*Ny*6;
    dims_o[1].n=Ny;
    dims_o[1].is=(Nz/2+1)*6;
    dims_o[1].os=Nz*6;
    dims_o[2].n=Nz;
    dims_o[2].is=6;
    dims_o[2].os=6;

    
    int howmany_rank=1;

    howmany_dims=new fftw_iodim64 [howmany_rank];
    
    if(howmany_dims==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    howmany_dims[0].n=6;
    howmany_dims[0].is=1;
    howmany_dims[0].os=1;
    

    p_t1_t2_frw=fftw_plan_guru64_dft_r2c(rank, dims,  howmany_rank, howmany_dims,  t1, t2, FFTW_ESTIMATE);
    p_t1_t4_frw=fftw_plan_guru64_dft_r2c(rank, dims,  howmany_rank, howmany_dims,  t1, t4, FFTW_ESTIMATE);

    p_t3_t1_bkw=fftw_plan_guru64_dft_c2r(rank, dims_o,howmany_rank, howmany_dims,  t3, t1, FFTW_ESTIMATE);

}
//------------------------------------------------------------------------------

void FFT_Analysis3D::setGF()
{
	double c1=2*mu0*(3*la0+2*mu0)/(la0+2*mu0);
	double c2=2*mu0;
	double c3=(c1-c2)/2.;
	double c4=(c2)/2.;

	double d1=1/(4*mu0);
	double d2=(la0+mu0)/(mu0*(la0+2*mu0));

	MatrixIdd q(21,4);

	/*[G1111 G1122 G1133 G1123 G1113 G1112
		    G2222 G2233 G2223 G2213 G2212
				 G3333 G3323 G3313 G3312
					   G3123 G3113 G3112
			sym              G3213 G3212
						          G1212]*/
			
  q<<   0, 0, 0, 0, 
		0, 0, 1, 1,  
		0, 0, 2, 2,
		0, 0, 1, 2, 
		0, 0, 0, 2,  
		0, 0, 0, 1,
		1, 1, 1, 1,	
		1, 1, 2, 2,
		1, 1, 1, 2,
		1, 1, 0, 2,
		1, 1, 0, 1,
		2, 2, 2, 2,
		2, 2, 1, 2,	
		2, 2, 0, 2,
		2, 2, 0, 1,
		1, 2, 1, 2,
		1, 2, 0, 2,
		1, 2, 0, 1,
		2, 0, 0, 2,
		2, 0, 0, 1,
		0, 1, 0, 1;  

	int x1, y1, z1;

	for (int x=0; x<Nx; x++)
		for (int y=0; y<Ny; y++)
			for (int z=0; z<Nz/2+1; z++)
			{  
				if(x<=(Nx-1)/2) {x1=x*1.;}
				else            {x1=-(Nx-x)*1.;}
						  
				if(y<=(Ny-1)/2) {y1=y*1.;}
				else            {y1=-(Ny-y)*1.;}

				if(z<=(Nz-1)/2) {z1=z*1.;}
				else            {z1=-(Nz-z)*1.;}

				double zet [3]={z1*1., y1*1., x1*1.}; // frequency vector for each voxel
				 
				int cc= (z+ y*(Nz/2+1)+x*((Nz/2+1)*Ny))*21;
				for (int s=0; s<21; s++)
				{
					int i= q(s,0);
					int j= q(s,1);
					int k= q(s,2);
					int h= q(s,3);

					if (strainBase==true)

						GF[cc+s] = (d1/nrm2(zet))*(dlt(k,i)*zet[h]*zet[j]+dlt(h,i)*zet[k]*zet[j]+dlt(k,j)*zet[h]\
				                  *zet[i]+dlt(h,j)*zet[k]*zet[i])-(d2*zet[i]*zet[j]*zet[k]*zet[h])/(nrm2(zet)*nrm2(zet));
					else
						GF[cc+s]=  c3*(dlt(i,j)-zet[i]*zet[j]/nrm2(zet))*(dlt(k,h)-zet[k]*zet[h]/nrm2(zet))+\
								   c4*(dlt(i,k)-zet[i]*zet[k]/nrm2(zet))*(dlt(j,h)-zet[j]*zet[h]/nrm2(zet))+\
								   c4*(dlt(i,h)-zet[i]*zet[h]/nrm2(zet))*(dlt(j,k)-zet[j]*zet[k]/nrm2(zet));

				}	       
			}

}

//------------------------------------------------------------------------------
void FFT_Analysis3D::multipCx()
{
	if (strainBase)

		for(int i=0; i<xyz; ++i)
			sigma.col(i)=C[RVE(i)]*epsilon.col(i);
	else 

		for(int i=0; i<xyz; ++i)
	        epsilon.col(i)=S[RVE(i)]*sigma.col(i);

}

//------------------------------------------------------------------------------
void FFT_Analysis3D::FFT_sigma()
{
	memcpy(t1, sigma.data(), 6*xyz*sizeof t1);
	fftw_execute(p_t1_t2_frw);
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::FFT_epsilon()
{
	memcpy(t1, epsilon.data(), 6*xyz*sizeof t1);

	if (mater_behav=="e")
		fftw_execute(p_t1_t2_frw);
	else 
		fftw_execute(p_t1_t4_frw);
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::FFTBCKWRD_sigma()
{
	fftw_execute(p_t3_t1_bkw);
	//sigma <- epsilon_c (epsilon_c in accelerated shceme)
	memcpy(sigma.data(), t1, 6*xyz*sizeof t1);
	sigma=sigma/xyz;
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::FFTBCKWRD_epsilon()
{
	fftw_execute(p_t3_t1_bkw);
	//epsilon <- sigma_c (sigma_c in accelerated shceme)
	memcpy(epsilon.data(), t1, 6*xyz*sizeof t1);
	epsilon=epsilon/xyz;		
}

//------------------------------------------------------------------------------
void FFT_Analysis3D::update()
{
	if (strainBase)

		for (int i=0; i<xyz; ++i)
			epsilon.col(i)= epsilon.col(i)+	C_acc[RVE(i)]*(sigma.col(i)-epsilon.col(i));

	else

		for (int i=0; i<xyz; ++i)
			sigma.col(i)= sigma.col(i)+	S_acc[RVE(i)]*(epsilon.col(i)-sigma.col(i));

}
//------------------------------------------------------------------------------
void FFT_Analysis3D::multiplyGF()
{
    MatrixIdd  w(6,6);
    w<< 0,  1,   2,  3,  4,   5,
        1,  6,   7,  8,  9,  10,
        2,  7,  11, 12, 13,  14,
        3,  8,  12, 15, 16,  17,
        4,  9,  13, 16, 18,  19,
        5, 10,  14, 17, 19,  20;

	if (mater_behav=="e") 
   
		for (int i=1; i<Nx*Ny*(Nz/2+1); i++)
		{  
			int u=i*6;
			int d=i*21;
			for(int j=0; j<6; j++)
			{ 
				t3[u+j][0]=-(GF[d+w(j,0)]*t2[u+0][0]  +GF[d+w(j,1)]*t2[u+1][0]  +GF[d+w(j,2)]*t2[u+2][0]
				            +GF[d+w(j,3)]*t2[u+3][0]*2+GF[d+w(j,4)]*t2[u+4][0]*2+GF[d+w(j,5)]*t2[u+5][0]*2);
				
				t3[u+j][1]=-(GF[d+w(j,0)]*t2[u+0][1]  +GF[d+w(j,1)]*t2[u+1][1]  +GF[d+w(j,2)]*t2[u+2][1]
				            +GF[d+w(j,3)]*t2[u+3][1]*2+GF[d+w(j,4)]*t2[u+4][1]*2+GF[d+w(j,5)]*t2[u+5][1]*2);
			}
		}

	else

		for (int i=1; i<Nx*Ny*(Nz/2+1); i++)
		{  
			int u=i*6;
			int d=i*21;
			for(int j=0; j<6; j++)
			{ 

				t3[u+j][0]=t4[u+j][0]-(GF[d+w(j,0)]*t2[u+0][0]  +GF[d+w(j,1)]*t2[u+1][0]  +GF[d+w(j,2)]*t2[u+2][0]
				                      +GF[d+w(j,3)]*t2[u+3][0]*2+GF[d+w(j,4)]*t2[u+4][0]*2+GF[d+w(j,5)]*t2[u+5][0]*2);
				
				t3[u+j][1]=t4[u+j][1]-(GF[d+w(j,0)]*t2[u+0][1]  +GF[d+w(j,1)]*t2[u+1][1]  +GF[d+w(j,2)]*t2[u+2][1]
				                      +GF[d+w(j,3)]*t2[u+3][1]*2+GF[d+w(j,4)]*t2[u+4][1]*2+GF[d+w(j,5)]*t2[u+5][1]*2);
			}
		}
			
//----------------------
	for (int i=0; i<6; ++i)
	{
		t3[i][0]=E(i)*xyz;
		t3[i][1]=0;
	}


}

//------------------------------------------------------------------------------
void FFT_Analysis3D::convergenceCheck()
{
    memcpy(t1, sigma.data(), 6*xyz*sizeof t1);
    fftw_execute(p_t1_t2_frw); 
    itNo+=1;
    double cc=0;
    double dd=sqrt((t2[0][0] *t2[0][0])+(t2[0][1] *t2[0][1])
                  +(t2[1][0] *t2[1][0])+(t2[1][1] *t2[1][1])
                  +(t2[2][0] *t2[2][0])+(t2[2][1] *t2[2][1])
                +2*(t2[3][0] *t2[3][0])+(t2[3][1] *t2[3][1])*2
                +2*(t2[4][0] *t2[4][0])+(t2[4][1] *t2[4][1])*2
                +2*(t2[5][0] *t2[5][0])+(t2[5][1] *t2[5][1])*2);

    
    double x1=0, y1=0, z1=0; 
    for (int x=0; x<Nx; x++)
        for (int y=0; y<Ny; y++)
            for (int z=1; z<Nz/2+1; z++)
            {  
                int i= (z+ y*(Nz/2+1)+x*((Nz/2+1)*Ny))*6;
                
                if(x<=(Nx-1)/2) {x1=x*1.;}
                else             {x1=-(Nx-x)*1.;}
                
                if(y<=(Ny-1)/2) {y1=y*1.;}
                else             {y1=-(Ny-y)*1.;}
                
                if(z<=(Nz-1)/2) {z1=z*1.;}
                else             {z1=-(Nz-z)*1.;}
                
                double zet [3]={z1, y1, x1}; // frequency vector or each voxel
                cc+= sqrt(abs ((t2[i+0][0]*t2[i+0][0]+t2[i+0][1]*t2[i+0][1])*zet[0]*zet[0]
                              +(t2[i+1][0]*t2[i+1][0]+t2[i+1][1]*t2[i+1][1])*zet[1]*zet[1]
                              +(t2[i+2][0]*t2[i+2][0]+t2[i+2][1]*t2[i+2][1])*zet[2]*zet[2]
                              +(t2[i+3][0]*t2[i+3][0]+t2[i+3][1]*t2[i+3][1])*((zet[1]+zet[1])+(zet[2]+zet[2]))
                              +(t2[i+4][0]*t2[i+4][0]+t2[i+4][1]*t2[i+4][1])*((zet[0]+zet[0])+(zet[2]+zet[2]))
                              +(t2[i+5][0]*t2[i+5][0]+t2[i+5][1]*t2[i+5][1])*((zet[0]+zet[0])+(zet[1]+zet[1]))));
    
            }


    avr_eqVar=cc/dd/xyz;

    cout<<"Itter="<<itNo<<"   eqVar= "<<avr_eqVar<<"   ";
	
	if (strainBase==true)
	{

        diff_sum_field=abs((sigma.sum()-sum_field_old)/(xyz*(0.5*C[0](0,0)+0.5*C[1](0,0))*E.sum()));

		diff_avr_eqVar=abs(avr_eqVar-aver_eqVarOld);

		cout<<"diff_sum_field= "<<diff_sum_field<<"    "<<"diff_avr_eqVar= "<<diff_avr_eqVar<<"    "<<endl<<endl;
	}
	else if (stressBase==true)
	{

        diff_sum_field=abs((epsilon.sum()-sum_field_old)/(xyz*(0.5*S[0](0,0)+0.5*S[1](0,0))*E.sum()));

		diff_avr_eqVar=abs(avr_eqVar-aver_eqVarOld);

		cout<<"diff_sum_field= "<<diff_sum_field<<"    "<<"diff_avr_eqVar= "<<diff_avr_eqVar<<"    "<<endl<<endl;
	}

	aver_eqVarOld=avr_eqVar;
}


//------------------------------------------------------------------------------
void FFT_Analysis3D::material (MatrixDdd &mat, MatrixIdd &RVE1, int noPhases1)
{ 
    strainBase=false;
	stressBase=false;

	noPhases=noPhases1;

	El=new double [noPhases];
	v =new double [noPhases];
	mu=new double [noPhases];
	la=new double [noPhases];
	K =new double [noPhases];
	H =new double [noPhases];
	n =new double [noPhases];
	Sy=new double [noPhases];
	P =new double [noPhases];
	Y0=new double [noPhases];

   	C=new MatrixDdd [noPhases1];
   	S=new MatrixDdd [noPhases1];

	for (int i=0; i<noPhases; ++i)
	{
		El[i]=mat(i,0);
		 v[i]=mat(i,1);
		Sy[i]=mat(i,2);	
		 H[i]=mat(i,3);
		 n[i]=mat(i,4);
		 P[i]=mat(i,5);
		Y0[i]=mat(i,6);

		mu[i]=El[i]/(1+v[i])/2.;
		la[i]=v[i]*El[i]/(1+v[i])/(1-2*v[i]);
		K[i]=la[i]+(2./3.)*mu[i];

	    C[i].resize(6,6);	

		C[i]<< 	la[i]+2*mu[i],la[i],          la[i],     0,     0,       0,
				la[i],        la[i]+2*mu[i],  la[i],     0,     0,       0,
				la[i],        la[i],  la[i]+2*mu[i],     0,     0,       0,
				0 ,           0,               0,   2*mu[i],    0,       0, 
				0 ,           0,               0,        0,  2*mu[i],    0,
				0 ,           0,               0,        0,     0,  2*mu[i];

		S[i].resize(6,6);
		S[i]=C[i].inverse();
	}
//------------------------
	RVE=RVE1;  
	double ctr;
 	vol_fractions=new double [noPhases];

	for (int k=0; k<noPhases; ++k)
	{
		ctr=0;
		for (int i=0; i<xyz; ++i)
			if (int(RVE(i))==k)
			{	ctr+=1;}

		vol_fractions[k]=ctr/xyz;		
	}

	double mu_max, mu_min, la_max, la_min, K_min, K_max, vol_frac_max, index_mu_max, temp[4];
    max_min_Array(mu, noPhases, temp);             
	mu_max=temp[0];  	mu_min=temp[2];  
	index_mu_max=int(temp[1]);

	max_min_Array(la , noPhases, temp);             
	la_max =temp[0],   la_min =temp[2];

    max_min_Array(K , noPhases, temp);             
	K_max =temp[0],   K_min =temp[2];

    max_min_Array(vol_fractions, noPhases, temp);   vol_frac_max=temp[0]; 


	if (mater_behav=="e")
	{		
		double K0;
		K0 =sqrt(K_min*K_max);
		mu0=sqrt(mu_min*mu_max);
		la0=K0-mu0*2./3.;
		
		if (vol_fractions[int(index_mu_max)]==vol_frac_max) 
		{
			strainBase=true;
			std::cout<< " Strain-Based Algorithm is Running"<<std::endl;
		}
		else
		{
			stressBase=true; 
			std::cout<< " Stress-Based Algorithm is Running"<<std::endl;
		}
	}

	else
	{
		strainBase=true;

		double alfa=1;
		mu0=alfa*(mu_min+mu_max)/2;
		la0=alfa*(la_min+la_max)/2;
	}


	C0.resize(6,6);
	C0<<la0+2*mu0, la0,         la0,     0,      0,    0,
		la0,        la0+2*mu0,  la0,     0,      0,    0,
		la0,        la0,  la0+2*mu0,     0,      0,    0,
		0 ,         0,          0,   2*mu0,      0,    0, 
		0 ,         0,          0,       0,  2*mu0,    0,
		0 ,         0,          0,       0,      0, 2*mu0;

	S0.resize(6,6);
	S0=C0.inverse();

 	C_acc=new MatrixDdd [noPhases];
   	S_acc=new MatrixDdd [noPhases];

	for (int i=0; i<noPhases; ++i)
	{
		C_acc[i].resize(6,6);
		C_acc[i]=2*(C[i]+C0).inverse()*C0;

		S_acc[i].resize(6,6);
		S_acc[i]=2*(S[i]+S0).inverse()*S0;
	}

	setGF();
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::setLoad(MatrixDdd &E1, MatrixDdd &X01)
{ 
    E=E1;
    X0=X01;
	if (strainBase==true)
		for (int i=0; i<xyz; ++i)
		    epsilon.col(i)=E1; 

	else if (stressBase==true)
		for (int i=0; i<xyz; ++i)
		    sigma.col(i)=E1; 
		

}
//------------------------------------------------------------------------------
double FFT_Analysis3D::prod(MatrixDdd &a, MatrixDdd &b)
{
  return a(0)*b(0)+a(1)*b(1)+a(2)*b(2)+2*a(3)*b(3)+2*a(4)*b(4)+2*a(5)*b(5);

}
//------------------------------------------------------------------------------
void FFT_Analysis3D::old_iter_fields_update()
{

	for (int i=0; i<6; ++i)
	{
		sigma_iter_old(i)=sigma.row(i).sum()/xyz;
		epsilon_iter_old(i)=epsilon.row(i).sum()/xyz;
	}
	
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::updateE()
{ 
	old_iter_fields_update();
	double k;

	MatrixDdd x;
	MatrixDdd y;
	x=S0*sigma_iter_old-epsilon_iter_old;
	y=S0*X0;
	k=(dE.sum()*incNo+prod(x, X0))/prod(y, X0);
	E=k*S0*X0-S0*sigma_iter_old+epsilon_iter_old;
}

//------------------------------------------------------------------------------
void FFT_Analysis3D::sigma_eff_calculator(int _i)
{
	double sm=0;

	sm=(sigma(0,_i)+sigma(1,_i)+sigma(2,_i))/3.;
	sigma_d=sigma.col(_i)-sm*I;
	sigma_eff(_i)=sqrt(3./2.)*sqrt(sigma_d(0)*sigma_d(0)  +sigma_d(1)*sigma_d(1)  +sigma_d(2)*sigma_d(2)+\
					             2*sigma_d(3)*sigma_d(3)+2*sigma_d(4)*sigma_d(4)+2*sigma_d(5)*sigma_d(5));

}

//------------------------------------------------------------------------------
//constitutive 
void FFT_Analysis3D::constitutive()
{  
	if (mater_behav=="ep")
		elastoPlastic();

	else if(mater_behav=="ed")
		elastoDamage();
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::elastoDamage()
{
	sum_field_old=sigma.sum();
	MatrixDdd temp(6,1), temp1(6,1);
	int k;
	for (int i=0; i<xyz; ++i)
	{	
		k=RVE(i);
		temp1=epsilon.col(i);
		temp=C[k]*temp1; temp(3)*=2; temp(4)*=2; temp(5)*=2;
		Y(i)=sqrt((temp.transpose()*temp1)(0));
			
		d_old(i)=d(i);

		if (Y(i)>Y0[k])
			d(i)=max((1-exp(-P[k]*(Y(i)-Y0[k]))), d_old(i));		

		sigma.col(i)=(1-d(i))*C[k]*(epsilon.col(i));
	}
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::elastoPlastic()
{   
	sum_field_old=sigma.sum();

	double  x, x_o=0, res, sm=0; int k;
	for(int i=0; i<xyz; ++i)
	{
		k=RVE(i);

      	sigma.col(i)=sigma_old.col(i)+C[k]*(epsilon.col(i)-epsilon_old.col(i));

		sigma_eff_calculator(i);

	 	phi(i)=sigma_eff(i)-Sy[k]-H[k]*pow(epsilon_p_eff(i), n[k]);
	
		if (phi(i)>0)
		{
			//Linear Hardening	
			//x=(sigma_eff(i)-Sy[k]-H[k]*epsilon_p_eff(i))/(3*mu[k]+H[k]);

			// Newton-Raphson
			x_o=0;
			while (true)
			{	
				//x	:d_epsilon_p	
				x=x_o-(sigma_eff(i)-3*mu[k]*x_o-Sy[k]-H[k]*pow((epsilon_p_eff(i)+x_o),n[k]))\
                                 /(-3*mu[k]-H[k]*n[k]*pow((epsilon_p_eff(i)+x_o),n[k]-1));
				x_o=x;
				res=(sigma_eff(i)-3*mu[k]*x-Sy[k]-H[k]*pow((epsilon_p_eff(i)+x),n[k]));
				
				//cout<< (x)<<"********"<<sigma_eff(i)<<"********"<<epsilon_p_eff(i)<<"******"<<RVE(i)<<"*****"<<res<<endl;
				
				if (isnan(x))
				{throw; cout<<" Newton-Raphson failed to converge in voxel i= "<<i<<endl;}

				if (abs (res)<1e-5)
					break;

			}

			d_ep_p(i)=x;
			sigma.col(i)-=(1.5*x/sigma_eff(i))*C[k]*sigma_d;
			sigma_eff_calculator(i);

		}

		else
			d_ep_p(i)=0;

	}
		
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::update_old_Variable()
{
	if (mater_behav=="ep")
	{
		for (int i=0; i<xyz; ++i)
		{
			if (phi(i)>0)
			{
				epsilon_p_eff(i)+=d_ep_p(i);
				//epsilon_p.col(i)+=(1.5*d_ep_p(i)/sigma_eff(i))*sigma_d.col(i);	
			}

			sigma_old.col(i)=sigma.col(i);
		}
	}

	MatrixDdd temp;
	temp.resize(6,1);
	for (int i=0; i<xyz; ++i)
	{
		temp=epsilon_old.col(i);	
		epsilon_old.col(i)=epsilon.col(i);
		epsilon.col(i)+=(epsilon.col(i)-temp);
	}

	
}

//------------------------------------------------------------------------------
void FFT_Analysis3D::export_homogen_var(int _i)
{
	int t;
	for (int j=0; j<6; ++j)
		if (X0(j)!=0)
		{	t=j; break;}

	prop_homogenized.resize(nInc,9);
	prop_homogenized(_i,0)=epsilon.row(t).sum()/xyz;

	for (int j=1; j<7;++j)
		prop_homogenized(_i,j)=sigma.row(j-1).sum()/xyz;

	if (mater_behav=="ep")
	{	
		prop_homogenized(_i,7)=sigma_eff.sum()/xyz;
		prop_homogenized(_i,8)=epsilon_p_eff.sum()/xyz;
	}

	if (mater_behav=="ed")
	{	
		prop_homogenized(_i,7)=d.sum()/xyz;
		prop_homogenized(_i,8)=0;
	}
}

//------------------------------------------------------------------------------
void FFT_Analysis3D::runAnalysis(double tol_field, double tol_eq, int nInc1)
{ 
	if (mater_behav=="e")
	{
		multipCx();
		while (true)
			{
				
				if (strainBase==true)
				{   //sigma<-tau
					sum_field_old=sigma.sum();
					sigma=sigma-C0*epsilon;
				}

				if (stressBase==true)
				{	//epsilon<-tau
					sum_field_old=epsilon.sum();
					epsilon=epsilon-S0*sigma;
				}

				if (strainBase==true)
					FFT_sigma();
				if (stressBase==true)
					FFT_epsilon();	
				
				multiplyGF();

				if (strainBase==true)
					FFTBCKWRD_sigma();
				if (stressBase==true)
					FFTBCKWRD_epsilon();

				update();
				multipCx();
				convergenceCheck();

				if (itNo>1 && abs(diff_avr_eqVar)<tol_eq  &&  abs(diff_sum_field)<tol_field )
					break;
				if (itNo>2000)
				{
					std::cout<<" iteration excede from the allowable number"<<std::endl;
					break;
				}

			}
	}
//-----------------------------------------------------------------------------
	else if (mater_behav=="ep" or mater_behav=="ed")
	{
		nInc=nInc1;
		dE=E/(nInc*1.);  E=dE;  setLoad(E, X0);
		incNo=0; int nx=0;
		for (int _i=0; _i<nInc; ++_i)
		{
			incNo+=1;
			itNo=0;
			cout<<"**********************"<<endl;
			cout<<"Increment= "<<_i<<endl;
			cout<<"**********************"<<endl;

			constitutive();
			while (true)
			{
				FFT_sigma();
				FFT_epsilon();
				updateE();
				multiplyGF();
				FFTBCKWRD_epsilon();
				constitutive();
				convergenceCheck();

				if (itNo>1 && abs(diff_avr_eqVar)<tol_eq  &&  abs(diff_sum_field)<tol_field)
					break;
				if (itNo>2000)
				{
					std::cout<<" iteration excede from the allowable number"<<std::endl;
					break;
				}
			}

			if (mater_behav=="ed")
				if (_i!=0){

					if (X0(0)!=0 and prop_homogenized(_i,1)<prop_homogenized(_i-1,1))
						nx+=1;
					if (X0(1)!=0 and prop_homogenized(_i,2)<prop_homogenized(_i-1,2))
						nx+=1;
					if (X0(2)!=0 and prop_homogenized(_i,3)<prop_homogenized(_i-1,3))
						nx+=1;
					if (X0(3)!=0 and prop_homogenized(_i,4)<prop_homogenized(_i-1,4))
						nx+=1;
					if (X0(4)!=0 and prop_homogenized(_i,5)<prop_homogenized(_i-1,5))
						nx+=1;
					if (X0(5)!=0 and prop_homogenized(_i,6)<prop_homogenized(_i-1,6))
						nx+=1;
				}

			if (nx==3)
				{nInc=_i; break;}	
						
			export_homogen_var(_i);
			update_old_Variable();
			E+=dE;
		}
	}
}
//------------------------------------------------------------------------------
double FFT_Analysis3D::homogenization(double tol_field, double tol_eq, string outputName)
{  
	mater_behav="e";
    ofstream homo(outputName+".txt");
    epHomo.resize(6,6);
    CHomo.resize(6,6);
    SHomo.resize(6,6);
    sigHomo.resize(6,6);


    for (int w=0; w<6; ++w)
	{
		std::cout<<"load case: "<<w+1<<std::endl;
	    itNo=0;
	    E.setZero();

		if (strainBase)
			E(w)=0.0005;
		else
			E(w)=1.e6;

        setLoad(E, E);
	    runAnalysis(tol_field, tol_eq, 1);
	    
	    for (int t=0; t<6; ++t)
	    {
			sigHomo(w,t)=  sigma.row(t).sum()/xyz;
			epHomo(t,w) =epsilon.row(t).sum()/xyz;
	    }
	    

	}

    for (int t=0; t<6; ++t)
	{
		if(strainBase)
		{
			CHomo.col(t)=epHomo.inverse()* sigHomo.col(t);
			SHomo=CHomo.inverse();
		}

		else
		{
			SHomo.col(t)=sigHomo.inverse()* epHomo.col(t);
			CHomo=SHomo.inverse();
		}
	}    
  
	E11=1./SHomo(0,0); 	E22=1./SHomo(1,1);  	E33=1./SHomo(2,2); 
/*-------------------------------------------------------------------------------
%Hashin-shtrikman Bounds For Effective Elastic Properties of the n-phases Composites
%Based on the:
%"A VARIATIONAL APPROACH TO THE THEORY OFTHE ELASTIC BEHAVIOUR OF MULTIPHASE
% MATERIALS, JMPS, 1963"
---------------------------------------------------------------------------------*/
	double temp[4];
	max_min_Array(mu, noPhases, temp);  double Gn =temp[0];   double G1 =temp[2];
	max_min_Array(K , noPhases, temp);  double Kn =temp[0];   double K1 =temp[2]; 

	double alfa1=-3/(3*K1+4*G1);
	double alfan=-3/(3*Kn+4*Gn);
	double beta1=-(3*(K1+2*G1))/((3*K1+4*G1)*5*G1);
	double betan=-(3*(Kn+2*Gn))/((3*Kn+4*Gn)*5*Gn);

	double A1=0, An=0,  B1=0, Bn=0;

	for (int i=0; i<noPhases; ++i)
	{
		if (K[i]!=K1)
			A1+=vol_fractions[i]/((1./(K[i]-K1))-alfa1);
		if (K[i]!=Kn)
			An+=vol_fractions[i]/((1./(K[i]-Kn))-alfan);

		if (mu[i]!=G1)
			B1+=vol_fractions[i]/((0.5/(mu[i]-G1))-beta1);
		if (mu[i]!=Gn)
			Bn+=vol_fractions[i]/((0.5/(mu[i]-Gn))-betan);
	}

	double K_min=K1+A1/(1+alfa1*A1);
	double K_max=Kn+An/(1+alfan*An);

	double G_min=G1+0.5*B1/(1+beta1*B1);
	double G_max=Gn+0.5*Bn/(1+betan*Bn);

	double E_min=9*K_min*G_min/(3*K_min+G_min);
	double E_max=9*K_max*G_max/(3*K_max+G_max);

    homo<<"*************************************************************************"<<endl<<endl;
 	
    homo<<"        F F T     H O M O G E N I Z A T I O N    B Y   M A T 3 D"<<endl<<endl;
    if (strainBase)
		homo<< "                      Stran-Based Algorithm Is Used"<<endl;
	else 
		homo<< "                      Stress-Based Algorithm Is Used"<<endl;

	homo<<"                           Nx="<<Nx<<"   "<<"Ny="<<Ny<<"   "<<"Nz="<<Nz<<endl<<endl;

    homo<<"      ----------------------------------------------------------------"<<endl;
	for (int i=0; i<noPhases; ++i)
	{
		string k=to_string(i)+" = ";
		homo<<"                      Volume Fraction Phase "<<k<<vol_fractions[i]<<endl;
		homo<<"     E"<<k<<El[i]<<"     "<<"v"<<k<<v[i]<<"     "<<"K"<<k<<K[i]<<"     "<<"G"<<k<<mu[i]<<"     "<<endl<<endl;

	}
    homo<<"      ----------------------------------------------------------------"<<endl<<endl;
	 homo<<"                        Hashin-Shtrikman Bounds"<<endl;
	homo<<"         E_min="<<E_min<<"   "<<"K_min="<<K_min<<"   "<<"G_min="<<G_min<<endl;
	homo<<"         E_max="<<E_max<<"   "<<"K_max="<<K_max<<"   "<<"G_max="<<G_max<<endl<<endl<<endl;


	homo<<"************************************************************************"<<endl<<endl;



    homo<<"Homogenized Elastic Modules:"<<endl;
    homo<<"E11="<<1/SHomo(0,0)<<"     "<<"E22="<<1/SHomo(1,1)<<"     "<<"E33="<<1/SHomo(2,2)<<"     "<<endl<<endl;

    homo<<endl<<"C_Homogenized="<<endl<<endl;
    homo<<CHomo<<endl;
	homo<<"      ----------------------------------------------------------------"<<endl;
    homo<<endl<<"S_Homogenized="<<endl<<endl;
    homo<<SHomo<<endl;
	homo<<"      ----------------------------------------------------------------"<<endl;

	for(int i=0; i<noPhases; ++i)
	{
		homo<<endl<<"C"<<to_string(i)<<"="<<endl;
		homo<<C[i]<<endl;
		homo<<"      ----------------------------------------------------------------"<<endl;
		homo<<endl<<"S"<<to_string(i)<<"="<<endl;
		homo<<S[i]<<endl;
    	homo<<"      ----------------------------------------------------------------"<<endl;
	}

    //homo<<endl<<endl<<"--------------------------------------------------------------------"<<endl<<endl;
    homo<<endl<<"Homogenized Load Cases (Each row is one load case i.e (E11 E22 E33 E32 E31 E12) and Eij is summation of strain components in in  all voxels )"<<endl<<endl;
    homo<<epHomo<<endl;
    homo<<endl<<"Homogenized stess (Each coulemn is (SS11 SS22 SS33 SS32 SS31 SS12) and SSij is summation of stress components in in all voxels )"<<endl<<endl;
    homo<<sigHomo<<endl;

    return 1./SHomo(0,0)+1./SHomo(1,1)+1./SHomo(2,2);
}
//------------------------------------------------------------------------------
void FFT_Analysis3D::printer (std::string outFldr)
{
    mkdir(outFldr.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    
	std::vector<string> fileNames={ "ep1.txt",  "ep2.txt",  "ep3.txt",  "ep4.txt",  "ep5.txt", "ep6.txt", 
			                       "sig1.txt", "sig2.txt", "sig3.txt", "sig4.txt", "sig5.txt", "sig6.txt", "ep_p"};
	if (mater_behav=="ed")
		fileNames[12]="d";

	int j=0;
	for (auto &elem : fileNames) 
	{   
		ofstream output100(outFldr+"/"+elem);
 
		for (int i=0; i<Nx*Ny*Nz; i++)
		{ 
			if (j<6)
				output100 << epsilon(j,i)<<endl;

			else if (j<12 && j>=6)
		   	 	output100 << sigma(j-6,i)<<endl;

			else if (j==12 && mater_behav=="ep")
		   	 	output100 << epsilon_p_eff(i)<<endl;

			else if (j==12 && mater_behav=="ed")
		   	 	output100 << d(i)<<endl;
				
		}
		output100.close();
		j++;
	}


	if (mater_behav!="e")
	{
		ofstream op002("epsilon_sigma.txt");
 
		if (mater_behav=="ep")
			op002<<"epsilon"<<"  "<<"sigma1"<<"   "<<"sigma2"<<"   "<<"sigma3"<<"   "<<"tau23"<<"   "<<"tau13"<<\
				   "   "<<"tau12"<<"   "<<"sigma_eff"<<"   "<<"epsilon_p_eff"<<endl;
		else
			op002<<"epsilon"<<"  "<<"sigma1"<<"   "<<"sigma2"<<"   "<<"sigma3"<<"   "<<"tau23"<<"   "<<"tau13"<<\
				   "   "<<"tau12"<<"   "<<"d"<<"   "<<"epsilon_p_eff"<<endl;

		op002<<0<<"  "<<0<<"  "<<0<<"  "<<0<<"  "<<0<<"  "<<0<<"  "<<0<<"  "<<0<<"  "<<0<<endl;

		for (int i=0; i<nInc; ++i)
			op002<<prop_homogenized(i, 0)<<"  "<<prop_homogenized(i, 1)<<"  "<<prop_homogenized(i, 2)<<"  "<<\
				   prop_homogenized(i, 3)<<"  "<<prop_homogenized(i, 4)<<"  "<<prop_homogenized(i, 5)<<"  "<<\
				   prop_homogenized(i, 6)<<"  "<<prop_homogenized(i, 7)<<"  "<<prop_homogenized(i, 8)<<endl;
	}
}
//vtiwriter---------------------------------------------------------------------------

void FFT_Analysis3D::vtiWriter(std::string fileName)
{	
	  vtkSmartPointer<vtkImageData> imageData =
	  vtkSmartPointer<vtkImageData>::New();
	  imageData->SetDimensions(Nx, Ny, Nz);//
#if VTK_MAJOR_VERSION <= 5
	  imageData->SetNumberOfScalarComponents(14);
	  imageData->SetScalarTypeToDouble();
#else
	  imageData->AllocateScalars(VTK_DOUBLE, 14);
#endif
	  int* dims = imageData->GetDimensions();


	  // Fill every entry of the image data with "2.0"
	  for (int z = 0; z < dims[2]; z++)
		for (int y = 0; y < dims[1]; y++)
		  for (int x = 0; x < dims[0]; x++)
		  {
			double* pixel= static_cast<double*>(imageData->GetScalarPointer(x,y,z));
            int len=x+y*dims[0]+z*dims[0]*dims[1];
			for(int i=0; i<6; ++i)
			{
				pixel[i] =    epsilon(i, len);
		        pixel[i+6] =  sigma(i+6, len);
			}
				if (mater_behav=="ep"){
					pixel[12] =  epsilon_p_eff(len);
					pixel[13] =  sigma_eff(len);
				}

				if (mater_behav=="ed"){
					pixel[12] =  d(len);
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

//------------------------------------------------------------------------------
FFT_Analysis3D::~FFT_Analysis3D()
{
    delete [] t1;
    delete [] t2;
	delete [] t3;
	delete [] t4;

	delete [] El;
	delete [] v;
	delete [] mu;
	delete [] la;
	delete [] K;
	delete [] H ;
	delete [] n;
	delete [] Sy;
	delete [] P ;
    delete [] Y0;
	//delete [] mater_behav;
   	delete [] C;
   	delete [] S;

   	delete [] C_acc;
   	delete [] S_acc;

    delete [] GF;

	delete [] dims_o;
	delete [] dims;
	delete [] howmany_dims;

    fftw_destroy_plan(p_t1_t2_frw);
    fftw_destroy_plan(p_t1_t4_frw);
    fftw_destroy_plan(p_t3_t1_bkw);

    fftw_cleanup();
    fftw_cleanup_threads();  

    delete [] vol_fractions; 
    epsilon.resize(0,0);
    sigma.resize(0,0);
}
