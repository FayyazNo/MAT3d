#include "FFT_analysis2D.hpp"
 
FFT_Analysis2D::FFT_Analysis2D (int Nx1, int Ny1)
{  
    Nx=Nx1;
    Ny=Ny1;
    xy=Nx*Ny;
    t1  =   new fftw_complex [3*xy];
    t2  =   new fftw_complex [3*xy];

    GF= new double [6*xy];
    cout<<endl<<6*xy;
    RVE.resize(1, xy);
    epsilon.resize(3, xy);
    sigma.resize(3, xy);
    sigmaOld.resize(3, xy);
    tau.resize(3, xy);
    E.resize(3,1);

    for(int i=0; i<3*xy; i++)
     {
      t1[i][1]=0;
      t2[i][1]=0;
     //b1[i][1]=0;
      //b2[i][1]=0;
     }

    int rank=2;
    fftw_iodim64 *dims;
    dims= new fftw_iodim64[rank];

    if(dims==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    dims[0].n=Nx;
    dims[0].is=Ny*3;
    dims[0].os=Ny*3;
    dims[1].n=Ny;
    dims[1].is=3;
    dims[1].os=3;

    int howmany_rank=1;
    fftw_iodim64 *howmany_dims;
    howmany_dims=new fftw_iodim64 [howmany_rank];

    if(howmany_dims==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    howmany_dims[0].n=3;
    howmany_dims[0].is=1;
    howmany_dims[0].os=1;

    p3=fftw_plan_guru64_dft(rank, dims,howmany_rank, howmany_dims,t1, t2,FFTW_FORWARD, FFTW_ESTIMATE);
    p4=fftw_plan_guru64_dft(rank, dims,howmany_rank, howmany_dims,t1, t2,FFTW_BACKWARD, FFTW_ESTIMATE);
}
 
void FFT_Analysis2D::setGF()
{
    double d1=1/(4*mu0);
    double d2=(la0+mu0)/(mu0*(la0+2*mu0));
    MatrixIdd q(6,4);

    q<<   0, 0, 0, 0, 
    0, 0, 1, 1,    
    0, 0, 0, 1,
    1, 1, 1, 1,	
    1, 1, 0, 1,
    0, 1, 0, 1;
    
    int x1, y1;
    for (int x=0; x<Nx; x++)
     {
        for (int y=0; y<Ny; y++)
         {
          //cout<<x<<"    "<<y<<endl;
          int cc= (y+x*Ny)*6;
          if(x<=(Nx-1)/2) {x1=x*1.;}
          else            {x1=-(Ny-x)*1.;}
               
          if(y<=(Ny-1)/2) {y1=y*1.;}
          else            {y1=-(Ny-y)*1.;}
        
          double zet [2]={x1*1.0, y1*1.0}; // frequency vector or each voxel
        
          for (int s=0; s<6; s++)
           {
            int i= q(s,0);
            int j= q(s,1);
            int k= q(s,2);
            int h= q(s,3);
            GF[cc+s] = (d1/nrm2(zet))*(dlt(k,i)*zet[h]*zet[j]+dlt(h,i)*zet[k]*zet[j]+dlt(k,j)*zet[h]*zet[i]+dlt(h,j)*zet[k]*zet[i])-\
                                            (d2*zet[i]*zet[j]*zet[k]*zet[h])/(nrm2(zet)*nrm2(zet));
            
            }
         }
      } 
} 
        

void FFT_Analysis2D::multipCx()
{
    sigmaOld=sigma;
    for(int i=0; i<xy; ++i)
    {
        if (RVE(i)==1)
          sigma.col(i)=Cm*epsilon.col(i);
        else
          sigma.col(i)=Ci*epsilon.col(i); 
    }
}
  
void FFT_Analysis2D::FFT_tau()
{
    memcpy(*t1, tau.data(), 2*3*xy*sizeof t1);
    fftw_execute(p3);
}

void FFT_Analysis2D::multiplyGF()
{
    MatrixIdd  w(3,3);
    w<< 0,  1,   2,  
        1,  3,   4,  
        2,  4,   5; 

    for (int i=1; i<xy; i++)
    {  
     int u=i*3;
     int d=i*6;
     for(int j=0; j<3; j++)
        { 
         t1[u+j][0]=-(GF[d+w(j,0)]*t2[u+0][0]  +GF[d+w(j,1)]*t2[u+1][0]  +GF[d+w(j,2)]*t2[u+2][0]*2);
         t1[u+j][1]=-(GF[d+w(j,0)]*t2[u+0][1]  +GF[d+w(j,1)]*t2[u+1][1]  +GF[d+w(j,2)]*t2[u+2][1]*2);
        }

    }
    t1[0][0]=E(0).real()*xy;
    t1[0][1]=E(0).imag()*xy;
               
    t1[1][0]=E(1).real()*xy;
    t1[1][1]=E(1).imag()*xy;
                   
    t1[2][0]=E(2).real()*xy;
    t1[2][1]=E(2).imag()*xy;

}

void FFT_Analysis2D::FFTBCKWRD_epsilon()
{
    fftw_execute(p4);
    memcpy(epsilon.data(), *t2, 2*3*xy*sizeof t1);
    epsilon=epsilon/xy;
}

void FFT_Analysis2D::convergenceCheck()
{
    memcpy(*t1, sigma.data(), 2*3*xy*sizeof t1);
    fftw_execute(p3); 
    itNo+=1;
    double cc=0;
    double dd=sqrt((t2[0][0] *t2[0][0])+(t2[0][1] *t2[0][1])
                  +(t2[1][0] *t2[1][0])+(t2[1][1] *t2[1][1])
                +2*(t2[2][0] *t2[2][0])+(t2[2][1] *t2[2][1])*2);
               

    double x1=0, y1=0;
    for (int x=0; x<Nx; x++)
	    {
	     for (int y=0; y<Ny; y++)
	      {
		      int i= (y+ x*Ny)*3;

              if(x<=(Nx-1)/2) {x1=x*1.;}
              else            {x1=-(Nx-x)*1.;}
                 
              if(y<=(Ny-1)/2) {y1=y*1.;}
              else            {y1=-(Ny-y)*1.;}

		      double zet [2]={x1, y1}; // frequency vector or each voxel
 
               cc+= sqrt ((t2[i+0][0] *t2[i+0][0])*zet[0]*zet[0]
                         +(t2[i+1][0] *t2[i+1][0])*zet[1]*zet[1]
                         +(t2[i+2][0] *t2[i+2][0])*(zet[1]+zet[2])*(zet[1]+zet[2]));
   
		   }
	     }
    
    convrg=cc/dd/xy;
    cout<<"Itter="<<itNo<<"---------"<<"Equilibrium is ="<<convrg<<endl<<endl;
    cout<<"Diffrence in two steps="<<(sigma.sum().real()-sigmaOld.sum().real())/xy<<"    "<<endl<<endl;

}

void FFT_Analysis2D::material (double *properties, MatrixDdd &RVE1)
{ 
    Cm.resize(3,3);
    Ci.resize(3,3);
    C0.resize(3,3);

    Em=properties[0];
    vm=properties[1];
    Ei=properties[2];
    vi=properties[3];

    mum=Em/(1+vm)/2;
    lam=vm*Em/(1+vm)/(1-2*vm);
    mui=Ei/(1+vi)/2;
	lai=vi*Ei/(1+vi)/(1-2*vi);
			
	Cm<<   lam+2*mum,   lam,           0,    
	       lam,         lam+2*mum,     0,    
	       0 ,          0,             2*mum;    
	
    Ci<<   lai+2*mui,   lai,           0,    
	       lai,         lai+2*mui,     0,    
           0 ,          0,             2*mui;    
						

	double alf=1;
	C0=alf*(Cm+Ci)/2.;

	la0=alf*(lai+lam)/2.;
	mu0=alf*(mui+mum)/2.;

    RVE=RVE1;
}

void FFT_Analysis2D::setLoad(MatrixCdd &E1)
{ 
    E=E1;
      
    for (int i=0; i<xy; ++i)
      epsilon.col(i)=E1; 

}

void FFT_Analysis2D::runAnalysis()
{ 
    setGF();

    for (int i=0; i<25; ++i)
    {
        multipCx();
        tau=sigma-C0*epsilon;
        FFT_tau();
        multiplyGF();
        FFTBCKWRD_epsilon();
        convergenceCheck();
    }
}


void FFT_Analysis2D::homogenization()
{  
    ofstream homo("Homogenized_Properties");
    epHomo.resize(3,3);
    CHomo.resize(3,3);
    SHomo.resize(3,3);
    Si.resize(3,3);
    Sm.resize(3,3);
    sigHomo.resize(3,3);

    for (int w=0; w<3; ++w)
    {
       E.setZero();
       E(w)=0.005;
       runAnalysis();
   
       double e1=0, e2=0, e3=0;
       double s1=0, s2=0, s3=0;

       for (int t=0; t<xy; ++t)
       {
         e1+=epsilon(0,t).real();
         e2+=epsilon(1,t).real();
         e3+=epsilon(2,t).real();

         s1+=sigma(0,t).real();
         s2+=sigma(1,t).real();
         s3+=sigma(2,t).real();

       }
      
       epHomo(w,0)=e1/xy;
       epHomo(w,1)=e2/xy;
       epHomo(w,2)=e3/xy;

       sigHomo(0,w)=s1/xy;
       sigHomo(1,w)=s2/xy;
       sigHomo(2,w)=s3/xy;
    }

    CHomo.col(0)=epHomo.inverse()* sigHomo.col(0);
    CHomo.col(1)=epHomo.inverse()* sigHomo.col(1);
    CHomo.col(2)=epHomo.inverse()* sigHomo.col(2);

     
    homo<<"Volume Fraction of the White Phase (or 1)="<<RVE.sum()/xy<<endl<<endl;
    
    homo<<"E phase matrix ( or 1)="<<Em<<"     "<<"v phase  matrix (or 1)="<<vm<<endl;
    homo<<"E phase inclusion(or 0)="<<Ei<<"     "<<"v Phase  iinclusion(or 0)="<<vi<<endl;
    homo<<"Nx="<<Nx<<"   "<<"Ny="<<Ny<<"   "<<endl<<endl;
    
    homo<<endl<<"CHomogenized="<<endl<<endl;
    homo<<CHomo<<endl;
    
    homo<<endl<<"Cm="<<endl<<endl;
    homo<<Cm<<endl;
    
    homo<<endl<<"Ci="<<endl<<endl;
    homo<<Ci<<endl;
    
    homo<<endl<<"SHomogenized="<<endl<<endl;
    homo<<CHomo.inverse()<<endl;
    
    homo<<endl<<"Sm="<<endl<<endl;
    homo<<Cm.inverse()<<endl;
    
    homo<<endl<<"Si="<<endl<<endl;
    homo<<Ci.inverse()<<endl;
    
    homo<<endl<<endl<<"--------------------------------------------------------------------"<<endl<<endl;
    
    homo<<endl<<"Homogenized Load Cases (Each row is one load case i.e (E11 E22 E33 E32 E31 E12) and Eij is summation of strain components in in  all voxels )"<<endl<<endl;
    homo<<epHomo<<endl;
    
    homo<<endl<<"Homogenized stess (Each coulemn is (SS11 SS22 SS33 SS32 SS31 SS12) and SSij is summation of stress components in in all voxels )"<<endl<<endl;
    homo<<sigHomo<<endl;

}


void FFT_Analysis2D::printer ()
{

   int ss= mkdir("Results",  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
   cout<<ss;

   ofstream op0("epsilon1.txt");
   for (int i=0; i<xy; ++i)
   op0<<epsilon(0,i).real()<<endl;

   ofstream op1("epsilon2.txt");
   for (int i=0; i<xy; ++i)
   op1<<epsilon(1,i).real()<<endl;

   ofstream op2("epsilon3.txt");
   for (int i=0; i<xy; ++i)
   op2<<epsilon(2,i).real()<<endl;

			//printer surfaces----------------------------------------------
   ofstream op00("sigma1.txt");
   for (int i=0; i<xy; ++i)
   op00<<epsilon(0,i).real()<<endl;

   ofstream op01("sigma2.txt");
   for (int i=0; i<xy; ++i)
   op01<<epsilon(1,i).real()<<endl;

   ofstream op02("sigma3.txt");
   for (int i=0; i<xy; ++i)
   op02<<epsilon(2,i).real()<<endl;


}

FFT_Analysis2D::~FFT_Analysis2D()
{ 
    delete [] t1;
    delete [] t2;
    delete [] GF;
    //delete [] dims;
    // delete [] howmany_dims;
    
    fftw_destroy_plan(p3);
    fftw_destroy_plan(p4);
    fftw_cleanup();
    
    RVE.resize(0,0);   
    epsilon.resize(0,0);
    sigma.resize(0,0);
    tau.resize(0,0);
    sigmaOld.resize(0,0);
		  
}
