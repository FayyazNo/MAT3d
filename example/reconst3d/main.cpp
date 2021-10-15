#include "reconstruction.hpp"
using namespace recon;
using namespace std;

int main ()
{
    /*
    // Recnstuction base on the two cut section and estimate of the TPCF
    CImg<unsigned char> C1 ("5.tiff");
    CImg<unsigned char> C2;
    C2=C1;
    recPlan a (C1, C2);               // Set reconstraction plan. The arguments are
                                      // two perpendicular cut sections of the sample
                                      
    a.estimateNormFFT_RVE();          // Estimate norm of teh Fourier transform of the RVE
                                      // based on the 2 cut section c1 and c2
    a.runRec(300, .8);                // Run reconstruction the first argument is number of the
                                      // allowable iteration and seconed one is beta
                                      // where used in phase recovery algorithm and 0<b<1
   
    a.saveRecImages();                // Export output results.*/
    
    
    //**************************************************************************
    // Reconstruction directly based on the TPCF
    

  	 		



   /*
   The scaled autocovariance function define as:
   
        f(r)=(S_2 (r)-fi_1^2)/(fi_1*fi_2)    (1)
   
   where S_2(r) id TPCF of phase 1 and fi_1 and fi_2 are valume fractio of the 
   phase 1 and 2 respectively and fi_1+fi_2=1
   
   Diffrent definition is presented for f(r) based on the Torquato (Necessary Conditions on Realizable Two-Point
   Correlation Functions of Random Media,arXive 2014) on(e?) form of f(r) for statistically isotropic
   heterogeneos material is:
   
       f(r)=exp(-r/A)* sin (Q*r)/(Q*r)      (2)
   
   where a nd q are constants and r is vetor length (here distance from center of rve)
   
   Combining (1) and (2) gives S_2(r) which can be used in reconstruction.
   */
	int Nx=171, Ny=171, Nz=171;
	int xyz=Nx*Ny*Nz;

	int i_mid,   j_mid,   k_mid;
	i_mid= Nx/2; j_mid= Ny/2;  k_mid= Nz/2;

	fftw_complex* mTilda = new fftw_complex[Nx*Ny*Nz];
	fftw_complex *MTilda = new fftw_complex[Nx*Ny*Nz];
	fftw_plan p1 = fftw_plan_dft_3d(Nx, Ny, Nz, mTilda, MTilda, FFTW_FORWARD, FFTW_ESTIMATE); //fft3 plan 
	double *tpcf=new double [Nx*Ny*Nz];
	double *TPCF=new double [Nx*Ny*Nz];

	double A=25;         // constant in eq. (2) (see above note)
	double Q=.4;        // constant in eq. (2) (see above note)
	double fi_1=.50;      // volume fraction phase 1
	double fi_2=1-fi_1;  // volume fraction phase 2
	double fr=0;
	double c1, c2;
	double r1=0, r2=0;
	int _i=0;

	for (int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{

				r1=sqrt(                           (i-i_mid)*(i-i_mid)+\
		                       1/((1.+_i)*(1.+_i))*(j-j_mid)*(j-j_mid)+\
						       1/((1.+_i)*(1.+_i))*(k-k_mid)*(k-k_mid));

				r2=sqrt(                (i-i_mid)*(i-i_mid)+\
									    (j-j_mid)*(j-j_mid)+\
					          1/((1*1))*(k-k_mid)*(k-k_mid));

				c1=0.3; c2=.7;
				fr= exp(-r1/A)*sin(Q*r1)/(Q*r1);
				//fr= (1+pow(1.3,3)*.006)/(1+.006*pow((r1/1.3+1.3),3)); 
				//fr=exp(-pow((r1/8.+.13),2))/exp(-pow(.13,2));

				if (!r1 )
					fr=1;

				tpcf[k+(Nz)*j+(Nz*Ny)*i]= (fr*fi_1*fi_2)+(fi_1*fi_1);
			}


	recon::ifftshift(tpcf, TPCF, Nx,  Ny, Nz);  



	for (int i=0; i<xyz; ++i)
		mTilda[i][0]=TPCF[i];

	fftw_execute(p1);

	for (int i=0; i<xyz; ++i)
		TPCF[i]=MTilda[i][0];

   //------------------------------
	recon::recPlan aa(TPCF, Nx, Ny, Nz);     // Constructor which take TPCF and RVE dimensions
	aa.runRec(100, .8);                       // Run reconstruction. the first argument is number of the
		                                    // allowable iteration and seconed one is beta
				                           // where used in phase recovery algorithm and 0<b<1

	aa.saveRecImages("ar"+to_string((_i+1)));   
	aa.vtiWriter("ar"+to_string((_i+1))+"/rec.vti", tpcf);

	if (TPCF)   delete [] TPCF;
	if (tpcf)   delete [] tpcf;
	if (mTilda)	delete [] mTilda;
	if (MTilda)	delete [] MTilda;    
    fftw_destroy_plan(p1);

    return 0;
}
