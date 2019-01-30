/*

This program is for calculating the kernel, prior and prior precision function


//--------------------------------------------------------------------------
April 22:

Using Guassian Kernel
Using constant prior precision function
Using logarithm of Guassian  as prior precision

*/


//#include<random>
#include<fstream>
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include <Eigen/Dense>
using namespace Eigen;
#include "configuration.h"
using namespace std;

double absolute(double x) {if(x<0)x=-x; return x;}

int compare (const void * a, const void * b)
{
  if ( *(double*)a <  *(double*)b ) return -1;
  if ( *(double*)a == *(double*)b ) return 0;
  if ( *(double*)a >  *(double*)b ) return 1;
}

void swap1(double *a, double *b)
{
	double temp=*a;
	*a=*b;
	*b=temp;
}


//---------------------------------standard Guassian Kernel:
double stdK(double *x1, double *x2, double l)
{
	double sum=0;
	for(int i=0; i<NK; i++)
		sum+=pow(x1[i]-x2[i],2);

	return exp(-sum/l/l/2.0);
}

//------------------------------------------protein-distance Guassian Kernel
double K(double *x1, double *x2, double l)
{
				double sum;

				sum = pow(RMSD_distance(x1, x2),2);
			

				return exp(-sum/l/l/2.0);
}


//--------constant prior precision function
//----------------------------------------------------Ortega's subroutinues------------------
double precision_prior(double *x_new)
{
				return 0.1;
}


//--------------logarithm of Guassian prior:

double prior(double *x_new)
{

	double sum=0;
	for(int i=0; i<6+2*NK; i++)
		 sum+=(x_new[i]-Mu)*(x_new[i]-Mu);

	return -sum/sd_prior/sd_prior/2.0;
}


double kernel0(double *x)
{
			//return 5;
			return 0;
}


double gradient_K(double *xi, double l, double *x, int n)
{
		return (xi[n]-x[n])/l/l*K(xi, x, l);
}

double gradient_prior(double *x, int n)
{
		return -(x[n]-Mu)/sd_prior/sd_prior;
}

double gradient_k0(double *x, int n)
{
	return 0;
}
//----------------------------------------------------Ortega's subroutinues------------------


//-------------------------------------------random generation
//default_random_engine generator1, generator2;
//uniform_real_distribution<double> distribution1(0.0,1.0);
//normal_distribution<double> distribution2(0.0,1.0);

double uniform_sampling_oneD(double lower, double upper)
{
	return ((double)rand())/((double)RAND_MAX)*(upper-lower)+lower;
//	return distribution1(generator1)*(upper-lower)+lower;
}



double normal_sampling(double u, double o)
{
	double s;
	s = (sqrt(-2.0*log(uniform_sampling_oneD(0, 1)))) * cos (2*pi_1*uniform_sampling_oneD(0, 1));
//	s = distribution2(generator2);
	return s*o+u;
} 

double distance(double *x, double *y, int d)
{
	double dis = 0;
	for(int i=0; i<d; i++)
	{
		dis += pow(x[i]-y[i], 2);
	}
	return sqrt(dis);
}



double rho_calc(Samples sample)
{
	double s1;
	double s2;

	s1=s2=0;
	for (int i=0; i<sample.n; i++)
		for(int j=0; j<sample.n; j++)
		{
			double temp = K(sample.x[i], sample.x[j], sample.l);
			if (i==j) s2 +=temp;
			s1 +=temp;
			//cout<<temp<<endl;
		}
	cout<<10*sample.n<<"---- "<<s2<<" "<<s2/s1<<endl;
	if (sample.n==0) 
		return rho_zero * epsilon_1;
	else
		return rho_zero * (epsilon_1 + (double)(sample.n)*s2/s1);

}

double NW_estimator(Samples sample, double *x_new)
{
	double s1;
	double s2;

	s1=s2=0;
		for(int i=0;i<sample.n;i++)
		{
			s1+=K(sample.x[i],x_new, sample.l)*sample.y[i];
			s2+=K(sample.x[i],x_new, sample.l);
		}
	
	s1+=kernel0(x_new)*prior(x_new);
	s2+=kernel0(x_new);

	//cout<<abs1(s1/s2-prior(x_new))<<endl;
	
	return s1/s2;

}

//-------------------------------------------------------------------calculating the rho
double rho_entropy(double (*pos_sample)[NK], int n, Samples sample)
{
	// calculate the entropy:
	fstream out;
	out.open("entropy_sample", ios::out);
	

	double H;
	double p[n];
	double scale = pow(2*pi_1*sample.l*sample.l, -NK/2.0);
	for(int i=0; i<n; i++)
	{
		double s=0;
		for(int j=0; j<n; j++)
			s+= stdK(pos_sample[i], pos_sample[j], sample.l)*scale;

		p[i] = s/n;


	}
	for(int i=0; i<n; i++) out<<pos_sample[i][0]<<" "<<p[i]<<endl; out.close();

	double Z=0;

	for(int i=0; i<n;i++)
		Z += exp(sample.rho*kriging(sample, pos_sample[i])) / p[i];


	Z/=n;

	H=0;
	for(int i=0; i<n; i++)
		H+= -log(exp(sample.rho*kriging(sample, pos_sample[i]))/Z) /log(2.0);
	
	H/=n;

		//cout<<Z<<" "<<H<<"RHO: "<<sample.rho<<endl;
	H=0;
	for(int i=0; i<n; i++)
		H+= -log(p[i]) /log(2.0);
	H/=n;

	cout<<"entropy:"<<H<<endl;
	//for(int nn=sample.n; nn<=500; nn+=10)
	//cout<<H<<" hehe "<<rho_zero*exp(1.0/H*pow(double(nn), 1.0/(dimen+2)))<<endl;
	return rho_zero*exp(1.0/H*pow(double(sample.n), 1.0/(NK/2.0+2.0)));
}


double coeff_cal(Samples *sample)
{
	MatrixXd KK((*sample).n, (*sample).n), y(1, (*sample).n), coeff(1,(*sample).n);
	for (int i=0; i<(*sample).n ;i++)
	{
			for(int j=0; j<(*sample).n ;j++)
			{
				if(i == j)
					KK(i,j) = K((*sample).x[i], (*sample).x[j], (*sample).l) + 0.5;
				else
					KK(i,j) = K((*sample).x[i], (*sample).x[j], (*sample).l);

			}
			y(0,i) = (*sample).y[i];
	}

	coeff = y * KK.inverse();
	//cout<<KK.inverse();
	
	for(int i=0; i<(*sample).n; i++)
	{
		(*sample).coeff[i] = coeff(0, i);
	}

}

double kriging(Samples sample, double *x_new)
{
	double sum=0.0;
	
	for (int i=0; i<sample.n ;i++)
	{//cout<<K(sample.x[i], x_new, sample.l)<<'a'<<endl;
		sum += K(sample.x[i], x_new, sample.l) * sample.coeff[i];
		//cout<<K(sample.x[i], x_new, sample.l) * sample.coeff[i]<<'a'<<" "<<K(sample.x[i], x_new, sample.l)<<" "<<sample.coeff[i]<<endl;;
	}
	

	return sum;
}

//----------------------------------------------convert NK-dimension cartisian coordinates to spherical coordinates
void cart_sph_conver(double *coor, double *phi)
{
	double sum=pow(coor[NK-1], 2);
	
	for(int i=NK-2; i>=0; i--)
	{

		sum+=pow(coor[i],2);

		if(sum<10E-10)
			phi[i]=0;
		else
			phi[i]=acos(coor[i]/sqrt(sum));

	}
	
	if(coor[NK-1]<0)
		phi[NK-2]*=-1;

}

//-----------------------------------------------convert NK-dimension spherical coordinates to cartisian coordinates
void sph_cart_conver(double *coor, double *phi)
{
	coor[0]=1;
	for(int i=1; i<NK; i++)
	{
		coor[i]= coor[i-1]*sin(phi[i-1]);
	}
	for(int i=0; i<NK-1; i++)
		coor[i]*=cos(phi[i]);
}



double scoring(double *feature)
{

	

	double c_bond, c_angl, c_dihe, c_impr, c_urey;
	double c_vdw;
	double c_sasa;
	double c_ele_npo;	
	double scores;

	fstream fin;
	system(current_directory"/minmiz.sh "current_directory" "pdb_crd" "charmm_setup" "charmm_exe" "crd_pdb);	
	fin.open(current_directory"/minmiz/com_energy.dat", ios::in);

	fin>>c_bond>>c_angl>>c_dihe>>c_impr>>c_urey>>c_vdw>>c_sasa>>c_ele_npo;

	fin.close();
	

	feature[0] = c_bond - un_bond;
	feature[1] = c_angl - un_angl;
	feature[2] = c_dihe - un_dihe;
	feature[3] = c_impr - un_impr;
	feature[4] = c_urey - un_urey;
	feature[5] = c_vdw  - un_vdw;
	feature[6] = c_sasa - un_sasa;
	feature[7] = c_ele_npo - un_ele_npo;

	char run[charsize];
	char path[charsize]=scoring_path;
	snprintf(run, sizeof(run), "python scoring_rmsd.py %lf %lf %lf %lf %lf %lf %lf %lf %d %s", feature[0], feature[1],
		feature[2], feature[3], feature[4], feature[5], feature[6], feature[7], nresir, path);
	system(run);

	fin.open("scores", ios::in);
	fin>>scores;
 	fin.close();
 	return scores;
}


double true_RMSD(char *unbound_protein_path)
{
        FILE *Irmsd;

        double true_rmsd;
        char run[charsize];

        snprintf(run,sizeof(run),current_directory"/true_rmsd "current_directory" "protein_name" "cnma_path" "protein_path" %s", unbound_protein_path);
     	system(run);

        Irmsd = fopen(current_directory"/cNMA_output/bound/iRMSD.txt", "r");

        fscanf(Irmsd, "%lf", &true_rmsd);

        fclose(Irmsd);

        return true_rmsd;

}
