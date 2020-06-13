//  This program is sampling the next data point form the posterior distribution 

//Version0.1, 5.17.2018
// rs8 Normal in ball (0,2.5)

#include<fstream>
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <sstream>
#include "configuration.h"
using namespace std;

//------------------------------------------------------------adjust the stepsize according to the acceptance ratio
void stepsize_adaption(double *stepsize, double *acratio, int t2)
{


	if ( (t2+1) % (burning_MCMC/100) == 0)
			{
					
					if ( (*acratio) / double(burning_MCMC/100) > 0.55)
					{
						(*stepsize) = (*stepsize) * 1.1;
					}
					else
					if ((*acratio) / double(burning_MCMC/100) < 0.45)
					{
						(*stepsize) = (*stepsize) * 0.9;
					}
					*acratio = 0.0;

			}
}

//----------------------------------------------------------------The classic Metropolis criteria----------------------
int Metropolis_criteria(double *x_new, Samples sample, double *current, double stepsize)
{

		//	double phi[NK-1];
			//double cand_phi[NK-1];
			double cand[NK];
			
			//current = sample.rho*kriging(sample, xin);

			//cart_sph_conver(x_new, phi);    //---------In MCMC, we use spherical coordinates


				for(int i=0;i<NK;i++)  //-------------------------------------generate proposal
					cand[i]=normal_sampling(x_new[i], stepsize);
	
		//	sph_cart_conver(cand, cand_phi);


			//--------------------------------------------------- calculate acceptance probability
			
				double candidate = sample.rho * kriging(sample, cand);
	         

				double u=uniform_sampling_oneD(0 , 1);


			
				if(judgebound(cand) && u<exp(candidate - *current))
				{
						*current = candidate;
						memcpy(x_new, cand,  sizeof(cand));
						return 1;
				}

			
			return 0;
}


void bayesiandock()
{

	srand(time(NULL));	


	Samples sample;

	double feature[8];
	double parameter;
	double current;
	double ans[NK];
	double x_new[NK];
	double cand[NK-1];
	double cand_phi[NK-1];
	double acratio;
	double acac;
	double temp;
	double ppp1, ppp2;
	double stepsize;
	int n_pos;

	double x_current_best[NK];
	double f_current_best;
	double pos_sample[n_sample*pick_MCMC][NK];
	
	int plotn=0;
	sample.n = 0;
	f_current_best = -oo;

	// -----------------------generate normal random points;
	//FILE *ac=fopen(current_directory"/Result/acratio.dat","w");
	fstream zhongjian;
	zhongjian.open(current_directory"/Result/acratio.dat", ios::out);
	


	for(int t1=0; t1<30; t1++)
	{
		
		sample_coefficient(x_new);
		
		if(t1==0) for(int i=0; i<NK; i++) x_new[i]=0;

		do_sampling(x_new,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
				current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");

		for(int i=0;i<NK;i++)
				{ 
					sample.x[sample.n][i]=x_new[i]; //cout<<sample.x[(sample.n) + t2][i]<<" ";
				}  
			
				sample.y[sample.n]=-scoring(feature)-20;

				if(sample.y[sample.n]>f_current_best)  //-----------------find the current best points
				{
					f_current_best = sample.y[sample.n];
					memcpy(x_current_best,x_new,sizeof(x_new));
				}

			sample.n++;
	}

	sample.l = h0 * pow(sample.n, alfa);
	sample.rho = rho_zero*1;
	coeff_cal(&sample);



//-------------------------------------------------------
//--------------------------------------shuchu
	fstream fin1, fin2;
	fin1.open("x", ios::out);
	fin2.open("y", ios::out);
	for(double i=-2.5; i<=2.5; i+=0.01)
	{ 
		for(int j=0; j<NK; j++) x_new[j]=1;
		x_new[0]=i;
		
		fin1<<i<<endl;
		
		fin2<<kriging(sample, x_new)<<" "<<NW_estimator(sample, x_new)<<endl;;
	}
	fin1.close();
	fin2.close();

	//exit(0);
//---------------------------------------jieshu shuchu
//------------------------------------------------------

	FILE *hehe;

	for(int t1=0;t1<iteration;t1++)
	{


		

		sample_coefficient(x_new);
		


		acac = 0;
		acratio = 0.0;
		stepsize = Metropolis_Steplenth * 4.0 ;///sample.rho ;

		current = sample.rho*kriging(sample, x_new);

//------------------------------------------------------------------burning MCMC-------------------------
		for(int t2=0; t2<burning_MCMC; t2++)
		{	

			int augment = Metropolis_criteria(x_new, sample, &current, stepsize);

			acac+=augment;
			acratio+=augment;

			
			//------------------------------------------tuning the stepsize of MCMC
			stepsize_adaption(&stepsize, &acratio, t2);
	
		}

	//	fprintf(ac, "%lf\n", acac / burning_MCMC);

//--------------------------------------------------------------------------pick MCMC

		n_pos=0;
		acratio = 0.0;
		ppp1 = 0;
		ppp2 = 0;
		for(int t2=0; t2<n_sample; t2++)
		{

			for (int t3=0; t3<pick_MCMC; t3++)
			{
				int augment = Metropolis_criteria(x_new, sample, &current, stepsize);

				//acac+=augment;
				acratio+=augment;


				memcpy(pos_sample[n_pos], x_new, sizeof(x_new));
				n_pos++;

				stepsize_adaption(&stepsize, &acratio, t3+t2*pick_MCMC);

			}



				//--------------------------------------pick a new sample point
				for(int i=0;i<NK;i++)
				{ 
					sample.x[(sample.n)+t2][i]=x_new[i]; //cout<<sample.x[(sample.n) + t2][i]<<" ";
				}  
					
					ppp1 +=x_new[0]; ppp2 +=x_new[1];
					
					
							
					do_sampling(x_new,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
								current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");


					sample.y[(sample.n)+t2]=-scoring(feature)-20;

					if(sample.y[(sample.n)+t2]>f_current_best)
					{
					
						f_current_best = sample.y[(sample.n)+t2];
					
						memcpy(x_current_best,x_new,sizeof(x_new));
					}


		}


		sample.n += n_sample;
		sample.l = h0 * pow(sample.n, alfa);
		sample.rho = rho_entropy(pos_sample, n_pos, sample);        // calculate the constant wrt xx
		coeff_cal(&sample);

		cout<<"mean: "<<x_current_best[0]<<" "<<x_current_best[1]<<" n: "<<sample.n<<endl;;
		cout<<"meansample:"<<ppp1/n_sample<<" "<<ppp2/n_sample<<endl;
		cout<< "ac:"<<acac / burning_MCMC<<endl; 
		cout<< "rho:"<<sample.rho<<"   "<<"l:"<<sample.l<<endl;

		zhongjian<<"mean: "<<x_current_best[0]<<" "<<x_current_best[1]<<" n: "<<sample.n<<endl;;
		zhongjian<<"meansample:"<<ppp1/n_sample<<" "<<ppp2/n_sample<<endl;
		zhongjian<< "ac:"<<acac / burning_MCMC<<endl; 
		zhongjian<< "rho:"<<sample.rho<<"   "<<"l:"<<sample.l<<endl; 
	}

	//fclose(ac);
	zhongjian.close();


	FILE *qual=fopen("Result/quality","w");
	fprintf(qual, "%lf\n", f_current_best);
	for(int i=0; i<NK; i++)
		fprintf(qual, "%lf ", x_current_best[i]);


	fclose(qual);

	qual=fopen("Result/samples","w");

	for(int i=0; i<sample.n; i++)
	{
		for(int j=0; j<NK; j++)
			fprintf(qual, "%8.3f", sample.x[i][j]);

		fprintf(qual, "%10.3f\n", sample.y[i]);
	}

	fclose(qual);
	
	uncertainty(sample, x_current_best, stepsize);
	posterior_anaylsis(x_current_best);
	

	exit(0);

	//	posterior_anaylsis();


	}
/*
	uncertainty(sample, x_current_best);


//----------------------------------------------------------------------End of data points collection and obtain the final posterior 











//---------------------------------------------------------------------------STARTING POSTERIOR ANALYSIS!!---------------------------------------------

	double starting_point[dimen];
	double ansp[dimen];
	double localmaxima[1000][dimen];
	double globalmaxima[dimen];


	int numlocalmaxima,k;
	double maxvalue;
	bool judge;
	
	maxvalue=-2147483640;
	numlocalmaxima=0;


//-------------------------------------------------------------------maxvalue is the global maximum of the posterior
	for(int i=1;i<=LBFGS_num; i++)
	{
		for(int j=0;j<dimen;j++)
			starting_point[j]=uniform_sampling(LOWERBOUND,UPPERBOUND);
		temp=LBFGS(starting_point,sample, ansp);

		if(temp>maxvalue)
			{
				maxvalue=temp;
				memcpy(globalmaxima, ansp, sizeof(ansp));
			}

		int j;
		for(j=1;j<=numlocalmaxima;j++)
		{
			judge=0;
			for(k=0;k<dimen;k++)
				if(!equalrule(ansp[k],localmaxima[j][k]))
				{     
					judge=1;
					break;
				}
			if(judge==0) break;
		}
		if(j==numlocalmaxima+1)
		{
			numlocalmaxima++;
			for(k=0;k<dimen;k++)
				localmaxima[numlocalmaxima][k]=ansp[k];
		}

	}


//----------------------------------------------------------find the best samples

	double best_sample[dimen];
	double best = -214748367;

	for (int i=0; i<sample.n; i++)
	{
		if (best < sample.y[i]){
			best =sample.y[i];
			memcpy(best_sample, sample.x[i], sizeof(sample.x[i]));
		}
	}

//	uncertainty (sample, best_sample, maxvalue, "best_smp");
//	uncertainty (sample, globalmaxima, maxvalue, "gmax");

}

*/
	double RMSD[rmsd_sample];
	double pmi[sample_pmi][NK];
	double numerical[sample_pmi];

void uncertainty(Samples sample, double *x_center, double stepsize)
{


	double current;
	double x_new[NK];
	double acratio;
	double acac;

//----------------------------------------------------------------------do Monte Carlo integration

	
	for(int i=0; i<NK; i++)x_new[i]=x_center[i];	

//		do_sampling(x_new,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
//				current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");

		current = sample.rho*kriging(sample, x_new);

		acac = 0;
		acratio = 0.0;

//------------------------------------------------------------------burning MCMC-------------------------
		for(int t2=0; t2<burning_MCMC; t2++)
		{	

			int augment = Metropolis_criteria(x_new, sample, &current, stepsize);

			acac+=augment;
			acratio+=augment;

			
			//------------------------------------------tuning the stepsize of MCMC
			stepsize_adaption(&stepsize, &acratio, t2);

			
	
		}
//---------------------------------------------------------------------end burning mcmc
	cout<< "ac:"<<acac / burning_MCMC<<endl; 
//----------------------------------------------------------------------------calculate gbest
	

//--------------------------------------------------------------------------pick MCMC

		acratio=0;

		for(int t2=0; t2<rmsd_sample; t2++)
		{

			int augment = Metropolis_criteria(x_new, sample, &current, stepsize);

			
			acratio+=augment;

			RMSD[t2] = RMSD_distance(x_new, x_center);
			if(t2<sample_pmi)
				memcpy(pmi[t2], x_new, sizeof(x_new));

			stepsize_adaption(&stepsize, &acratio, t2);

			
		}


//----------------------------------------------------------------------calculating Pmi---------------------------
	double Z=-oo;	
	double P_mi=0;

	double scale = pow(2*pi_1*sample.l*sample.l, -NK/2.0);

	for(int i=0; i<sample_pmi; i++)
	{
		double s=0;
		for(int j=0; j<sample_pmi; j++)
			s+= stdK(pmi[i], pmi[j], sample.l)*scale;


		double probability = s/sample_pmi;
		numerical[i] = sample.rho*kriging(sample, pmi[i]) - log(probability);  //here we do a trick to prevent the numerical issue.

		if(numerical[i] > Z) Z=numerical[i];
	}

	for(int i=0; i<sample_pmi; i++)
		P_mi += exp(numerical[i]-Z);

	P_mi = log(P_mi/sample_pmi) + Z;
//----------------------------------------------------------------------------------------





	FILE *rmsd;
	FILE *a2;
	FILE *a3;
	FILE *numerica;

	qsort(RMSD, rmsd_sample, sizeof(double), compare);


	numerica=fopen(current_directory"/Result/numerical.dat","w");

	for(int i=0; i<sample_pmi; i++)
		fprintf(numerica, "%.4lf\n", numerical[i]);
	fprintf(numerica, "%.4lf\n", Z);


	rmsd=fopen(current_directory"/Result/Rmsd_dis.dat","w");

	for(int i=0; i<rmsd_sample; i++)
		fprintf(rmsd, "%.4lf\n", RMSD[i]);

		
	a2 = fopen(current_directory"/Result/uncertainty.dat","w");

	int i;
	for(i=0; i<rmsd_sample; i++)
		if(RMSD[i]>4)
			break;

	fprintf(a2, "%10.4lf\n", (double(i))/double(rmsd_sample));
	

	a3 = fopen(current_directory"/Result/PMI_log","w");  //----In order to keep numerical stable, we record log(P(Mi))

	fprintf(a3, "%10.4lf\n", P_mi);


	fclose(numerica);
	fclose(a3);
	fclose(a2);
	fclose(rmsd);
}


