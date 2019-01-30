//  this program is to rotate the ligand at the center of interface of ligand

#include<fstream>
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <sstream>
#include "configuration.h"
using namespace std;


//--------------------------------------------------------------This program is for rigid-body sampling
double R_L_distance;
double rotation[3][3];
double center_rx, center_ry, center_rz;
double center_lx, center_ly, center_lz;
double center_int_ligx, center_int_ligy,center_int_ligz;
double rmsd_receptor_prediction;
//----------------------------------------------------only for interface heavy atoms, normal modes as same, xyz for ligand-ref
double sum_r_r[NK][NK];
double sum_int_c_c[NK][NK];
double sum_wt_l[12][NK];
double sum_x;
double sum_y;
double sum_z;
double sum_xx;
double sum_xy;
double sum_xz;
double sum_yy;
double sum_yz;
double sum_zz;
double sum_t1;
double sum_t2;
double sum_t3;

//double scaling_factor;

/*
void rigid_coor_transfer(double *coor, double *t, double *w)
{
	double stheta;
	double sphi;
	double lmx, lmy, lmz;
	double u1, u2, u3;
	double sumofu;
	double sumofb;

//-----------------------------------------------------------------translation
	stheta = atan2(-coor[1] , coor[2]);
	sphi =2*(asin(3.0/R_L_distance))*sqrt(coor[1]*coor[1]+coor[2]*coor[2]);   

	lmx = R_L_distance*cos(stheta)*sin(sphi);
    lmy = R_L_distance*sin(stheta)*sin(sphi);
    lmz = R_L_distance*cos(sphi);



    t[1]= rotation[0][0]*lmx + rotation[0][1]*lmy + rotation[0][2]*lmz + center_rx;
    t[2]= rotation[1][0]*lmx + rotation[1][1]*lmy + rotation[1][2]*lmz + center_ry;
    t[3]= rotation[2][0]*lmx + rotation[2][1]*lmy + rotation[2][2]*lmz + center_rz;
    


    u1 = t[1] - center_rx;
    u2 = t[2] - center_ry;
    u3 = t[3] - center_rz;
    sumofu = sqrt(u1*u1 + u2*u2 + u3*u3);
    u1 /= sumofu;
    u2 /= sumofu;
    u3 /= sumofu; 

    t[1] = t[1] + u1*coor[0];
    t[2] = t[2] + u2*coor[0];
    t[3] = t[3] + u3*coor[0];

    if(rigid_rotate==1){                            //--rotate around the center of interface of ligand
  	  	t[1] = t[1] - center_lx + center_int_ligx;
  	 	 t[2] = t[2] - center_ly + center_int_ligy;
  		 t[3] = t[3] - center_lz + center_int_ligz;
    					}
   //----------------------------------------------rotation------------------------

    double b0, b1, b2;
    double b4;

    sumofb = sqrt(coor[3]*coor[3] + coor[4]*coor[4] + coor[5]*coor[5]);
    if( sumofb ==0 ){ b0 = b1 = 0; b2 = 1;}
    else
    	{
    	b0 = coor[3]/sumofb;
    	b1 = coor[4]/sumofb;
    	b2 = coor[5]/sumofb;
    	}

    w[1] = (b0)*(b0)*(1-cos(sumofb))+cos(sumofb);
    w[2] = (b0)*(b1)*(1-cos(sumofb))-(b2)*sin(sumofb);
    w[3] = (b0)*(b2)*(1-cos(sumofb))+(b1)*sin(sumofb);
    w[4] = (b0)*(b1)*(1-cos(sumofb))+(b2)*sin(sumofb);
    w[5] = (b1)*(b1)*(1-cos(sumofb))+cos(sumofb);
    w[6] = (b1)*(b2)*(1-cos(sumofb))-(b0)*sin(sumofb);
    w[7] = (b0)*(b2)*(1-cos(sumofb))-(b1)*sin(sumofb);
    w[8] = (b1)*(b2)*(1-cos(sumofb))+(b0)*sin(sumofb);
    w[9] = (b2)*(b2)*(1-cos(sumofb))+cos(sumofb);


} 
*/
//---------------------This sampling is for complex_normal_modes_sampling
void do_sampling(double *coor_1, char *lpath, char *result_l, char *rpath, char *result_r)
{
	

	FILE *post_ligand;
	FILE *post_receptor;
	FILE *mean;

	double x,y,z;
	double scaling_factor;
	double displacement;
	double mean_dis;
	double coor[NK];

	post_ligand = fopen(result_l, "w");
	post_receptor = fopen(result_r, "w");
	mean= fopen(output_path"/Result/mean_dis","a");

	mean_dis = 0;

	for(int i=0; i<NK; i++) mean_dis+=pow(coor_1[i], 2); // cout<<sqrt(mean_dis)*rmsd_receptor_prediction<<endl; 

	if(mean_dis==0)
		scaling_factor=0;
	else
		scaling_factor = rmsd_receptor_prediction*sqrt(nresir/l2_norm(coor_1)*mean_dis) ;

	for(int i=0; i<NK; i++) coor[i] = coor_1[i]*scaling_factor/sqrt(eigenvalues[i]);
//		cout<<scaling_factor<<endl;

	//for(int i=0; i<NK; i++) cout<<coor_1[i]<<" "; cout<<endl;
//	for(int i=0; i<NK; i++) cout<<coor[i]<<" ";  cout<<endl;

//------------------------------------------start receptor----------------------------------------

	mean_dis = 0;

	for(int i=0; i<natomr; i++)
	{

		x = coor_rec[i][0];
		y = coor_rec[i][1];
		z = coor_rec[i][2];

 		displacement=0;                                              //---------------displacement for x
 		for(int j=0;j<NK; j++)
 			displacement += eigenvectors[j][3*residue_ind_rec[i]]*coor[j];

 		
 		x=x+displacement;
 		
 		if(atom_name_rec[i][0]=='C' && atom_name_rec[i][1]=='A')
 			 mean_dis +=displacement*displacement; //cout<<residue_ind_rec[i]<<endl;}

 		displacement=0;                                              //---------------displacement for y
 		for(int j=0;j<NK; j++)
 			displacement += eigenvectors[j][3*residue_ind_rec[i]+1]*coor[j];

 		
 		y=y+displacement;
 		
 		if(atom_name_rec[i][0]=='C' && atom_name_rec[i][1]=='A')
 			mean_dis +=displacement*displacement;

		displacement=0;                                              //---------------displacement for z
 		for(int j=0;j<NK; j++)
 			displacement += eigenvectors[j][3*residue_ind_rec[i]+2]*coor[j];


 		z=z+displacement;

 		if(atom_name_rec[i][0]=='C' && atom_name_rec[i][1]=='A')
 			mean_dis +=displacement*displacement;



		fprintf(post_receptor, "ATOM%7d %-5s%s %c%4d", i+1, atom_name_rec[i], residue_name_rec[i], chainid_rec[i], residue_ind_rec[i]+1);
 		fprintf(post_receptor, "%12.3lf%8.3lf%8.3lf%6s%6s\n", x, y, z, charge1_rec[i], charge2_rec[i]);

 		sample_rec[i][0] = x;
 		sample_rec[i][1] = y;
 		sample_rec[i][2] = z;
	}

	cout<< "mean_dis_move_rec:"<<sqrt(mean_dis/double(nresir))<<endl;
	fprintf(mean, "mean_dis_move_rec%.4lf\n", sqrt(mean_dis/double(nresir)));

	mean_dis = 0;
//-------------------------------------------------------start_ligand
	for(int i=0; i<natoml; i++)
	{

		x = coor_lig[i][0];
		y = coor_lig[i][1];
		z = coor_lig[i][2];


 		displacement=0;                                              //---------------displacement for x
 		for(int j=0;j<NK; j++)
 			displacement += eigenvectors[j][3*(residue_ind_lig[i]+nresir)]*coor[j];

 		x=x+displacement;

 		if(atom_name_lig[i][0]=='C' && atom_name_lig[i][1]=='A')
 			mean_dis +=displacement*displacement;

 		displacement=0;                                              //---------------displacement for y
 		for(int j=0;j<NK; j++)
 			displacement += eigenvectors[j][3*(residue_ind_lig[i]+nresir)+1]*coor[j];

 		y=y+displacement;

 		if(atom_name_lig[i][0]=='C' && atom_name_lig[i][1]=='A')
 			mean_dis +=displacement*displacement;

		displacement=0;                                              //---------------displacement for z
 		for(int j=0;j<NK; j++)
 			displacement += eigenvectors[j][3*(residue_ind_lig[i]+nresir)+2]*coor[j];

 		z=z+displacement;

 		if(atom_name_lig[i][0]=='C' && atom_name_lig[i][1]=='A')
 			mean_dis +=displacement*displacement;

		fprintf(post_ligand, "ATOM%7d %-5s%s %c%4d", i+natomr+1, atom_name_lig[i], residue_name_lig[i], chainid_lig[i], residue_ind_lig[i]+1);
 		fprintf(post_ligand, "%12.3lf%8.3lf%8.3lf%6s%6s\n", x, y, z, charge1_lig[i], charge2_lig[i]);

 	//	if(i==0) cout<<x<<" "<<y<<" "<<z<<endl;
 		sample_lig[i][0] = x;
 		sample_lig[i][1] = y;
 		sample_lig[i][2] = z;

	}
	cout<< "mean_dis_move_lig:"<<sqrt(mean_dis/double(nresil))<<endl;
	fprintf(mean, "mean_dis_move_lig%.4lf\n", sqrt(mean_dis/double(nresil)));
//	fstream qop;
//	qop.open("qop", ios::out | ios::app);
//	qop<<sqrt(mean_dis/double(nresil))<<endl;
//	qop.close();


	fclose(post_ligand);
	fclose(post_receptor);
	fclose(mean);
	
}


double l2_norm(double *y)
{
	double l2norm=0;
	for(int i=0; i<NK; i++)
	{
		for(int j=i; j<NK; j++) l2norm += y[i]*sum_r_r[i][j]*y[j];
	}
	return l2norm;
}

//-----------------------------------------sample y within [0,2.5]
void sample_coefficient(double *y)
{
	
		
	double l2norm;
	double radius;
	int in=0;
	double u=oo;
	double ratio;

	while(1) {

	uniform_sampling_multiD(y);		//--------------------------------------y is the unit vector at beginning


	do{
		ratio = normal_sampling(1.27, 0.46);
		}
	while(ratio>2.5 || ratio<0);    //--------------make sure the ratio is between [0,2.5], it is trucated Gaussian.


	radius=0.0;

	for(int i=0; i<NK; i++)
	{
		radius += y[i]*y[i]/eigenvalues[i];
	}

	l2norm = l2_norm(y);


	double thereshold = nresil * pow(range_ligand,2) / pow(rmsd_receptor_prediction*ratio,2) / nresir +1;

	if(radius <= l2norm*thereshold)
		break;
	in++;

	if(in>=100000)
		{
		rmsd_receptor_prediction/=2;
		break;
		}
	}

	fstream ac; 
	ac.open(output_path"/Result/step_uniform_sample", ios::out | ios::app);
	ac<<in<<endl;
	ac.close();

	for(int i=0; i<NK; i++)             //--------------------------------------scale y to become the non-unit vector
		y[i]*=ratio;

}



void uniform_sampling_multiD(double *temp)
{


//-----------------------------------------------sample on the shell of NK-ksional sphere for receptor

	double	s = 0;

	for(int i=0; i < NK; i++)
	{
		temp[i]=normal_sampling(0, 1);
		s+=temp[i]*temp[i];
	}

	//double u=pow(uniform_sampling_oneD(0,1), 1.0/NK);

	for(int i=0; i< NK; i++) 
		temp[i]=temp[i]/sqrt(s);


}


bool judgebound(double *y)
{

	double thereshold;
	double l2norm;
	double radius;
	double ratio=0;

	for(int i=0; i<NK; i++)
		ratio += y[i]*y[i];

	ratio=sqrt(ratio);

	if(ratio==0) return true;
	if(ratio>2.5) return false; //---------------------out of boundary

	thereshold = nresil * pow(range_ligand,2) / pow(rmsd_receptor_prediction*ratio,2) / nresir +1;
	
	
	radius=0.0;

	for(int i=0; i<NK; i++)
	{
		radius += y[i]*y[i]/eigenvalues[i];
	}

	l2norm = l2_norm(y);

	if(radius <= l2norm*thereshold)
		return true;

	return false;
	
}

//------------------------------------------fast calculate the RMSD of putative interface of two sets of coordinates
double RMSD_distance (double *coor1, double *coor2)
{

	double ratio1=0.0, ratio2=0.0, scaling_factor1, scaling_factor2;

	for(int i=0; i<NK ; i++)
	{
		ratio1+=pow(coor1[i],2);
		ratio2+=pow(coor2[i],2);
	}

	if(ratio1==0) scaling_factor1=0;
	else scaling_factor1 = rmsd_receptor_prediction*sqrt(nresir/l2_norm(coor1)*ratio1);
	
	if(ratio2==0) scaling_factor2=0;
	else scaling_factor2 = rmsd_receptor_prediction*sqrt(nresir/l2_norm(coor2)*ratio2);

	double coor[NK];
	for(int i=0; i<NK; i++)
	{
		coor[i] = (coor1[i]*scaling_factor1 - coor2[i]*scaling_factor2)/sqrt(eigenvalues[i]);
	}

	double rmsd=0;
	for(int i=0; i<NK; i++)
	for(int j=0; j<NK; j++)
		rmsd += coor[i] * coor[j]*sum_int_c_c[i][j];

	return sqrt(rmsd/double(int_nresidue));
}
/*
	for(int i=1; i<=3; i++)
		t[i] = t1[i] - t2[i];
	for(int i=1; i<=9; i++)
		w[i] = w1[i] - w2[i];

	for(int i=0; i<NK; i++)
	{
		r[i] = coor1[i+6] - coor2[i+6];
		l[i] = coor1[i+6+NK] - coor2[i+6+NK];
	}

	double rmsd = 0;

	//cout<<t[1]<<" "<<t[2]<<"-"<<t[3]<<endl;
	//for(int i=1; i<=9; i++)
	//	cout<<w[i]<<" "; cout<<endl;
//	cout<<sum_t1<<" "<<sum_t2<<" "<<sum_t3<<endl;
	rmsd +=   (w[1]*w[1] + w[4]*w[4] + w[7]*w[7]) * sum_xx;
	rmsd += 2*(w[1]*w[2] + w[4]*w[5] + w[7]*w[8]) * sum_xy;
	rmsd += 2*(w[1]*w[3] + w[4]*w[6] + w[7]*w[9]) * sum_xz;
	rmsd +=   (w[2]*w[2] + w[5]*w[5] + w[8]*w[8]) * sum_yy;
	rmsd += 2*(w[2]*w[3] + w[5]*w[6] + w[8]*w[9]) * sum_yz;
	rmsd +=   (w[3]*w[3] + w[6]*w[6] + w[9]*w[9]) * sum_zz;
	rmsd += 2*(t[1]*w[1] + t[2]*w[4] + t[3]*w[7]) * sum_x;
	rmsd += 2*(t[1]*w[2] + t[2]*w[5] + t[3]*w[8]) * sum_y;
	rmsd += 2*(t[1]*w[3] + t[2]*w[6] + t[3]*w[9]) * sum_z;
	rmsd += t[1]*t[1] * sum_t1;
	rmsd += t[2]*t[2] * sum_t2;
	rmsd += t[3]*t[3] * sum_t3;
//cout<<sum_l_l[0][0]<<endl;	
//---------------------------------ligand flexible itself
	for(int i=0; i<NK; i++)
	for(int j=0; j<NK; j++)
		rmsd += l[i]*l[j] * sum_l_l[i][j];
//cout<<l[1]<<" "<<l[2]<<" "<<l[0]<<rmsd<<endl;
//---------------------------------receptor flexible itself
	for(int i=0; i<NK; i++)
	for(int j=0; j<NK; j++)
		rmsd += r[i]*r[j] * sum_r_r[i][j];
//cout<<rmsd<<"asas"<<endl;
	for(int i=0; i<12; i++)
		for(int j=0; j<NK; j++)
		{
			if ( i < 9)
			rmsd += 2 * w[i+1] * l[j] * sum_wt_l[i][j];
			else
			rmsd += 2 * t[i-8] * l[j] * sum_wt_l[i][j];
		}
  //cout<<rmsd<<"asas"<<endl;   
	return sqrt(rmsd/double(int_natom));

}
*/

//-------------------------------------------------------------------------this is the feature-------------------

