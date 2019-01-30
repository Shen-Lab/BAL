//---------------------------attention!:  residue_ind starts from index 0 everywhere.
//-------------------------------------This is for the  prepocessing of complex_normal_modes;
#include<fstream>
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<sstream>
//#include<random>
#include "configuration.h"
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

int nresir;    //----------------------num of residues of receptor
int nresil;	  //-----------------------num of residues of ligand
int natomr;	  //-----------------------num of atoms of receptor
int natoml;	  //-----------------------num of atoms of ligand
int int_nresidue; //----------------------num of atoms of the interface (heavy)
//double nomo[NK][3*residue_max];  //---------rec normal modes
//double nomo_lig[NK][3*residue_max];  //---------lig normal modes

bool ptin_rec[residue_max];             //----------putative interface heavy atoms for receptor
bool ptin_lig[residue_max];				//----------putative interface heavy atoms for ligand
double coor_lig[atom_max][3];       //---------ligand atom coordinates
double coor_rec[atom_max][3];		//--------receptor atom coordinates
char charge1_lig[atom_max][6];
char charge1_rec[atom_max][6];
char charge2_lig[atom_max][6];
char charge2_rec[atom_max][6];
int residue_ind_lig[atom_max];  
int residue_ind_rec[atom_max];
char chainid_lig[atom_max];
char chainid_rec[atom_max];
char atom_name_lig[atom_max][5];
char atom_name_rec[atom_max][5];
char residue_name_lig[atom_max][5]; 
char residue_name_rec[atom_max][5];
double sample_lig[atom_max][3];
double sample_rec[atom_max][3];
double Kd;

double eigenvalues[numofnodes];
double eigenvectors[numofnodes][residue_max];
//double radius_soid[NK];
//double transform1[NK][NK];
///double ellipsoid_radius[NK];

//--------------------------------------------unbound energy-------------------------------------------
double un_bond, un_angl, un_dihe, un_impr, un_urey, un_vdw, un_ele_npo, un_sasa;        

void preprocess()
{


	double scaling_factor_r;

	
	char trash1;
	char trash2[charsize];
	
	FILE *rmsd_prediction;
	FILE *Eigenvalue;
	FILE *Kdd;
	string str;


	system("./extract_Kd "protein_name);
    Kdd=fopen(current_directory"/Kd_value","r");
    fscanf(Kdd,"%lf",&Kd);

    fclose(Kdd);

//----------------------------------------------------------------Hbuild by CHARMM-----------------------------------
//
//
//

   	system(current_directory"/charmm_run "current_directory" "protein_name" "protein_path" "pdb_crd" "charmm_setup" "charmm_exe" "crd_pdb);

//------------------------------------------------------------------calculate unbound energy--------------------------
   	double t1,t2,t3,t4,t5,t6,t7,t8;
   	fstream fin;
   	fin.open(current_directory"/prepocess/rec_energy.dat", ios::in);
   	fin>>un_bond>>un_angl>>un_dihe>>un_impr>>un_urey>>un_vdw>>un_sasa>>un_ele_npo;
   	fin.close();

   	fin.open(current_directory"/prepocess/lig_energy.dat", ios::in);
   	fin>>t1>>t2>>t3>>t4>>t5>>t6>>t7>>t8;
   	fin.close();

   	un_bond += t1;
   	un_angl += t2;
   	un_dihe += t3;
   	un_impr += t4;
   	un_urey += t5;
   	un_vdw  += t6;
   	un_sasa += t7;
   	un_ele_npo += t8;

   	cout<<"unbound_energy:"<< un_bond+un_angl+un_dihe+un_impr+un_urey<<" "<<un_vdw<<" "<<un_sasa<<" "<<un_ele_npo<<endl;
//------------------------------------------------------------------end calculate unbound energy--------------------------


//--------------------------------------------obtain num of atoms-----------------------------------
   	
   	fin.open(current_directory"/prepocess/receptor.crd", ios::in);
   	for(int i=0; i<3; i++)
   		getline(fin, str);
   	fin>>natomr;
   	fin.close();

   	fin.open(current_directory"/prepocess/ligand.crd", ios::in);
   	for(int i=0; i<3; i++)
   		getline(fin, str);
   	fin>>natoml;
   	fin.close();



 //---------------------------------------------------------save the pdb information from "prepocess/ligand(receptor).pdb "

   	FILE *receptor, *ligand;

   	char atom[4];
   	int atom_ind;

   	receptor = fopen(current_directory"/prepocess/receptor.pdb", "r");
	ligand   = fopen(current_directory"/prepocess/ligand.pdb", "r");

	for(int i=0; i<natoml; i++)
	{
		fscanf(ligand, "%s %d %s %s %c %d", atom, &atom_ind, atom_name_lig[i], residue_name_lig[i], &chainid_lig[i], &residue_ind_lig[i]);
		fscanf(ligand, "%lf %lf %lf %s %s\n", &coor_lig[i][0], &coor_lig[i][1], &coor_lig[i][2], charge1_lig[i], charge2_lig[i]);
		residue_ind_lig[i]--;
	
	}
	
	for(int i=0; i<natomr; i++)
	{
		fscanf(receptor, "%s %d %s %s %c %d", atom, &atom_ind, atom_name_rec[i], residue_name_rec[i], &chainid_rec[i], &residue_ind_rec[i]);
		fscanf(receptor, "%lf %lf %lf %s %s\n", &coor_rec[i][0], &coor_rec[i][1], &coor_rec[i][2], charge1_rec[i], charge2_rec[i]);
		residue_ind_rec[i]--;

	}

	fclose(receptor);
	fclose(ligand);



//-----------------------------------------------------------------------------------------------

	nresir = residue_ind_rec[natomr-1]+1;
	nresil = residue_ind_lig[natoml-1]+1;

//
//
//
//
//
//




//-------------------------------------------------------calculate normal modes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


	system(current_directory"/runcnma_com "current_directory" "protein_name" "cnma_path);

	rmsd_prediction = fopen(current_directory"/cNMA_output/Rec/Output/RMSDPrediction.txt", "r");
	Eigenvalue = fopen(current_directory"/cNMA_output/Rec/Output/eigenvaluesComplex.txt", "r");   //The complex eigenvalues
	fin.open(current_directory"/cNMA_output/Rec/Output/eigenvectorsComplex.txt", ios::in);        //The complex eigenvectors


//-------------------------------------------------------------calculate receptor normal modes


	
	fscanf(rmsd_prediction, "%c %s", &trash1, trash2);
	fscanf(rmsd_prediction, "%lf", &rmsd_receptor_prediction);
	

	
	fscanf(Eigenvalue, "%c %s", &trash1, trash2);


	//scaling_factor_r=rmsd_receptor*rmsd_up;  //-------------------------------------------------------------Here is the key on the shell
	//*rmsd_up;

	double rescale[numofnodes], index[numofnodes];

//----------------------------------------------sampling through receptor---------------------------------------------------------


     getline(fin, str);
     for(int i=0; i<numofnodes; i++)
     {
     	
     	if(i<6)
     		fscanf(Eigenvalue, "%lf", &eigenvalues[i]);
     	else
     		fscanf(Eigenvalue, "%lf", &eigenvalues[i-6]);


     	getline(fin, str);
     	istringstream ss(str);
     	double value;
    	int n=0;
     	double sum_norm=0;
     	while(ss >> value)
     	{
     			if(i>=6)
     			eigenvectors[i-6][n] = value;

    		n++;
    		if(n<=3*nresir)
    			sum_norm+=value*value;
     	}

     	if(i>=6)
     		{
     			rescale[i-6]= eigenvalues[i-6]/sum_norm;
     			index[i-6]=i-6;
     			cout<<"receptor contributions: "<<sum_norm<<endl;
     	    }

     }

     for(int i=3+NK/2; i<numofnodes-6; i++)
     	for(int j=i; j<numofnodes-6; j++)
     		if(rescale[i]>rescale[j])
     		{
     			swap1(&rescale[i], &rescale[j]);
     			swap1(&eigenvalues[i], &eigenvalues[j]);
     			swap1(&index[i], &index[j]);

     			double candidate[residue_max];
     			memcpy(candidate, eigenvectors[i], 		 sizeof(candidate));
     			memcpy(eigenvectors[i], eigenvectors[j], sizeof(candidate));
     			memcpy(eigenvectors[j], candidate,       sizeof(candidate));
     		}



     //	cout<<"receptor contributions: "<<sum_norm<<endl;
//------------------------------------------------------------------calculate rescaling factor  	

	fin.close();
	fclose(rmsd_prediction);
   	fclose(Eigenvalue);

   	cout<<"normal modes chosen as basis and rescale:"<<endl;
   	for(int i=0; i<NK; i++)
   		cout<<index[i]<<" "<<rescale[i]<<endl;
     

   	cout<<"#residues of receptor:"<<nresir<<endl;
   	cout<<"#residues of ligand:"<<nresil<<endl;
   	cout<<"#atoms of receptor:"<<natomr<<endl;
   	cout<<"#atoms of ligand:"<<natoml<<endl;

  
}


void constant_setting()
{

	double sumr;
	double x[NK];

 	MatrixXd Ra(NK, NK);
//--------------------------------------------------------------calculate interface of ligand-------------------------

	for(int i=0; i<NK; i++)
		for(int j=0; j<NK; j++)
		{
			sumr=0;
			for(int i1=0; i1<nresir*3; i1++)
				sumr += eigenvectors[i][i1]*eigenvectors[j][i1]/sqrt(eigenvalues[i] * eigenvalues[j]);

			Ra(i,j)=sumr;

			if(i!=j) sumr*=2;
			sum_r_r[i][j] =sumr;
						

		}

//-------------------------------------------------------Test the feasibility of the search space
/*
	SelfAdjointEigenSolver <MatrixXd> es(Ra);
	double vol=1;
	for(int i=NK-1; i>=0; i--)
	{
		double a= 1/es.eigenvalues().col(0)[i];
		cout<< range_ligand*sqrt(nresil/((a-1.0)*nresir))<<" "<<rmsd_receptor_prediction<<endl;

		vol *= sqrt(a);
	}
	cout<<pow(pow(range_ligand,2)*nresil + pow(rmsd_receptor_prediction,2)*nresir, NK/2)


	rmsd_receptor_prediction/=2;
*/

//------------------------------------------------done with constant for sampling
	choo_puta_inte_atom();//--------------------------choosing interface atoms
//--------------------------------------------------start setting constants for rapid RMSD calculation	

	for(int i=0; i<nresir; i++)
	{

		if(ptin_rec[i] == 1)
		{
			
			for(int i1=0; i1<NK; i1++)
				for(int i2=0; i2<NK; i2++)
				for(int i3=0; i3<3; i3++)
				  sum_int_c_c[i1][i2] += eigenvectors[i1][3*i+i3] * eigenvectors[i2][3*i+i3];
		}

	}


	for(int i=0; i<nresil; i++)
	{

		if(ptin_lig[i] == 1)
		{
			
			for(int i1=0; i1<NK; i1++)
				for(int i2=0; i2<NK; i2++)
				for(int i3=0; i3<3; i3++)	
				sum_int_c_c[i1][i2] += eigenvectors[i1][3*(i+nresir)+i3] * eigenvectors[i2][3*(i+nresir)+i3];

		}

	}

	
	for(int i1=0; i1<NK; i1++){
				for(int i2=0; i2<NK; i2++)
					cout<<sum_int_c_c[i1][i2]<<" "; cout<<endl;}


}

//-------------------------------------------choossing putative interface atoms
// randomly perturb starting structures 100 times 
void choo_puta_inte_atom()
{
	cout<<"choosing putative interface atoms..."<<endl;

	memset(ptin_rec, 0 ,sizeof(ptin_rec));
	memset(ptin_lig, 0 ,sizeof(ptin_lig));

	double x_new[NK];
	double phi[NK-1];
	int a1;	
	int a2;
	int hl;
	int hr;


	for(int ncall=0; ncall<50; ncall++)
	{	
		
		sample_coefficient(x_new);
		//cout<<rmsd_receptor_prediction<<endl;
		//memset(x_new,0,sizeof(x_new));
		do_sampling(x_new,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
				current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");
		

		hl = 0;
		hr = 0;
		
		for(int i=0; i<natomr; i++)
		if(atom_name_rec[i][0]!='H')
			hr++;

		for(int i=0; i<natoml; i++)
		if(atom_name_lig[i][0]!='H')
		    hl++;


		for(int i=0; i<natomr; i++)
			for(int j=0; j<natoml; j++)
				if (atom_name_rec[i][0]!='H' && atom_name_lig[j][0]!='H' && distance(sample_rec[i], sample_lig[j], 3) < cutoff_iRMSD)
				{
					ptin_rec[residue_ind_rec[i]] = 1;
					ptin_lig[residue_ind_lig[j]] = 1;
				}
	//		break;
	
	}

	cout<<hr<<" "<<hl<<"#heavy_atoms"<<endl;

	a1 = a2 =0;
	for(int i=0; i<nresir; i++)
		if(ptin_rec[i] ==1 ) a1++;
	for(int i=0; i<nresil; i++)
		if(ptin_lig[i] ==1 ) a2++;

	cout<<"total "<<a1<<" residues choosing as interface residues for receptor"<<endl;
	cout<<"total "<<a2<<" residues choosing as interface residues for ligand"<<endl;
	
	int_nresidue=a1+a2;

	a1 = a2 = 0;
	for(int i=0; i<natomr; i++)
		if(atom_name_rec[i][0]!='H' && ptin_rec[residue_ind_rec[i]]==1)a1++;
	for(int i=0; i<natoml; i++)
		if(atom_name_lig[i][0]!='H' && ptin_lig[residue_ind_lig[i]]==1)a2++;

	cout<<"total "<<a1<<" heavy atoms choosing as interface atoms for receptor"<<endl;
	cout<<"total "<<a2<<" heavy atoms choosing as interface atoms for ligand"<<endl;

	

	cout<<"done with interface-atom choosing process..."<<endl;

	
}



int main()
{
        srand(time(NULL));
        preprocess();
        constant_setting();
	    bayesiandock();


        double x[12]={-0.0677,  0.0327, -0.3143,  0.4218, -0.3891,  0.4027,  0.4881, -0.6794, -0.3804 , 0.5015,  0.4759 , 0.0162};
        double y[12]={0,0,0,0,0,0,0,0,0,0,0,0};
        double feature[1000][10];
        double label[1000];
        double scores[1000];
        int ncall=20;


        FILE *coor_save = fopen(output_path"/Result/samples", "w");
//      fstream fout;
//      fout.open("/home/cyppsp/x_y",ios::out);
        for(int i=0; i<ncall; i++)
        {

              
                sample_coefficient(x);

                if(i==0) for(int ii=0; ii<NK; ii++) x[ii]=0;

                do_sampling(x,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
                                current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");

//              system(("cp "current_directory"/sampling/ligand.pdb "output_path"/Result/smpl_"+to_string(i+1)).c_str());
//              system(("cp "current_directory"/sampling/receptor.pdb "output_path"/Result/smpr_"+to_string(i+1)).c_str());
                                scores[i] = scoring(feature[i]);
                label[i] = true_RMSD(current_directory"/minmiz/");
                //fout<<"beforermsd"<<label[i]<<" "<<true_RMSD(current_directory"/minmiz/")<<endl;
        //      cout<<"beforermsd"<<label[i]<<" "<<true_RMSD(current_directory"/minmiz/")<<endl;
                        //-----------------------------------------------------save coordinate-------------------------


                for(int j=0; j<NK; j++)
                        fprintf(coor_save, "%8.4lf", x[j]);
                fprintf(coor_save, "\n");

                if(i==0)
                        if(label[i]>4) ncall=10;
                        else  ncall=100;


        }
        fclose(coor_save);

   		fstream fout;
        fout.open(output_path"/Result/training.dat",ios::out);

        for(int i=0; i<ncall; i++)
        {
                for(int j=0; j<8; j++)
                        fout<<feature[i][j]<<" ";
                fout<<nresir<<"       "<<label[i]<<endl;
        }

        fout.close();

        fout.open(output_path"/Result/test_score.dat",ios::out);
         for(int i=0; i<ncall; i++)
        {
                fout<<scores[i]<<endl;
        }

        fout.close();
}
