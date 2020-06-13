#ifndef configuration_H
#define configuration_H

#define protein_name "80r"
#define protein_path "/scratch/user/cyppsp/project_bayesian/covid_19/ncov_80r/Input/"

#define current_directory "/scratch/user/cyppsp/project_bayesian/covid_19/ncov_80r/"
#define output_path "/scratch/user/cyppsp/project_bayesian/covid_19/ncov_80r/"

#define cnma_path "/scratch/user/cyppsp/project_bayesian/covid_19/ncov_80r/cNMA/"   //---------------------------Here is complex cNMA!!
#define scoring_path "~/random_forest.sav"

#define charmm_exe   "/home/cyppsp/bin/charmm36a1.exe"
#define charmm_para   "/home/cyppsp/param/par_all27_prot_na_drv.v4.1.inp"
#define charmm_top     "/home/cyppsp/param/top_all27_prot_na_drv.v4.1.inp"
#define charmm_setup  "/home/cyppsp/scripts/charmm_setup.prl"
#define charmm_gbsw  "/home/cyppsp/param/radius_gbsw.str"
#define pdb_crd        "/home/cyppsp/scripts/pdb_crd.prl"
#define crd_pdb        "/home/cyppsp/scripts/crd_pdb.prl"


#define charsize 10000

#define NK  12         // Number of normal modes used for sampling
#define numofnodes   400 //  Number of nodes for use
#define residue_max   5000
#define atom_max 20000

//#define lig_thereshold 0.001




//--------------------------------------------------Bayesian and Metropolis algorithm parameters

//------------------------------------------------constant pi and E
#define pi_1 3.1415926
#define E 2.7183
#define oo 2147483647  //----------------------infinite

//#define l 0.55   			   //--------------------------------------width of the Guassian kernel ,    The most important parameter here!!!!

// -----------------------------------------new update for l: l=(h0*n^\alpha)

#define h0 2
#define alfa -1.0/(12.0+5.0)
#define Mu 0                 //----------------------The location parameter of the prior
#define sd_prior 5				   //----------------------The  scaling parameter of the prior
#define epsilon_1  1                 //----------------------- The number of the prior points, 
#define rho_zero 1				  //----------------------The scaling parameter of the posterior,       influence the ac ratio of MC large

//--------------------------------------------------------



#define nosp  700              //----------------------number of sampling points used for generating posterior 

#define Metropolis_Steplenth 0.5 //-------------------- The length of the metropolis_step

#define burning_MCMC 100000   //-----------------The number of the metropolis step
#define pick_MCMC 1000
#define n_sample 20				  //  how many samples are being picked at each iteration
#define iteration 15              

#define sample_LBFGS 10000  //----------------------------The number of points for L-BFGS 
#define equal_precise  1E-3  //----------------------------The thereshold for equal judgement

#define sample_pmi  10000   //-----------------------------The number of points for Importance sampling of P(Mi)
#define rmsd_sample 1000000
//#define number_of_data_points_in_each_interation  10
#define pdb_file 0
#define rigid_rotate 0      //---------------------------rigid rotate over the interface of ligand(1),,,,or, center of ligand(0)

#define rmsd_low  0.7             //-------------These two are the lower and upper bound of RMSD of the shell
#define rmsd_up  1.2			//----------------------------------------------------

//#define range_of_dis_up    0.2
//#define range_of_dis_low  -0.2
#define range_ligand 6.0



class Samples 
{
public:
	double x[nosp][NK];
	double y[nosp];
	double l, rho; 
	int n;
	double coeff[nosp];
	//double current_entropy;
};

extern int nresir;    //-----------------------num of residues of receptor
extern int nresil;	  //-----------------------num of residues of ligand
extern int natomr;	  //-----------------------num of atoms of receptor
extern int natoml;	  //-----------------------num of atoms of ligand
extern int int_nresidue; //-----------------------num of residues of the interface 

//extern double nomo[NK][3*residue_max];  //---------complex normal modes (first 6 normal modes)

extern bool ptin_rec[residue_max];         //---------putative interface residues of receptor
extern bool ptin_lig[residue_max];         //---------putative interface residues of ligand
extern double coor_lig[atom_max][3];       //---------ligand atom coordinates
extern double coor_rec[atom_max][3];		//--------receptor atom coordinates
extern char charge1_lig[atom_max][6];
extern char charge1_rec[atom_max][6];
extern char charge2_lig[atom_max][6];
extern char charge2_rec[atom_max][6];
extern int residue_ind_lig[atom_max];  
extern int residue_ind_rec[atom_max];
extern char chainid_lig[atom_max];
extern char chainid_rec[atom_max];
extern char atom_name_lig[atom_max][5];
extern char atom_name_rec[atom_max][5];
extern char residue_name_lig[atom_max][5]; 
extern char residue_name_rec[atom_max][5];
extern double sample_lig[atom_max][3];
extern double sample_rec[atom_max][3];

//-------------------------------------------------starting individual protein internal energy
extern double un_bond, un_angl, un_dihe, un_impr, un_urey, un_vdw, un_ele_npo, un_sasa;        



//extern double R_L_distance;                      //---------------The distance of two individual proteins
//extern double rotation[3][3];                   //----------------translational rotation matrix
//extern double center_rx, center_ry, center_rz;   //---------------coor center of receptor
//extern double center_lx, center_ly, center_lz;   //---------------coor center of ligand
//extern double center_int_ligx, center_int_ligy,center_int_ligz;
//----------------------------------------------------only for interface heavy atoms, normal modes as same, xyz for ligand-ref      

extern double sum_r_r[NK][NK];               //---------\mu_i^R \cdot \mu_i^R.  this is very important
extern double sum_int_c_c[NK][NK];
//extern double sum_x;
//extern double sum_y;
//extern double sum_z;
//extern double sum_xx;
//extern double sum_xy;
//extern double sum_xz;
//extern double sum_yy;
//extern double sum_yz;
//extern double sum_zz;
//extern double sum_t1;
//extern double sum_t2;
//extern double sum_t3;




extern double rmsd_receptor_prediction;
extern double eigenvalues[numofnodes];
extern double eigenvectors[numofnodes][residue_max];
//extern double radius_soid[NK];
//extern double transform1[NK][NK];
//extern double ellipsoid_radius[NK];


void preprocess();
void sampling(double *x_new);
void constant_setting();
void rigid_coor_transfer(double *coor, double *t, double *w);
void do_sampling(double *coor, char *lpath, char *result_l, char *rpath, char *result_r);
double scoring(double *feature);
double true_RMSD(char *a);
double RMSD_distance (double *coor1, double *coor2);
void bayesian_op();
double uniform_sampling_oneD(double lower, double upper);
void uniform_sampling_multiD(double *temp);
double normal_sampling(double u, double o);
double distance(double *x, double *y, int d);
double rmsd_distance1(double *x1, double *x2);
void choo_puta_inte_atom();
void charmm_minmiz();
bool judgebound(double *temp);
void swap1(double *a, double *b);

double K(double *x1, double *x2, double l);
double stdK(double *x1, double *x2, double l);
double precision_prior(double *x_new);
double prior(double *x_new);
double gradient_K(double *xi, double l, double *x, int n);
double gradient_prior(double *x, int n);
double gradient_k0(double *x, int n);
double uniform_sampling_oneD(double lower, double upper);
double normal_sampling(double u, double o);
double distance(double *x, double *y, int d);
void uniform_sampling_multiD(double *temp);
double absolute(double x);
double NW_estimator(Samples sample, double *x_new);
double rho_calc(Samples sample);
double kriging(Samples sample, double *x_new);
double coeff_cal(Samples *sample);
double rho_entropy(double (*pos_sample)[NK], int n, Samples sample);
void bayesiandock();

#define cutoff_iRMSD 10
//double rmsd_prediction(double *rigid_x,  double *lc);

void uncertainty(Samples sample, double *x_center, double stepsize);
void posterior_anaylsis(double *x);
void post_sdy();
int compare (const void * a, const void * b);



void sample_coefficient(double *temp);
double l2_norm(double *y);
void cart_sph_conver(double *coor, double *phi);
void sph_cart_conver(double *coor, double *phi);

#define confidence_level 0.90

#endif
