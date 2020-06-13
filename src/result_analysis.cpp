

#include<fstream>
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <sstream>
#include "configuration.h"
using namespace std;

void posterior_anaylsis(double *x)
{
	double feature[8];
	double label;



//---------------------------------------------------------end structure---------------------------------------------
				
	
	do_sampling(x,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
				current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");
	



	scoring(feature);
	system("mkdir Result/end");
        system("cp minmiz/receptor.pdb Result/end/");
        system("cp minmiz/ligand.pdb Result/end/");
        system("python pre_bal.py Result/end/receptor.pdb Result/end/ligand.pdb out");

//--------------------------------------starting structure--------------------------------
	for(int i=0; i<NK; i++) x[i]=0;

	do_sampling(x,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
				current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");
	

	scoring(feature);
	
	system("mkdir Result/start");
	system("cp minmiz/receptor.pdb Result/start/");
	system("cp minmiz/ligand.pdb Result/start/");
        system("python pre_bal.py Result/start/receptor.pdb Result/start/ligand.pdb out");
}
