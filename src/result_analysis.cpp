

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
	system("mkdir "output_path"/Result/end");
        system("cp minmiz/receptor.pdb "output_path"Result/end/");
        system("cp minmiz/ligand.pdb "output_path"Result/end/");
	system("cp scores "output_path"Result/");

//--------------------------------------starting structure--------------------------------
	for(int i=0; i<NK; i++) x[i]=0;

	do_sampling(x,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
				current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");
	

	scoring(feature);
	
	system("mkdir "output_path"Result/start");
	system("cp minmiz/receptor.pdb "output_path"Result/start/");
	system("cp minmiz/ligand.pdb "output_path"Result/start/");

}
