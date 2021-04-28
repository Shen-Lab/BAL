

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
	double x_start[NK];
	double e;


//--------------------------------------starting structure--------------------------------
        for(int i=0; i<NK; i++) x_start[i]=0;

        do_sampling(x_start,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
                                current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");


        e = scoring(feature);

	FILE *printfeat=fopen(output_path"/Result/energy_features","a");	
        for (int j=0; j<8; j++)
                                        fprintf(printfeat, "%9.2f", feature[j]);
        fprintf(printfeat, "%9.2f\n", e);


        system("mkdir "output_path"/Result/start");
        system("cp minmiz/receptor.pdb "output_path"/Result/start/");
        system("cp minmiz/ligand.pdb "output_path"/Result/start/");
        system("python pre_bal.py "output_path"/Result/start/receptor.pdb "output_path"/Result/start/ligand.pdb out");





//---------------------------------------------------------end structure---------------------------------------------
					
	do_sampling(x,current_directory"/prepocess/ligand.pdb", current_directory"/sampling/ligand.pdb",
				current_directory"/prepocess/receptor.pdb", current_directory"/sampling/receptor.pdb");
	



	e = scoring(feature);

	        for (int j=0; j<8; j++)
                                        fprintf(printfeat, "%9.2f", feature[j]);
                             fprintf(printfeat, "%9.2f\n", e);
	fclose(printfeat);

	system("mkdir "output_path"/Result/end");
        system("cp minmiz/receptor.pdb "output_path"/Result/end/");
        system("cp minmiz/ligand.pdb "output_path"/Result/end/");
        system("python pre_bal.py "output_path"/Result/end/receptor.pdb "output_path"/Result/end/ligand.pdb out");


}
