#include<fstream>
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

using namespace std;

int main(int argc, char **argv)
{
	fstream Kd_list, Kd_zero, out;

	Kd_list.open("kd_list", ios::in);
	Kd_zero.open("kd_zero", ios::in);
	out.open("Kd_value", ios::out);

	char protein[4];
	double kd;
	int i=1;
	while(i<200)
	{
		Kd_list>>protein;
		Kd_zero>>kd;

		if(!strcmp(protein,argv[1])){

				out<<kd;
				break;
		}
		i++;

	}
	out<<kd;

        Kd_list.close();
        Kd_zero.close();
        out.close();
	

	return 0;
}
