/*
This file is to modify the starting pdb files from ZDOCK  and generate charmm input scripts for hbuild and minmiz.

1.  multiple residues at the same locations
2.  Non-20 common residues 
3.  Missing residues (reindex)
4.  #residue>1000

1. GLYP and PROP for specific chains
2. receptor and ligand chain definition

*/

#include<fstream>
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<sstream>
#include<iomanip>
using namespace std;
#include "configuration.h"

void modify_charmm(char *path_hbuild, char *path_minimiz, char rec_chain, char lig_chain, string *nstep);
void modify_charmm_unbound(char *path_hbuild, char *path_minimiz, char rec_chain, char lig_chain, string *nstep);
double str_num(string s)
{
	int ans=0;
	for(int i=0; i<s.length(); i++)
	{
		if(s[i]==' ') continue;
		ans = ans*10 + (int(s[i]))- (int('0'));
	}
	return ans;
}
int main(int argc, char **argv)
{

	//----------------------------------------------------------------------generate working dir----------------
	system("rm -rf "current_directory"/prepocess");
	system("rm -rf "current_directory"/sampling");
	system("rm -rf "current_directory"/minmiz");
	system("rm -rf "current_directory"/cNMA_output");
	system("rm -rf "output_path"/Result");
	system("rm -rf "current_directory"/stdout");
	system("mkdir "current_directory"/prepocess");
	system("mkdir "current_directory"/sampling");
	system("mkdir "current_directory"/minmiz");
	system("mkdir "current_directory"/cNMA_output");
	system("mkdir "output_path"/Result");
	system("sed -i %s/aaa/"depd_dir"/g compile");
//----------------------------------------------------------------------------
	char file[1000];
	string nstep[1000];
	char rec_chain;
	string ss;
	char c, p ; 
	int residue_ind, t;
	

		system("grep 'ATOM' "protein_path""protein_name"_r_u.pdb | egrep 'ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|HSE|HSD|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL' > temp.dat");
	
		ifstream rec("temp.dat", ios::in);
		

		ofstream outrec(protein_path""protein_name"_r_u.pdb", ios::out);
		
		t = -1;
		residue_ind = 0;
		c='~';
		p='A'-1;

		while( getline(rec, ss) )
		{	
			if (ss[26]!=' ') continue;

			int sub = str_num(ss.substr(22,4));
			if(sub != t)    //--------------------------------for residue index
			{
				t = sub;
				residue_ind ++ ;
			}

			if(ss[21]!= c)   //--------------------------------for chain ID
			{ 
				c=ss[21]; p++; 
					//---------------------------GLYP PROP NSTEP----------------------
				if(ss.substr(17,3)=="GLY")
					nstep[int(p-'A')] = "GLYP";
				else
					if(ss.substr(17,3)=="PRO")
					 nstep[int(p-'A')] = "PROP";
				else
					 nstep[int(p-'A')] = "NTER";
			}
			ss[21]=p;
			
			outrec<<ss.substr(0,22);
			outrec<<right<<setw(4)<<residue_ind;
			outrec<<ss.substr(26, -1)<<endl;
		}
             
		rec_chain = p;


//-------------------------------------------------------modify ligand file-----------------
//--------------------------------------------------------------------------------------------

		system("grep 'ATOM' "protein_path""protein_name"_l_u.pdb | egrep 'ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|HSE|HSD|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL' > temp.dat");
	
		ifstream lig("temp.dat", ios::in);
		

		ofstream outlig(protein_path""protein_name"_l_u.pdb", ios::out);
		
		t = -1;
		residue_ind = 0;
		c='~';


		while( getline(lig, ss) )
		{
			if (ss[26]!=' ') continue;

			int sub = str_num(ss.substr(22,4));
			if(sub != t)
			{
				t = sub;
				residue_ind ++ ;
			}

			if(ss[21]!= c) 
			{ 
				c=ss[21]; p++; 

				if(ss.substr(17,3)=="GLY")
					nstep[int(p-'A')] = "GLYP";
				else
					if(ss.substr(17,3)=="PRO")
					 nstep[int(p-'A')] = "PROP";
				else
					 nstep[int(p-'A')] = "NTER";
			}
			ss[21]=p;
			
			outlig<<ss.substr(0,22);
			outlig<<right<<setw(4)<<residue_ind;
			outlig<<ss.substr(26, -1)<<endl;
		}


		rec.close();
		outrec.close();
		lig.close();
		outlig.close();
		system("rm -f temp.dat");

		modify_charmm(current_directory"/prepocess/hbuild.inp", current_directory"/minmiz/minmiz.inp", rec_chain, p, nstep);
		modify_charmm_unbound(current_directory"/prepocess/urec.inp", current_directory"/prepocess/ulig.inp", rec_chain, p, nstep);
 	return 0;
}

void modify_charmm_unbound(char *path_hbuild, char *path_minimiz, char rec_chain, char lig_chain, string *nstep)
{
 //----------------------------------------modify path_hbuild-----------------------------------------

	ofstream fout(path_hbuild);
//	cout<<path_hbuild<<" "<<path_minimiz<<endl;
//	ofstream f1("hehe");
//	f1<<1;
//	f1.close();
	fout<<"! =========================================="<<endl;
	fout<<"! Read in topology and parameter files"<<endl;
	fout<<"! ===================================="<<endl;
	fout<<endl;
	fout<<"BOMLEV -2"<<endl;
	fout<<endl;
	fout<<"open read unit 1 card name \""<<charmm_top<<"\""<<endl;
	fout<<"read rtf card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;
	fout<<endl;
	fout<<"open read unit 1 card name \""<<charmm_para<<"\""<<endl;
	fout<<"read para card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;

//------------------------------------------read chain files
	for(char c='A'; c<=rec_chain; c++)
	{
	fout<<"open read unit 1 card name "<<c<<".SEQ"<<endl;
	fout<<"read sequence unit 1 card"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<"generate "<<c<<" first "<<nstep[int(c-'A')]<<" last CTER setup"<<endl;
	fout<<"stream \""<<c<<"_FIXRES.INP\""<<endl;
	fout<<endl;
	}
//---------------------------------------------------------

	fout<<"open read unit 1 card name \"un_receptor.crd\""<<endl;
	fout<<"read coor card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;
	fout<<endl;
	fout<<"ic fill preserve"<<endl;
	fout<<"ic param"<<endl;
	fout<<"ic build"<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"coor init sele hydrogen end"<<endl;
	fout<<endl;
	fout<<"! Re-build hydrogens"<<endl;
	fout<<endl;

//---------------------------------rebuild hydrogens-----------------------------------------------
	for(int i=0; i<2; i++)
	{
			fout<<"hbuild select hydrogen end -"<<endl;
			fout<<"       electrostatic atom cdiel eps 4.0 -"<<endl;
			fout<<"       switch vdw vatom vswitch cutnb 299.0 ctofnb 298.0 ctonnb 297.0 -"<<endl;
			fout<<"       wmin 1.5 e14fac 1 nbxmod 5"<<endl;
			fout<<endl;
	}
	fout<<endl;

	fout<<"coor copy comp"<<endl;
	fout<<endl;


	fout<<"! calculate starting energy"<<endl;
	fout<<"update -"<<endl;
	fout<<"energy elec atom rdiel switch eps 4.0 e14fac 1.0 -"<<endl;
	fout<<"                cutnb 14. ctonnb 10. ctofnb 12. -"<<endl;
	fout<<"                nbxmod 5 vswitch wmin 1.5 vatom vdistance"<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"skip all excl bond angl urey dihe impr"<<endl;
	fout<<"mini sd nstep 100"<<endl;
	fout<<"skip all excl bond angl urey dihe impr vdw vatt vrep"<<endl;	
	fout<<endl;

	fout<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<endl;

	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;

	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;
	
	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;
	
	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;

	fout<<"mini abnr nstep 100"<<endl;
	fout<<endl;
	fout<<"skip none"<<endl;

	fout<<"prnlev 0"<<endl;
	fout<<"stream "<<charmm_gbsw<<endl;
	fout<<"prnlev 5"<<endl;
	fout<<"scalar wmain statistics select .not. type H* end"<<endl;
	fout<<"GBSW sw 0.3 sgamma 0.03 dgp 1.5 GBenergy"<<endl;
	fout<<"Energy"<<endl;

	fout<<"calc cinternal ( ?bond + ?angl + ?dihe + ?impr + ?urey )"<<endl;
        fout<<"calc bond ?bond"<<endl;
        fout<<"calc angl ?angl"<<endl;
        fout<<"calc dihe ?dihe"<<endl;
        fout<<"calc impr ?impr"<<endl;
        fout<<"calc urey ?urey"<<endl;
        fout<<"calc cvdw ?vdw"<<endl;
        fout<<"calc elec ?gben"<<endl;
        fout<<"coor surf acce rprobe 1.4 sele all END"<<endl;
        fout<<"calc csasa ( 0.005 * ?area + 0.860)"<<endl;
        fout<<endl;


        fout<<"open write unit 10 card name rec_energy.dat"<<endl;
        fout<<"write title unit 10"<<endl;
        fout<<"* @bond @angl @dihe @impr @urey @cvdw @csasa @elec"<<endl;
        fout<<"*"<<endl;
        fout<<"close unit 10"<<endl;
        fout<<endl;
        fout<<endl;

	fout.close();


//-------------------------------------------------------------------------write minimiz.inp----------------------
//-------------------------------------------------------------------------write minimiz.inp----------------------
//-------------------------------------------------------------------------write minimiz.inp----------------------

	fout.open(path_minimiz);

	fout<<"! =========================================="<<endl;
	fout<<"! Read in topology and parameter files"<<endl;
	fout<<"! ===================================="<<endl;
	fout<<endl;
	fout<<"BOMLEV -2"<<endl;
	fout<<endl;
	fout<<"open read unit 1 card name \""<<charmm_top<<"\""<<endl;
	fout<<"read rtf card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"open read unit 1 card name \""<<charmm_para<<"\""<<endl;
	fout<<"read para card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;

//------------------------------------------read chain files
	for(char c=char(rec_chain+1); c<=lig_chain; c++)
	{
	fout<<"open read unit 1 card name "<<c<<".SEQ"<<endl;
	fout<<"read sequence unit 1 card"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<"generate "<<c<<" first "<<nstep[int(c-'A')]<<" last CTER setup"<<endl;
	fout<<"stream \""<<c<<"_FIXRES.INP\""<<endl;
	fout<<endl;
	}
//---------------------------------------------------------

	fout<<"open read unit 1 card name \"un_ligand.crd\""<<endl;
	fout<<"read coor card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;
	fout<<endl;
	fout<<"ic fill preserve"<<endl;
	fout<<"ic param"<<endl;
	fout<<"ic build"<<endl;
	fout<<endl;
	fout<<endl;


	fout<<"coor init sele hydrogen end"<<endl;
	fout<<endl;
	fout<<"! Re-build hydrogens"<<endl;
	fout<<endl;

//---------------------------------rebuild hydrogens-----------------------------------------------
	for(int i=0; i<2; i++)
	{
			fout<<"hbuild select hydrogen end -"<<endl;
			fout<<"       electrostatic atom cdiel eps 4.0 -"<<endl;
			fout<<"       switch vdw vatom vswitch cutnb 299.0 ctofnb 298.0 ctonnb 297.0 -"<<endl;
			fout<<"       wmin 1.5 e14fac 1 nbxmod 5"<<endl;
			fout<<endl;
	}
	fout<<endl;

	fout<<"coor copy comp"<<endl;
	fout<<endl;


	fout<<"! calculate starting energy"<<endl;
	fout<<"update -"<<endl;
	fout<<"energy elec atom rdiel switch eps 4.0 e14fac 1.0 -"<<endl;
	fout<<"                cutnb 14. ctonnb 10. ctofnb 12. -"<<endl;
	fout<<"                nbxmod 5 vswitch wmin 1.5 vatom vdistance"<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"skip all excl bond angl urey dihe impr"<<endl;
	fout<<"mini sd nstep 100"<<endl;
	fout<<"skip all excl bond angl urey dihe impr vdw vatt vrep"<<endl;	
	fout<<endl;

	fout<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<endl;

	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;

	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;
	
	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;
	
	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;

	fout<<"mini abnr nstep 100"<<endl;
	fout<<endl;
	fout<<"skip none"<<endl;
	
	fout<<"prnlev 0"<<endl;
	fout<<"stream "<<charmm_gbsw<<endl;
	fout<<"prnlev 5"<<endl;
	fout<<"scalar wmain statistics select .not. type H* end"<<endl;
	fout<<"GBSW sw 0.3 sgamma 0.03 dgp 1.5 GBenergy"<<endl;
	fout<<"Energy"<<endl;

	fout<<"calc cinternal ( ?bond + ?angl + ?dihe + ?impr + ?urey )"<<endl;
	fout<<"calc bond ?bond"<<endl;
	fout<<"calc angl ?angl"<<endl;
	fout<<"calc dihe ?dihe"<<endl;
	fout<<"calc impr ?impr"<<endl;
	fout<<"calc urey ?urey"<<endl;
	fout<<"calc cvdw ?vdw"<<endl;
	fout<<"calc elec ?gben"<<endl;
	fout<<"coor surf acce rprobe 1.4 sele all END"<<endl;
	fout<<"calc csasa ( 0.005 * ?area + 0.860)"<<endl;
	fout<<endl;


	fout<<"open write unit 10 card name lig_energy.dat"<<endl;
	fout<<"write title unit 10"<<endl;
	fout<<"* @bond @angl @dihe @impr @urey @cvdw @csasa @elec"<<endl;
	fout<<"*"<<endl;
	fout<<"close unit 10"<<endl;
	fout<<endl;
	fout<<endl;

	fout.close();



}




void modify_charmm(char *path_hbuild, char *path_minimiz, char rec_chain, char lig_chain, string *nstep)
{
 //----------------------------------------modify path_hbuild-----------------------------------------

	ofstream fout(path_hbuild);

	fout<<"! =========================================="<<endl;
	fout<<"! Read in topology and parameter files"<<endl;
	fout<<"! ===================================="<<endl;
	fout<<endl;
	fout<<"BOMLEV -2"<<endl;
	fout<<endl;
	fout<<"open read unit 1 card name \""<<charmm_top<<"\""<<endl;
	fout<<"read rtf card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"open read unit 1 card name \""<<charmm_para<<"\""<<endl;
	fout<<"read para card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;

//------------------------------------------read chain files
	for(char c='A'; c<=lig_chain; c++)
	{
	fout<<"open read unit 1 card name "<<c<<".SEQ"<<endl;
	fout<<"read sequence unit 1 card"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<"generate "<<c<<" first "<<nstep[int(c-'A')]<<" last CTER setup"<<endl;
	fout<<"stream \""<<c<<"_FIXRES.INP\""<<endl;
	fout<<endl;
	}
//---------------------------------------------------------

	fout<<"open read unit 1 card name \"complex.crd\""<<endl;
	fout<<"read coor card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;
	fout<<endl;
	fout<<"ic fill preserve"<<endl;
	fout<<"ic param"<<endl;
	fout<<"ic build"<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"coor init sele hydrogen end"<<endl;
	fout<<endl;
	fout<<"! Re-build hydrogens"<<endl;
	fout<<endl;

//---------------------------------rebuild hydrogens-----------------------------------------------
	for(int i=0; i<2; i++)
	{
			fout<<"hbuild select hydrogen end -"<<endl;
			fout<<"       electrostatic atom cdiel eps 4.0 -"<<endl;
			fout<<"       switch vdw vatom vswitch cutnb 299.0 ctofnb 298.0 ctonnb 297.0 -"<<endl;
			fout<<"       wmin 1.5 e14fac 1 nbxmod 5"<<endl;
			fout<<endl;
	}
	fout<<endl;

	fout<<"coor copy comp"<<endl;
	fout<<endl;

	fout<<"define REC sele segid A";
	for(char c='B'; c<=rec_chain; c++)
		fout<<" .OR. segid "<<c;
	fout<<" "<<"end"<<endl;

	fout<<"define LIG sele segid "<<char(rec_chain+1);
	for(char c=char(rec_chain+2); c<=lig_chain; c++)
		fout<<" .OR. segid "<<c;
	fout<<" "<<"end"<<endl;
	fout<<endl;

	fout<<"define INT sele .byres. LIG .around. 10 end"<<endl;
	fout<<"define RINT sele INT .and. REC end"<<endl;
	fout<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"! calculate starting energy"<<endl;
	fout<<"update -"<<endl;
	fout<<"energy elec atom rdiel switch eps 4.0 e14fac 1.0 -"<<endl;
	fout<<"                cutnb 14. ctonnb 10. ctofnb 12. -"<<endl;
	fout<<"                nbxmod 5 vswitch wmin 1.5 vatom vdistance"<<endl;
	fout<<endl;
	fout<<endl;

	
	fout<<"skip all excl bond angl urey dihe impr"<<endl;
	fout<<"mini sd nstep 100"<<endl;
	fout<<"skip all excl bond angl urey dihe impr vdw "<<endl;	


	
	fout<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<endl;

	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;

	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;
	
	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;
	
	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;

	fout<<"mini abnr nstep 100"<<endl;
	fout<<endl;
	fout<<"!cons fix sele none end"<<endl;
	fout<<endl;
	fout<<"skip none"<<endl;

	fout<<endl;
	fout<<endl;

	fout<<"coor orie rms sele LIG end"<<endl;
	fout<<"inte sele REC end sele LIG end"<<endl;
	fout<<"scalar wmain = charge"<<endl;
	fout<<endl;
	fout<<endl;

//---------------------------------------output receptor file and ligand file-----------------------
	fout<<"open write unit 1 card name \"receptor.crd\""<<endl;
	fout<<"write coor card unit 1 sele REC end"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;

	fout<<"open write unit 1 card name \"ligand.crd\""<<endl;
	fout<<"write coor card unit 1 sele LIG end"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;

//---------------------------------------------------------------output energy-----------------------

	
	fout.close();


//-------------------------------------------------------------------------write minimiz.inp----------------------
//-------------------------------------------------------------------------write minimiz.inp----------------------
//-------------------------------------------------------------------------write minimiz.inp----------------------

	fout.open(path_minimiz);

	fout<<"! =========================================="<<endl;
	fout<<"! Read in topology and parameter files"<<endl;
	fout<<"! ===================================="<<endl;
	fout<<endl;
	fout<<"BOMLEV -2"<<endl;
	fout<<endl;
	fout<<"open read unit 1 card name \""<<charmm_top<<"\""<<endl;
	fout<<"read rtf card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"open read unit 1 card name \""<<charmm_para<<"\""<<endl;
	fout<<"read para card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;

//------------------------------------------read chain files
	for(char c='A'; c<=lig_chain; c++)
	{
	fout<<"open read unit 1 card name "<<c<<".SEQ"<<endl;
	fout<<"read sequence unit 1 card"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<"generate "<<c<<" first "<<nstep[int(c-'A')]<<" last CTER setup"<<endl;
	fout<<"stream \""<<c<<"_FIXRES.INP\""<<endl;
	fout<<endl;
	}
//---------------------------------------------------------

	fout<<"open read unit 1 card name \"complex.crd\""<<endl;
	fout<<"read coor card unit 1"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;
	fout<<endl;
	fout<<"ic fill preserve"<<endl;
	fout<<"ic param"<<endl;
	fout<<"ic build"<<endl;
	fout<<endl;
	fout<<endl;


	fout<<"coor copy comp"<<endl;
	fout<<endl;

	fout<<"define REC sele segid A";
	for(char c='B'; c<=rec_chain; c++)
		fout<<" .OR. segid "<<c;
	fout<<" "<<"end"<<endl;

	fout<<"define LIG sele segid "<<char(rec_chain+1);
	for(char c=char(rec_chain+2); c<=lig_chain; c++)
		fout<<" .OR. segid "<<c;
	fout<<" "<<"end"<<endl;
	fout<<endl;

	fout<<"define INT sele .byres. LIG .around. 10 end"<<endl;
	fout<<"define RINT sele INT .and. REC end"<<endl;
	fout<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"! calculate starting energy"<<endl;
	fout<<"update -"<<endl;
	fout<<"energy elec atom rdiel switch eps 4.0 e14fac 1.0 -"<<endl;
	fout<<"                cutnb 14. ctonnb 10. ctofnb 12. -"<<endl;
	fout<<"                nbxmod 5 vswitch wmin 1.5 vatom vdistance"<<endl;
	fout<<endl;
	fout<<endl;

	//fout<<"shape desc RIGILIG RIGI sele LIG end"<<endl;
//	fout<<"skip all excl bond angl urey dihe impr"<<endl;
//	fout<<"mini sd nstep 100"<<endl;
//	fout<<"skip all excl bond angl urey dihe impr vdw vatt vrep"<<endl;	

//-------------using "shape" will cause large rigid-body motion of receptor, which leading to the unaccpetable structure.
//	fout<<"shape desc RIGILIG RIGI sele LIG end"<<endl;
//	fout<<"shape desc RIGIREC RIGI sele REC end"<<endl;
	fout<<"mini sd nstep 50"<<endl;
//	fout<<"shape clear"<<endl;

	
	fout<<endl;
	fout<<"mini sd nstep 100"<<endl;
	fout<<endl;

	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;

	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;
	
	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;
	
	fout<<"if ?vdw GT -200 then"<<endl;
	fout<<"mini sd nstep 50"<<endl;
	fout<<"endif"<<endl;
	fout<<endl;

	fout<<"mini abnr nstep 100"<<endl;
	fout<<endl;
	fout<<"!cons fix sele none end"<<endl;
	fout<<endl;
	fout<<"skip none"<<endl;
	//fout<<"shape clear"<<endl;
	fout<<endl;
	fout<<endl;

	fout<<"coor orie rms sele LIG end"<<endl;
	fout<<"inte sele REC end sele LIG end"<<endl;
	fout<<"scalar wmain = charge"<<endl;
	fout<<endl;
	fout<<endl;

//---------------------------------------output receptor file and ligand file-----------------------
	fout<<"open write unit 1 card name \"receptor.crd\""<<endl;
	fout<<"write coor card unit 1 sele REC end"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;

	fout<<"open write unit 1 card name \"ligand.crd\""<<endl;
	fout<<"write coor card unit 1 sele LIG end"<<endl;
	fout<<"close unit 1"<<endl;
	fout<<endl;

//---------------------------------------------------------------output energy-----------------------

	fout<<"prnlev 0"<<endl;
	fout<<"stream "<<charmm_gbsw<<endl;
	fout<<"prnlev 5"<<endl;
	fout<<"scalar wmain statistics select .not. type H* end"<<endl;
	fout<<"GBSW sw 0.3 sgamma 0.03 dgp 1.5 GBenergy"<<endl;
	fout<<"Energy"<<endl;

	fout<<"calc cinternal ( ?bond + ?angl + ?dihe + ?impr + ?urey )"<<endl;
        fout<<"calc bond ?bond"<<endl;
        fout<<"calc angl ?angl"<<endl;
        fout<<"calc dihe ?dihe"<<endl;
        fout<<"calc impr ?impr"<<endl;
        fout<<"calc urey ?urey"<<endl;
        fout<<"calc cvdw ?vdw"<<endl;
        fout<<"calc elec ?gben"<<endl;
        fout<<"coor surf acce rprobe 1.4 sele all END"<<endl;
        fout<<"calc csasa ( 0.005 * ?area + 0.860)"<<endl;
        fout<<endl;


        fout<<"open write unit 10 card name com_energy.dat"<<endl;
        fout<<"write title unit 10"<<endl;
        fout<<"* @bond @angl @dihe @impr @urey @cvdw @csasa @elec"<<endl;
        fout<<"*"<<endl;
        fout<<"close unit 10"<<endl;
        fout<<endl;
        fout<<endl;

	fout.close();


}
