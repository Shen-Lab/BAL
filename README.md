# BAL
Bayesian Active Learning for Optimization and Uncertainty Quantification with Applications in Protein Docking
https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c00476


# Versions in branch:
* master: The standard BAL
* noclashv1: The version whose output structure will gurantee \delta vdw<0.


## Dependencies:
* C++ 4.8.5 or higher
* cNMA: Download and install the cNMA from "https://github.com/Shen-Lab/cNMA".
* Energy model:  Please download the random forest energy model from 
https://drive.google.com/open?id=17ByuNoYy0t1R8EjuTK_cMyul5K004MHa
* CHARMM: Download the executable 'CHARMM36a1.exe', GBSW: 'radius_gbsw.str', CHARMM27 topology and parameter files and put them into 'dependencies/'.
* Linux Environment

## Change directory path
In src/configuration.h, please change the macros as follow:
* src_dir:  change to your current "src" path.
* cnma_path: change to your cNMA path.
* scoring_path: change to your 'random_forest.sav' path.
* output_path: change to the directory where you want to place the 'Result' folder which contains all output files.
* charmm_exe, charmm_para, charmm_top, charmm_setup, charmm_gbsw, pdb_crd, crd_pdb:  change paths of these macros to their corresponding paths in 'dependencies/'. 

## Add unbound proteins:
In src/configuration.h, please change the macro 'protein_name' to your 4-letter/digit protein code. Also please change the macro 'protein_path" to the path where the unbound protein is.

## If your unbound protein is not in Protein Docking Benchmark 4.0 and you want to get the UQ results, please append its Kd value into 'src/kd_zero' and append the protein name into 'src/kd_list'. Otherwise, you will only get the refined structures and the area under the posterior.

## Compile and Run
* Go to 'BAL/src/'.
* Type './complie' to compile.
* Type './for_train' to run BAL.

## Output file:

* a.'end'(dir):
        containing the refined structures: receptor.pdb ligand.pdb

* b.'scores':
        The final energy (unit in Kcal/mol) of refined structures.

* c.'lgAUC':
        The log of area under the posterior.

* d.'cond_prob.dat':
        The conditional probability of refined structure: P(RMSD(x^,x*)<4 | x* \in Mi)

* e. 'Rmsd_dis':
        The RMSD(x^,x*)  distribution of the posterior

* f. 'UQ':
	The [lb,ub] values

	
## Citation:
```
@article{cao2020bayesian,
  title={Bayesian active learning for optimization and uncertainty quantification in protein docking},
  author={Cao, Yue and Shen, Yang},
  journal={Journal of chemical theory and computation},
  volume={16},
  number={8},
  pages={5334--5347},
  year={2020},
  publisher={ACS Publications}
}
```

## Contact:
Yang Shen: yshen@tamu.edu
