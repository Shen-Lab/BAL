# BAL
Bayesian Active Learning for Optimization and Uncertainty Quantification with Applications in Protein Docking

## Dependencies:
* C++ 4.8.5 or higher
* cNMA:  Download and install the cNMA from "https://github.com/Shen-Lab/cNMA".
* Energy model:  Please download the random forest energy model from 
"https://drive.google.com/file/d/17ByuNoYy0t1R8EjuTK_cMyul5K004MHa/view?usp=sharing"
* CHARMM27: Download the executable 'CHARMM36a1.exe', GBSW: 'radius_gbsw.str', CHARMM27 topology and parameter files and put them into 'dependencies/'.
* Linux Environment

## Change directory path
In src/configuration.h, please change the macros as follow:
* src_dir:  change to your current "src" path.
* cnma_path: change to your cNMA path.
* scoring_path: change to your 'random_forest.sav' path.
* output_path: change to the directory where you want to output.

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

* c.'PMI_log':
        The log of area under the posterior.

* d.'cond_prob.dat':
        The conditional probability of refined structure: P(RMSD(x^,x*)<4 | x* \in Mi)

* e. 'Rmsd_dis':
        The RMSD(x^,x*)  distribution of the posterior

* f. 'UQ':
	The [lb,ub] values
	
## Citation:

@article{Cao537035,
        author = {Cao, Yue and Shen, Yang},
        title = {Bayesian Active Learning for Optimization and Uncertainty Quantification in Protein Docking},
        elocation-id = {537035},
        year = {2019},
        doi = {10.1101/537035},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2019/01/31/537035},
        eprint = {https://www.biorxiv.org/content/early/2019/01/31/537035.full.pdf},
        journal = {bioRxiv}
	}

## Contact:
Yang Shen: yshen@tamu.edu
