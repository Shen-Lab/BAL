# Animation of protein motions in BAL

All files are in the format of mp4.  There are two types of movies: 

## complex normal modes (by cNMA) 

Files are named {PDB ID}_m{Starting Model Index}_N{cNMA Normal Mode Index}. Receptors are in cyan cartoons and ligands are in orange cartoons. 

In these examples for 4JW2 (Starting Model 7 from ZDOCK unbound docking), we see that the first 3 complex normal modes are dominated by ligand motions, whereas the 6th, 7th, and 9th normal modes include substantial receptor deformations along with ligand motions.  

How these movies of normal modes are made: We show the motion of each protein complex moving along each designated complex normal mode.  The extent is currently set at the predicted.   

In BAL sampling, we use the first K1 complex normal modes (ranked by the whole-complex eigen values) plus K2 more (ranked by those eigen values after being rescaled to the contribution of receptors alone; as in Eq. 7 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4765865/).  Currently K1 is set at 9 (to approximately include 6 involving ligand rigid motions and 3 involving ligand flexibility motions) and K2 is set at 3 (to include those involving receptor flexibility motions).  For instance, in the case of 4JW2_m7, the 6th, 7th, and 9th normal modes are actually the slowest normal modes for the receptor portion (after aforementioned re-scaling).  Of course, in this case, they would have been included in the top K1=9 complex normal modes anyway.  So the other K2=3 of the total 12 normal modes would be chosen further down the list of receptor-dominant normal modes.   


## BAL sampling trajectories 

Files are named {PDB ID}_m{Starting Model Index}.mp4. Receptors are in cyan cartoons and ligands are in orange cartoons. At the end of each move, bound conformations are shown in gray.  

In the four sampling examples given, 3 are about the refinement of near-native starting models.  The other (4CPA_m4) starts with a non-native starting model.   

How these sampling movies are made: We record the best sample (of the lowest energy) at the end of each BAL iteration and show the motions through linear interpolation between the ith and the (i+1)th such samples.  Sometimes the movies appear to be "frozen" but actually it is because not every iteration sees a lower-energy sample.  

PS:  The flashing part in the video is due to the style of Pymol visualization for secondary structures. 


