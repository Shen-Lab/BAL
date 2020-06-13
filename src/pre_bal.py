import sys
import pickle
import os


def encoder():


	pdb_rec = sys.argv[1]
	pdb_lig = sys.argv[2]
	mapping={}
	def encoder_individual(new_chain, path):

		

		pre_chain='hehe'
		pre_resi=-1

		
		new_resi=0
		with open("temp", "w") as fo:
			with open(path, "r") as f:
				for lines in f:
					if lines[0:4]!='ATOM':
						continue
					key = lines[21:27]
					chain= lines[21]
					resi_ind = lines[22:27]

					

					if chain!=pre_chain:
						pre_chain=chain
						pre_resi=resi_ind
						new_resi+=1
						new_chain =chr(ord(new_chain)+1)
					elif resi_ind!=pre_resi:
						pre_resi=resi_ind
						new_resi+=1


					fout = "%s%c%4d %s" %(lines[0:21],  new_chain, new_resi, lines[27:])
					fo.write(fout)
					#print (fout, end='')
					mapping[fout[21:27]] =lines[21:27]
					
		os.system("mv temp "+path)
		return new_chain

	new_chain=chr(ord('A')-1)
	new_chain = encoder_individual(new_chain, pdb_rec)
	encoder_individual(new_chain, pdb_lig)


	with open("pre_map.pickle", "wb") as f:
		pickle.dump(mapping, f)

def decoder():

	pdb_rec = sys.argv[1]
	pdb_lig = sys.argv[2]

	def decoder_individual(path):
		with open("pre_map.pickle", "rb") as f:
			mapping=pickle.load(f)	

		with open("temp", "w") as fo:
			with open(path, "r") as f:
				for lines in f:
					if lines[0:4]=='ATOM':
						fo.write(lines[0:21]+mapping[lines[21:27]]+lines[27:])
					else:
						fo.write(lines)

		os.system("mv temp "+path)

	decoder_individual(pdb_rec)
	decoder_individual(pdb_lig)

if sys.argv[3]=='in':
 encoder()
else:
 decoder()
