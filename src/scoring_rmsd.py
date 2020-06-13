from sklearn.ensemble import RandomForestRegressor
import sys
import numpy as np
import pickle
sample=[]

# 'random_forest.sav' is the random forest version.  'ridge_rbf.sav' is the ridge rbf version
if(sys.argv[10][-5]=='t'):

	for i in range(8):
    		sample.append(((float(sys.argv[i+1]))))
	print sample

	rf = pickle.load(open(sys.argv[10], 'rb'))
	np.savetxt("Result/scores", rf.predict(np.asarray(sample).reshape(1,-1)))


else:
	
	mean = [3.016821065854138340e+01,
2.333334821112641322e+01,
-3.005433855907200780e+01,
2.571341663947315226e+00,
7.150039753292698030e+00,
7.861239728720292419e+00,
-1.035823247493610744e+01,
-1.653124477884800854e+03
]
	scale =[2.333590815035907156e+02,
7.921142666067396476e+01,
4.433288510052415177e+01,
3.303284207555231333e+00,
4.028303297576421471e+01,
1.089866893350210262e+02,
3.298226217857801235e+00,
1.099759536416193214e+03,
]

	for i in range(8):
        	sample.append(((float(sys.argv[i+1]))-mean[i])/scale[i])
	print sample

	rf = pickle.load(open(sys.argv[10], 'rb'))
	np.savetxt("Result/scores", rf.predict(np.asarray(sample).reshape(1,-1)))	
