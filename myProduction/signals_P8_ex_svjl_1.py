flist = []

mPseudo = 8.0
mVector = 15
#mPiOverLambda = 1.6
#lambdaHV=5.0
for MZp in [1500, 2000, 3000, 4000, 5000]:
        for rinv in [0.1, 0.2, 0.3, 0.5, 0.5,0.7]:
          for mPiOverLambda in [1.6, 2.0, 3.0]:     
                for lambdaHV in [5, 10, 20]: 
                        flist.append({"channel": "s", "svjl": 1 , "mMediator": MZp, "mPseudo": mPseudo,"mVector": mVector,"rinv": rinv, "alpha": "peak", "lambdaHV": lambdaHV, "mPiOverLambda": mPiOverLambda })
