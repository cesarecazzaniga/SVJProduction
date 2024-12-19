flist = []

mPseudo = 8.0
mVector = 15
for MZp in [3000]:
        for rinv in [0.3]:
          for brgamma in [0.3]:       
            for mPiOverLambda in [1.0]:     
                 for lambdaHV in [2.0]: 
                        flist.append({"channel": "s", "svjgamma": 1 , "mMediator": MZp, "mPseudo": mPseudo,"mVector": mVector,"rinv": rinv, "alpha": "peak", "lambdaHV": lambdaHV, "mPiOverLambda": mPiOverLambda, "BRGamma": brgamma })
