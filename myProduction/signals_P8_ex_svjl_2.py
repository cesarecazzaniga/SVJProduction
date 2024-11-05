flist = []

mPseudo = 8.0
mVector = 15
mPiOverLambda = 1.6
lambdaHV=5.0
for MZp in [700 , 800, 1000, 1500]:    
        for rinv in [0.3, 0.5 , 0.7]:   
                flist.append({"channel": "s", "svjl": 1 , "mMediator": MZp, "mPseudo": mPseudo,"mVector": mVector,"rinv": rinv, "alpha": "peak", "lambdaHV": lambdaHV, "mPiOverLambda": mPiOverLambda, "scout": True })

