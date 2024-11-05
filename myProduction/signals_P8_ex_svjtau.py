flist = []

mPseudo = 8.0
mVector = 15
#mPioverLambda 1.0, 1.2

for MZp in [5000]: 
    for rinv in [0.1, 0.2, 0.3, 0.5, 0.7]:          
        for lambdaHV in [10.0, 20.0, 50.0]: 
            for mPiOverLambda in [0.8,1.0, 1.2]:
                for BRtau in [0.1,0.3,0.5, 0.7, 0.8]:    
                    flist.append({"channel": "s", "svjl": 1 , "mMediator": MZp, "mPseudo": mPseudo,"mVector": mVector,"rinv": rinv, "alpha": "peak", "lambdaHV": lambdaHV, "mPiOverLambda": mPiOverLambda, "BRtau": BRtau})

