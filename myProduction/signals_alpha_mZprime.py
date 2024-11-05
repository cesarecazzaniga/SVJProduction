flist = []

alpha = "peak"
m = 20
for z in [800 , 900, 1000, 1500]:
    for rinv in [0.3,0.5,0.7]:
        flist.append({"channel": "s", "mMediator": z, "mDark": m, "rinv": rinv, "alpha": alpha, "scout": True})
#flist.append({"channel": "s", "svjl": 0, "mMediator": 700, "mDark": m, "rinv": 0.3, "alpha": "peak", "scout": True})

