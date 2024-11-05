#!/bin/bash

YEAR=
NJOBS=
DIR=root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/cms_offline/mini_aod/svjtau/UL17_SVJtau_private_prod_scan/
USER=$(whoami)
case $USER in
	pedrok)
		YEAR=2016
		NJOBS=46
	;;
	bregnery)
		YEAR=2017
		NJOBS=52
	;;
	mekhpar)
		YEAR=2018
		NJOBS=75
	;;
	cazzanig)
		YEAR=2017
		NJOBS=75
	;;
esac

if [ -z "$YEAR" ] || [ -z "$NJOBS" ]; then
	echo "Unknown user $USER"
	exit 1
fi

(set -x;
python runProd.py -P P8 -G="-p -d signals_P8_ex_svjtau --svjl 1 -E 2000 -N 80 --cpus 4 --memory 8000 --args 'quiet=1 syst=1'" -y 2017 -n chain2017_ -o root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/cms_offline/mini_aod/svjtau/UL17_SVJtau_private_prod_scan/ -t root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/cms_offline/mini_aod/svjtau/UL17_SVJtau_private_prod_scan/ -c -s
)
