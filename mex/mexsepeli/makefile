# svn $Id$

linux:
	mex -f ./mexopts.sh -v -o mexsepeli mexgateway_sepeli.c sepeli.F

sol:
	mex -v -o mexsepeli mexgateway_sepeli.F mexsepeli.F

clean:
	rm -f *.o *.f mexsepeli.mex*
