## Makefile.am -- Process files to build Makefile for Gibbs
bin_PROGRAMS = Gibbs

## CFLAGS = -g -D_CYGWIN_ -O3 -ffloat-store -malign-double -Wall -Werror 

EXTRA_DIST=libml.a gibbsmpi.c gibbsmpi.h hier.c hier.h digammal.c digammal.h gpl.license.txt NonGUI.nompi.inc make.default hum.mouse.pr Gibbs.linux seq.weights.pl calc.nmer.prob.pl

if CYGWIN
    LDADD = libml.a -lm
    Gibbs_CFLAGS = -malign-double -D_CYGWIN_
else
    LDADD = -lm
    Gibbs_CFLAGS = 
endif

if MPI
   SRC = gibbsmpi.c gibbsmpi.h hier.c hier.h digammal.c digammal.h
   MPILDADD = -lmpe
else
   SRC = 
   MPILDADD =
endif

if SUNMATH
   SUNLDADD = -lsunmath
else
   SUNLDADD =
endif

Gibbs_LDADD = ${LDADD} ${MPILDADD} ${SUNLDADD}

Gibbs_SOURCES = Gibbs.c Gibbs.h Gibbs_NON_GUI.c Gibbs_NON_GUI.h calcstats.c calcstats.h centroid.c centroid.h column.c column.h common.c common.h counts.c counts.h expectmax.c expectmax.h fragadj.c fragadj.h fragment.c fragment.h freqbkgnd.c freqbkgnd.h initseq.c initseq.h karlin.c karlin.h maximize.c maximize.h mem_mgmt.c mem_mgmt.h motif.c motif.h motifsamp.c motifsamp.h nearopt.c nearopt.h parameters.c parameters.h printdata.c printdata.h probability.c probability.h process.c process.h sampling.c sampling.h sitesamp.c sitesamp.h subopt.c subopt.h width.c width.h wilcoxon.c wilcoxon.h perror.c perror.h recur.c recur.h rsample.c rsample.h prior.c prior.h motifpair.c motifpair.h footprint.c footprint.h exchange.c exchange.h groupsampler.c groupsampler.h weights.c weights.h trans.c trans.h position.c position.h seqmap.c seqmap.h phyllip.c phyllip.h fels.c fels.h impsamp.c impsamp.h rgamma.c rgamma.h model.c model.h bayes.c bayes.h blosum.h species_sample.c species_sample.h ${SRC}




