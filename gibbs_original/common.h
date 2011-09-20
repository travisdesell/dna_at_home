/****************************************************************/
/* Gibbs - A program for detecting subtle sequence signals      */
/* */
/* Please acknowledge the program authors on any publication of */
/* scientific results based in part on use of the program and   */
/* cite the following articles in which the program was         */
/* described.                                                   */
/* For data involving protein sequences,                        */
/* Detecting subtle sequence signals: A Gibbs sampling          */
/* strategy for multiple alignment. Lawrence, C. Altschul,      */
/* S. Boguski, M. Liu, J. Neuwald, A. and Wootton, J.           */
/* Science, 262:208-214, 1993.                                  */
/* and                                                          */
/* Bayesian models for multiple local sequence alignment        */
/* and Gibbs sampling strategies, Liu, JS. Neuwald, AF. and     */
/* Lawrence, CE. J. Amer Stat. Assoc. 90:1156-1170, 1995.       */
/* For data involving nucleotide sequences,                     */
/* Gibbs Recursive Sampler: finding transcription factor        */
/* binding sites. W. Thompson, E. C. Rouchka and                */
/* C. E. Lawrence, Nucleic Acids Research, 2003,                */
/* Vol. 31, No. 13 3580-3585.                                   */
/* */
/* Copyright (C) 2006   Health Research Inc.                    */
/* HEALTH RESEARCH INCORPORATED (HRI),                          */
/* ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.             */
/* Email:  gibbsamp@wadsworth.org                               */
/* */
/****************************************************************/
/* */
/* Changes Copyright (C) 2007   Brown University                */
/* Brown University                                             */
/* Providence, RI 02912                                         */
/* Email:  gibbs@brown.edu                                      */
/* */
/* For the Centroid sampler, please site,                       */
/* Thompson, W.A., Newberg, L., Conlan, S.P., McCue, L.A. and   */
/* Lawrence, C.E. (2007) The Gibbs Centroid Sampler             */
/* Nucl. Acids Res., doi:10.1093/nar/gkm265                     */
/* */
/* For the Phylogenetic Gibbs Sampler, please site,             */
/* Newberg, L., Thompson, W.A., Conlan, S.P., Smith, T.M.,      */
/* McCue, L.A. and Lawrence, C.E. (2007) A phylogenetic Gibbs   */
/* sampler that yields centroid solutions for cis regulatory    */
/* site prediction., Bioinformatics,                            */
/* doi:10.1093/bioinformatics/btm241.                           */
/* */
/****************************************************************/
/* */
/* This program is free software; you can redistribute it       */
/* and/or modify it under the terms of the GNU General Public   */
/* License as published by the Free Software Foundation;        */
/* either version 2 of the License, or (at your option)         */
/* any later version.                                           */
/* */
/* This program is distributed in the hope that it will be      */
/* useful, but WITHOUT ANY WARRANTY; without even the implied   */
/* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR      */
/* PURPOSE. See the GNU General Public License for more         */
/* details.                                                     */
/* */
/* You should have received a copy of the GNU General Public    */
/* License along with this program; if not, write to the        */
/* Free Software Foundation, Inc., 51 Franklin Street,          */
/* Fifth Floor, Boston, MA  02110-1301, USA.                    */
/****************************************************************/
/**************************************************************************/
/* $Id: common.h,v 1.25 2008/06/24 22:16:38 Bill Exp $            */
/* */
/* Author:       Eric C. Rouchka,1996,7.                                  */
/* Jun Zhu, 1996,8.                                         */
/* Bill Thompson 1/27/97				          */
/**************************************************************************/

#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <limits.h>
#include <sys/stat.h>
#include "perror.h"

/* ======================================================================== */
/* =================================CONSTANTS============================== */
/* ======================================================================== */
#define DEBUG 1
/* #define MAX_WIDTH_MULT 10 */
#define MAX_WIDTH_MULT 1.5
#define PLATEAU_PERIODS 20
#define NUM_SEEDS 10
#define MAX_ITERATIONS 500
#define PSEUDO_CNT_WT 0.1
#define PSEUDO_SITE_WT 0.8
#define NEIGHBOR_RESIDUES 5
#define NEAROPT_CUT 0.5
#define NEAROPT_CUTOFF_PERCNT 0.01
#define MIN_PSEUDOCNT 0.001
#define FALSE 0
#define TRUE 1
#define WT_TABLE_ID_SIZE 32

#define MIN_MOTIF_LEN    5	/* BT 5/7/97 */
#define MOTIF_SAMPLE_CNT 100	/* BT 5/7/97 - temporary */
#define RCUTOFF_FACTOR   0

#define MAX_ARG_SIZE    1024

#define EXCHANGE_ITERATIONS   20
#define MINTEMP               1.0
#define MAXTEMP               2.0

#ifndef INT32_MAX
#define INT32_MAX  2147483647
#endif

/* ======================================================================== */
/* ==================================MACROS================================ */
/* ======================================================================== */

#define max(a, b)	((a) > (b) ? (a) : (b))
#define min(a, b)	((a) < (b) ? (a) : (b))

/* BT 1/10/2000 */

#define MEW(x,n,t)      (( (x=(t*) malloc(((n)*sizeof(t))))==NULL) ? \
                         (t*) (ErrorMsg( "Out of memory", __FILE__, __LINE__),exit(1),(int *) 0):x)
#define NEW(x,n,t)      (( (x=(t*) calloc(max(n,1),sizeof(t)))==NULL) ? \
                         (t*) (ErrorMsg( "Out of memory", __FILE__, __LINE__),exit(1),(int *) 0):x)
#define NEWP(x,n,t)     (( (x=(t**) calloc(max(n,1),sizeof(t*)))==NULL) ? \
                        (t**) (ErrorMsg( "Out of memory", __FILE__, __LINE__),exit(1),(int *) 0):x)

#define NEWPP(x,n,t)    (( (x=(t***) calloc(max(n,1),sizeof(t**)))==NULL) ? \
                        (t***) (ErrorMsg( "Out of memory", __FILE__, __LINE__),exit(1),(int *) 0):x)

#define NEWP3(x,n,t)    (( (x=(t****) calloc(max(n,1),sizeof(t***)))==NULL) ? \
                        (t****) (ErrorMsg( "Out of memory", __FILE__, __LINE__),exit(1),(int *) 0):x)


/* Note for FREEP need local variable i defined  */
/* and for FREEPP need variables i and j defined */

/*
 * Note: It's too easy to get in trouble with local variables i,j. It's too
 * common to use them for loop counters and then include FREEP in the loop.
 * See motifsamp.c for an example. Instead we will define GLOBAL variables
 * GL_I__, GL_J__ for these macros.
 * 
 * DO NOT declare them anywhere but in common.c
 * 
 * BT 1/21/97
 */

/* BT 1/21/97  					 */

extern short    GL_GUI;		/* BT 6/17/97 */
extern int      GL_I__;
extern int      GL_J__;
extern int      GL_Finished;

extern int      glVerMajor;	/* BT 3/28/97 */
extern int      glVerMinor;
extern int      glRev;
extern char     glRevDate[];
extern int      nComp[];

extern int      GUI_argc;
extern char   **GUI_argv;


#define FREEP(x, n) 	{for(GL_I__=0;GL_I__ < (n);GL_I__++) \
                           if((x)[GL_I__] != NULL) free((x)[GL_I__]);\
                        if((x) != NULL) free((x));(x) = NULL;}

#define FREEPP(x,ii,jj)	{for(GL_I__=0;GL_I__ <(ii);GL_I__++){\
			   for(GL_J__=0;GL_J__<(jj);GL_J__++) \
                             free((x)[GL_I__][GL_J__]);\
                           free((x)[GL_I__]);}\
                           if((x) != NULL) free((x));(x)=NULL;}



#define NUMMOTIFS(x)  ((x)[FORWARD] + (x)[REVERSE])

#define CH2INT(n,B)	((int)(*(B)->Seq->R)[(n)] - 97)	/* BT 2/21/97 */

#define CH2INTCOMP(n,B)	((int)(complement((*(B)->Seq->R)[(n)])) - 97)	/* BT 2/21/97 */

/*
#define GETSTRING(S) i = 0; Finished = FALSE;\
                     while(!Finished) {\
                        scanf("%c", &(S)[i]);\
                        if((S)[i] == '\n') {\
                           Finished = TRUE;\
                           (S)[i] = '\0';\
                        }\
                        i++;\
                     }
*/

/* BT 1/27/97 */

#define GETSTRING(S) GL_I__ = 0; GL_Finished = FALSE;\
                     while(!GL_Finished) {\
                        scanf("%c", &(S)[GL_I__]);\
                        if((S)[GL_I__] == '\n') {\
                           GL_Finished = TRUE;\
                           (S)[GL_I__] = '\0';\
                        }\
                        GL_I__++;\
                     }



/* ================== define enumerate constants ======================= */

enum flags {
	DELETE, ADD
};
enum where {
	BG, MOTIF
};
enum indicator {
	POSSIBLE, ENDS, MARKED
};
enum direction {
	FORWARD, REVERSE
};
enum col_frag {
	COL_OFF, COL_ON, COL_BLOCKED, COL_DONT_USE
};
enum cl_params {
	cl_A, cl_B, cl_C, cl_D, cl_E, cl_F, cl_G, cl_H, cl_I, cl_J,
	cl_K, cl_L, cl_M, cl_N, cl_O, cl_P, cl_Q, cl_R, cl_S, cl_T,
	cl_U, cl_V, cl_W, cl_X, cl_Y, cl_Z,
	cl_a, cl_b, cl_c, cl_d, cl_e, cl_f, cl_g, cl_h, cl_i, cl_j,
	cl_k, cl_l, cl_m, cl_o, cl_p, cl_q, cl_r, cl_s, cl_t, cl_u,
	cl_v, cl_w, cl_x, cl_y, cl_z,
	cl_frag, cl_nopt, cl_nopt_disp, cl_sample_model, cl_hm,
	cl_wilcox, cl_bayes, cl_opt, cl_freq, cl_bayes_counts,
	cl_pred_update, cl_max, cl_align_centroid, cl_no_cred, cl_freq_background
};
enum print_types {
	SUBOPT, NEAROPT, MAPMAX, CENTROID, OTHER
};				/* for print_info to print control printing
				 * of cl_N option */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* cl_a : Concentrated Region           cl_t : Sequence type             */
/* cl_c : Collapsed Alphabet            cl_w : Pseudocount weight        */
/* cl_e : Use Expectation/Maximization  cl_x : Don't Xnu sequence        */
/* cl_h : Help flag                     cl_C : Near optimal cutoff       */
/* cl_i : Number of iterations          cl_F : Don't fragment            */
/* cl_m : Don't use map maximization    cl_I : Interactive mode          */
/* cl_o : Output file                   cl_P : Informed priors file      */
/* cl_p : Plateau periods               cl_R : palandromic sequence      */
/* cl_r : Reverse complement            cl_S : Number of seeds           */
/* cl_s : Seed Value                    cl_W : Pseudosite weight         */
/* cl_u : suboptimal output             cl_v : overlap value             */
/* cl_d : Allow width to vary           cl_l : Wilcoxon                  */
/* cl_a : Sample along length           cl_N : Output Scan file          */
/* cl_O : Output prior file             cl_E : Use recursive sampling    */
/* cl_f : recursion cutoff              cl_q : post. width sampling      */
/* cl_Q : Save probability samples      cl_b : sample background         */
/* cl_z : Ignore spacing model          cl_k : post plateau iter         */
/* cl_B : use bkgnd comp model          cl_A : init sample from prior    */
/* cl_T : use footprint data            cl_X : exchange Monte Carlo      */
/* cl_M : frag width multiplier         cl_Z : no progress info          */
/* cl_U : spacing model                 cl_J : fragment only middle      */
/* cl_V : verify mode                   cl_K : alternate sample on k     */
/* cl_y : don't print freq. soln        cl_Y : calc. def. pseud wt       */
/* cl_L : motif sample befor recursive                                   */
/* cl_phy : phylogenic sampling         cl_frag : old style frags        */
/* cl_nopt : supress near opt soln.     cl_sample_model : sample mod     */
/* cl_hm : hierarchical model           cl_wilcox: call wilcoxon         */
/* cl_bayes : Bayesian sampling         cl_opt : print near opt soln     */
/* cl_freq : print freq soln            cl_bayes_counts: counts in output */
/* cl_pred_update: use pred update      cl_max : Map Maximization        */
/* cl_align_centroid: alingn cent model cl_no_cred : no calc. cred limit */
/* cl_freq_background: use n-mer freqs                                   */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ======================================================================= */
/* ================================DATA TYPES============================= */
/* ======================================================================= */

typedef struct M1 {
	int             Mtype;
	int             pos;
	int             seq_num;
	int             left, right;
	char           *sMotif;
	short           RevComp;
	struct M1      *next;
}               MlistEl;

typedef struct {
	int             nNumMotifs;
	int             nMotifLen;
	MlistEl        *Motifs;
}               Mlist_struct;

typedef Mlist_struct **Mlist;

typedef struct {
	int             nPossStartPos;
	int             nInMotif;
	int             nMotifStartPos;
	short           RevComp;/* Is Motif Reverse complement? */
	short           bEndOfSeq;	/* Is it the end of a sequence? - BT
					 * 4/18/97 */
	int             nSeq;	/* Sequence number - BT 4/21/97 */
}               PoSition;

typedef struct {
	double          dResidueInBGProb;
	double          dResidueInMotifProb;
	double          dinBGProb;
	double          dinMotifProb;
	double          dBGProbDenom;
	double          dMotifProbFact;
	double       ***dvInMotifProb;	/* For Each individual residue */
	double        **dvInBGProb;
	unsigned short  update;
}               ProbStruct;

typedef struct {
	double       ***fCounts;/* Observed count      *//* BT 4/18/2000 */
	/* fCounts[Motiftype]  */
	/* [PositioninMotif] */
	/* [Alphabet]    */
	double       ***dPseudoCounts;	/* Pseudo Counts        */
	/* dPseudoCounts[Motiftype]  */
	/* [PositioninMotif] */
	/* [Alphabet]    */
	double       ***dSumPseudo;	/* sum pseudo counts    */
	/* dSumPseudo[Motiftype] */
	/* [Background/Motif]  */
	/* [Alphabet]    */
	/* arrange memory in this way */
	/* is not very efficient, but */
	/* it is clear when used with */
	/* above arrays               */

	double          dSumBGPseudo;
	int             nTotBack;	/* Total Observed BG Residues    */
	int             nTotMotif;	/* Total Observed Motif Residues */
	double          dTotSumMotifPseudo;

	double        **dmodel_sites;	/* See protein science paper */
	double        **dmodel_pseudo;	/* Appendix and p 1627 (3)   */
	double         *dtot_sites;	/* for the description       */
	double         *dtot_pseudo;
	double         *dbg_pseudo;
	double         *dTot;

	int             nMaskCount;	/* count of mask characters - BT
					 * 4/16/97 */
	double          dBackCnt;	/* actual backgrnd count for weights */
	double       ***wCounts;/* Observed count      *//* BT 4/18/2000 */
	/* fCounts[Motiftype]  */
	/* [PositioninMotif] */
	/* [Alphabet]    */
}               Counts;

typedef Counts *Ctype;

typedef struct {
	FILE           *fpt;	/* sequence file pointer    */
	FILE           *prior_fpt;	/* beta priors file pointer */
	FILE           *out_fpt;/* output file pointer      */
	FILE           *sn_fpt;	/* BT 9/17/97 */
	FILE           *prout_fpt;	/* prior output file *//* BT 9/26/97 */
	FILE           *control_fpt;	/* Control file for Wilcoxon */
	FILE           *sample_fpt;	/* Sample count file */
	FILE           *occur_fpt;	/* occurence count file */
	FILE           *near_occur_fpt;	/* near opt occurence count file */
	FILE           *sites_fpt;	/* sites distribution file */
	FILE           *sitespr_fpt;	/* sites prior distribution file */
	FILE           *dist_fpt;	/* distance distribution file */
	FILE           *prob_fpt;	/* alignment prob */
	FILE           *bkgnd_fpt;	/* background composition */
	FILE           *spacing_fpt;	/* spacing file     */
	FILE           *weight_fpt;	/* sequence weights */
	FILE           *mpiTemp_fpt;	/* temp output for MPI */
	FILE           *mpiSample_fpt;	/* temp output for MPI centroid
					 * samples */
	FILE           *hier_fpt;	/* temp file for hierarchical model */

	char           *filename;
	char           *prior_filename;
	char           *output_filename;
	char           *InputFilename;	/* BT 7/2/97 */
	char           *WilcoxonFilename;	/* BT 7/2/97 */
	char           *ScanFilename;	/* BT 9/17/97 */
	char           *PriorOutFilename;	/* BT 9/26/97 */
	char           *ControlFileName;	/* BT 1/28/98 */
	char           *SampleFileName;	/* BT 3/25/98 */
	char           *OccurFileName;	/* BT 4/3/98 */
	char           *NearOccurFileName;	/* BT 4/3/98 */
	char           *SitesFileName;
	char           *SitesPrFileName;
	char           *DistFileName;
	char           *ProbFileName;
	char           *BkgndFileName;
	char           *SpacingFileName;
	char           *WeightFileName;
	char           *mpiTempFileName;
	char           *mpiSampleFileName;
	char           *hierTempFileName;
}               files;

typedef struct {
	int            *SeqLen;	/* length of each sequence  */
	short         **Orig;	/* sequence in number       */
	char          **R;	/* all strings concatenated */
	char          **ProcessedSTR;	/* sequence after Xnu       */
	int           **nvEndLocs;	/* end location of each sequence in
					 * **R                  */
}               Stringstruct;

typedef Stringstruct *Stype;

typedef struct {
	short         **Concentrated;
	short          *NumConcen;
}               Concenstruct;

typedef Concenstruct *ConType;

typedef struct {
	int           **Collapsed;	/* Array indicating positions that
					 * are collapsed   */
	int           **Palandromic;	/* Array indicating positions that
					 * are palandromic */
	int           **Repeat;	/* Array indicating positions that are direct
				 * repeats */
	ConType         Concen;
}               AltModelStruct;

typedef AltModelStruct *AltModelType;


typedef struct {
	short           is_defined[66];
	short           site_samp;
	int            *nMotifLen;	/* Motif Residue Lengths       */
	int             nSeqLen;/* Sequence Residue Length     */
	int             nAlphaLen;	/* # of characters in alphabet */
	int           **nNumMotifs;	/* number of initial motifs    */
	int             nNumSequences;	/* number of different seq.    */
	files          *Datafiles;	/* datafiles structure         */
	double          dPseudoCntWt;
	double          dPseudoSiteWt;
	int             nMaxIterations;
	long            lSeedVal;
	int             nPlateauPeriods;
	int             nSeeds;
	double         *dposterior_prob;
	int            *nPossSites;
	double          dnull_map;
	int            *col_shift;
	int             nNumMotifTypes;
	double          dCutoff;
	double          dNearOptCutoff;
	short           RevComplement;
	int             nNumProcessed;
	int            *DOF;	/* degree of freedom for each motif type                    */
	AltModelType    AltModel;
	double          Logl;
	double          glOverlapParam;	/* % of overlap at ends - BT 4/21/97 */

	/* Values for Metropolis Length Sampler */
	int             nMotifSampleCnt;	/* BT 5/7/97 */
	int            *nMinMotifLen;
	int            *nMaxMotifLen;
	int            *nInputMotifLen;	/* Original Motif Residue Lengths       */
	double       ***nWidthCnts;	/* Histogram for post. dist. of width
					 * sampling */
	int             nOrigSeqCount;	/* BT 7/2/97 *//* For Wilcoxon calc. */
	int             nSeqAllocCnt;	/* BT 8/27/97 */

	int             nMaxBlocks;	/* BT 2/11/98 *//* max blocks for
					 * recursion */
	int             nMinBlocks;	/* min sites/seq  */
	double          dRCutoff;
	int             nSampleCnt;	/* sample count for width dist. */

	int         ****nAlignCnts;	/* sample counts for alignment */
	int             nMaxSeqLen;	/* max sequence length */
	int             nRank;	/* process num for multiprocessing */
	int             nAdjustPeriod;	/* Fragment/Shift Period */
	int             nPostPlateauIter;	/* Post plateau sample
						 * iterations */
	int            *nMaxFragWidth;	/* frag widths */
	int             nMPI;	/* using MPI ? */
	int             nMPIProcesses;	/* # of MPI processes */
	int             nSeedRun;	/* current seed */
	int             nVerifySeq;	/* number of sequences to verify */
	char           *comment;/* optional comment printed at start of
				 * output */
	long            lSeedVal0;	/* MPI - seed for node zero */
	double          dKSampleMap;	/* temp storage for dKSampleMap */
	char          **fastaHeader;	/* storage for the fasta headers */
	int             hier_iter;	/* number of iterations between
					 * exchange for hier model */
	int             burnInPeriod;	/* number of burnin iterations for
					 * Bayesian sampling */
	int             bayesSamplePeriod;	/* number of sampling
						 * iterations for Bayesian
						 * sampling */
	int             bayesSampleIter;	/* sampling iteration , -1
						 * for burn-in */
	int             inCentroidAlign;	/* are we aligning the
						 * centroid solution? */
	int             argc;	/* initial arguments */
	char          **argv;
	char           *programName;	/* for align_centroid */
	double        **freqSeqs;	/* for n-mer background */
	int            *nMerSize;	/* size of n-mer for each sequence */
	long            ticks;	/* clock ticks per second */
}               InputParams;

typedef InputParams *IPtype;


typedef struct PHYLONODE {
	char           *species;/* species name */
	int             offset;	/* offset of species sequence in groups of
				 * sequences */
	double          length;	/* distance to parent */
	double         *likelihood;	/* vector of probabilities for this
					 * node */
	double        **subs;	/* the substitution matrix */
	struct PHYLONODE *leftChild;
	struct PHYLONODE *rightChild;
}               PhyloNode;

typedef PhyloNode *PhyloTree;

typedef struct {
	int             treeCount;	/* number of trees */
	PhyloTree      *phyloTree;	/* phylogenic tree structure array
					 * for Felsenstein algorithm */
	int            *phyloSpecies;	/* no of aligned species for each
					 * tree */
	int            *phyloSeq;	/* number of sequences covered by
					 * each tree */
	int            *phyloIndex;	/* which tree describes each sequence */
	int             bCalcPhylo;	/* use Fels algorithm for map */
	int             phyloSpeciesCount;	/* count of SeqInc */
	double        **phyloBkgnd;	/* background Felsenstein sample */
	double          phyloNull;	/* sum of logs of phyloBkgnd */
	int             maxPhyloSeq;	/* max sequence to apply the trees to */
	int             phyloAllocCnt;	/* allocation sise of phyloTree array */
	int            *phyloSpeciesSample;	/* should tree use species
						 * sampling */
	int            *phyloTreeStartSeq;	/* starting sequence
						 * described by each tree */
}               PhyloStruct;

typedef PhyloStruct *PhyloType;


typedef struct {
	int            *nMaxLen;/* 5 * width of motif         */
	int           **nColMask;	/* Column Mask for each motif */
	int           **nOldColMask;
	int            *FragWidth;
	int            *nOldFragWidth;	/* BT 02/05/2001 */
	double       ***nvFragCnts;
	int            *shift;	/* how far first column has shifted */
	int           **fragPos;/* distance from 1st pos */
	int           **fragInitMask;	/* initial frag mask */
}               FragStruct;

typedef FragStruct *Ftype;


typedef struct {
	int             nMotifStart;
	int             nRevComp;
	double          dSpacingProb;
	double          dPartialSum;
	double          dSpacingProbSave;
	double          dPartialSumSave;
	double          addVal;
}               SiteStruct;

typedef SiteStruct *SiteType;


typedef struct {
	int             pos;
	int             width;
	int             t;
}               SampleStruct;

typedef SampleStruct *SampleType;


typedef struct {
	int             nMaxBlocks;	/* Maximum number of blocks per
					 * sequence */
	int             nMinBlocks;	/* Minimum number of blocks per
					 * sequence */
	int             nMaxGap;/* max size of Gap probaility array */
	int             nSeq;	/* The current sequence */
	int             nSeqLen;/* length of sequence nSeq  */
	int             nMaxSeqLen;
	double          dPriorSites;	/* alpha - prior on number of sites */
	double          dPriorSeq;	/* beta - prior on number of seqs */
	double          dSiteProbCons;	/* Prob of a site given that it is
					 * conserved */
	double          dSiteProbNonCons;	/* Prob of a site given that
						 * it is not conserved */
	double          dSiteProbWt;	/* Amount to weight pseudocnts for
					 * site prob */
	int           **nMaxGapWidth;	/* Max Gap, if > 0 then beyond this
					 * width, pr = 0 */
	double          prTotalSum;	/* Total sum of pr matrix */
	double         *prSum;	/* denominator for probability calc */
	double      ****dPTrans;/* Probability of transition between motif
				 * types */
	double       ***dPGap;	/* Gap probability */
	double        **dPEndGap;	/* Gap probability to end of sequence */
	double        **dPSpacing;	/* Probability of a site a loc i */
	double      ****dProb;	/* P(R | t, k, dir )   (unnormalized) */
	double         *dExpBlkCnt;	/* Expected number of sites/sequence */
	double        **dProbCons;	/* Prob that site is conserved w
					 * similar species */
	double        **dProbSite;	/* Prob that site starts at pos,
					 * given footprint */
	double        **dProbNonSite;	/* Prob that site is bkg at pos,
					 * given footprint */
	int            *nSites;	/* sites in each seq */
	double         *alphaCons;	/* Parameters for footprint
					 * calculations */
	double         *alphaNonCons;
	double         *betaCons;
	double         *betaNonCons;
	int             nUseFixedBlockProbs;	/* use fixed probs rather
						 * than neg binomial */
	int             nUseSpacingProb;	/* use spacing model,
						 * dPSpacing, for sites */
	SiteStruct   ***sitePos;
	int             nGroupCnt;	/* groups for block sampler */
	int            *groups;
	int             nScaleZeroBlocksProb;	/* scale fixed block prob ? */
	int            *nSeqPos;/* get the seq from position */
	int            *nInMotif;	/* is position in Motif? */
	int             nonCutoffSites;	/* sites >= FOOTPRINT_CUTOFF */
	double        **dPFootprint;	/* normalized footprint prob for
					 * current seq */
	int             nUseGap;/* Is there an explict spacing prob? */
	int             bUseFlatGap;	/* spacing prob is flat an dthe same
					 * for all motif types */
	int             nFlatGapLen;	/* length of flat spacing widnow */
	double          dFlatGapProb;	/* constant gap prob - override other
					 * gap probs */
	double         *dEndSiteProb;	/* prob of sampling each site type as
					 * kmax site */
	double         *dInitialEndSiteProb;	/* Initial prob of sampling
						 * each site type as kmax
						 * site */
	double         *dEndCnts;	/* pseudocnts for each site type as
					 * kmax site */
	double         *dPriorEndCnts;	/* Initial prior counts for each site
					 * type as kmax site */
	double          dTotalPriorEndCnts;	/* sum of the prior cnts for
						 * each site type as kmax
						 * site */
	double         *dBeginSiteProb;	/* prob of sampling each site type as
					 * site # 1 */
	double         *dInitialBeginSiteProb;	/* Initial prob of sampling
						 * each site type as 1st site */
	double         *dBeginCnts;	/* pseudocnts for each site type as
					 * 1st site */
	double         *dPriorBeginCnts;	/* Initial prior counts for
						 * each site type as 1st site */
	double          dTotalPriorBeginCnts;	/* sum of the prior cnts for
						 * each site type as 1st site */
	double      ****dInitialPTrans;	/* Initial Probability of transition
					 * between motif types */
	double        **dPriorTransCnts;	/* Prior on the number of
						 * transitions between motif
						 * types */
	double        **dTransCnts;	/* Updated number of transitions
					 * between motif types */
	int             bUseTrans;	/* boolean indicating whether trans
					 * matrix is used */
	int             bUpdateTrans;	/* shoud trans matrix be updated? */
	int             bSymTrans;	/* shoud trans matrix be symmetric? */
	double          dTransWt;	/* weight factor for trans
					 * pseudocounts */
	double          dTransTotalPriorCnts;	/* the sum of dPriorTransCnts */
	double          dTransTotalCnts;	/* the sum of dTransCnts */
	int             bUsePosMatrix;	/* use the position matrix method for
					 * finding modules */
	double        **dInitialPosMatrix;	/* Initial probability of
						 * each motif model in each
						 * pos of module */
	double        **dPriorPosMatrixCnts;	/* prior on number of each
						 * motif model in each pos of
						 * module */
	double        **dPosMatrixCnts;	/* updated number of each motif model
					 * in each pos of module */
	double        **dPosMatrix;	/* Probability of each motif model in
					 * each pos of module */
	int             nModuleSites;	/* number of sites in module */
	int             bUpdatePosMatrix;	/* should the matrix be
						 * updated? */
	double          dPosWt;	/* weight for pos pseudocounts */
	double          dPosTotalPriorCnts;	/* the sum of dPriorPosCnts */
	double          dPosTotalCnts;	/* the sum of dPosCnts */
	double          dExpBlkWt;	/* weight on Expected number of
					 * sites/sequence */
	double          dMinSiteMap;	/* min map value for including P(A)
					 * in modular sampling */
	int             nIncludeAlign;	/* should we include P(A) in modular
					 * sampling */
	int             nMinSiteMotifs;	/* min number of motifs before
					 * incuding P(A) in modular sampling */
	double          dInitProb;	/* initial MAP value - set if
					 * B->InitPos != NULL */
	int             nDontAdjustPrior;	/* don't recalculate prior on
						 * k */
	int           **nPriorMotifSites;	/* initial estimate of number
						 * of sites for each motif
						 * type */
	double          dKSampleMap;	/* min map before sampling on all k */
	double         *priorSitesPerSeq;	/* pseudoCounts for sites/seq */
	double         *probK;	/* probability of sites/seq */
	int             useProbK;	/* flag - set by ExchangeK */
	double          sitesPerSeqNullMap;	/* site per seq portion of
						 * Null MAP */
	int           **kCounts;/* keep track of the counts of the number of
				 * times each k is sampled */
	int          ***totalPosCounts;	/* counts for each position for
					 * Bayesian sampling */
	int           **kTotalCounts;	/* counts of the number of times each
					 * k is sampled */
	SampleStruct ****samples;	/* samples for credibility limits */
	int             acceptCount;	/* acceptance count for M-H test in
					 * model.c */
	int             bayesSampleCount;	/* total samples for M-H test
						 * in model.c */
}               RProbStruct;

typedef RProbStruct *RPType;	/* BT 12/12/97 */


typedef struct {
	int             nMaxBlocks;	/* Maximum number of blocks per
					 * sequence */
	double        **dAlignCount;	/* Count of number of possible
					 * alignments  */
	double        **dAlignCountSave;
	double          dTotalAlgnCnt;	/* total alignments */
	double         *dAlignWt;	/* weights on alignments for priors */
}               AlignStruct;

typedef AlignStruct *ALType;	/* BT 1/16/98 */


typedef struct {
	double       ***dBkgndProb;	/* bkgndProb[seq][n][pos] n =
					 * 0...alphalen-1 */
}               BkgndStruct;

typedef BkgndStruct *BkgType;


typedef struct {
	double          frequency;	/* sampling frequency of motif  */
	double          rev_freq;	/* Reverse Complement frequency */
	double          fwd_freq;	/* Forward frequency */
}               Frequency;


typedef struct {
	double          dProbability;
	int            *nNumMotifs;
	int           **nMotifLoc;
	short         **RevComp;
	long            nseed;
	int             nIterationNum;
	double        **dvMotifProb;
	Frequency     **frequency;	/* frequency structure */
	Ftype           F;
	int            *nMotifLen;	/* BT 5/23/97 */
	double         *dMap;	/* BT 3/4/98 */
	int             nMaxLen;
	double          dTotalMap;
	double          dTotalNonPalinMap;
	double          dTotalNonPhyloMap;
	double          dBkgndMap;
	double          dBetaMap;
	double         *dFragMap;
	double          dSeqMap;
	int             nSuboptSeed;
}               MaxResults;


typedef struct {
	int             nTemps;
	double         *dTemp;
	MaxResults     *results;
	double          dMinTemp;
	double          dMaxTemp;
	double          currTemp;
	int             nExchangeIterations;
	double          dTempStep;	/* BT 12/20/2000 */
	int             bExchange;	/* BT 2/20/2001 */
	clock_t         exchangePeriod;
	clock_t         lastExchangeTime;
}               AnealStruct;

typedef AnealStruct *AnlType;


typedef struct {
	int             code;
	char            id[WT_TABLE_ID_SIZE];
	double         *weight;
}               WeightStruct;


typedef struct {
	int             nTableSize;
	int             nSeqTypes;
	int            *seqCodes;
	char            seqIdStr[WT_TABLE_ID_SIZE * 2 + 1];
	WeightStruct   *weightTable;
	double         *weight;
}               SeqWeights;


typedef struct {
	IPtype          IP;
	Ctype           C;
	Ctype           First;
	ProbStruct     *P;
	files          *Datafiles;
	Ftype           F;
	Stype           Seq;
	RPType          RP;	/* BT 12/12/97 */
	ALType          AP;	/* BT 1/16/98 */
	PoSition      **InitPos;/* BT 8/5/98 */
	BkgType         BP;
	AnlType         AN;	/* BT 3/30/99 */
	SeqWeights     *WT;	/* BT 4/17/00 */
	PhyloType       Phylo;
}               Modelstruct;

typedef Modelstruct *Model;

typedef struct {
	long            begin;
	long            end;
	long            total;
}               SimTime;

typedef struct {
	double          abs_logl;
	double          pct_logl;
	double          abs_tot;
	double          pct_tot;
	double          last_logl;
	double          last_tot;
	double          new_tot;
	double          logl;
	double          likelihood;
	double          cost;
}               EMStats;

typedef struct MC {
	char           *desc;
	char           *cmp_desc;
	double          logl;
	double          cmp_logl;
	double          dof;
	double          lrt;
	double          pval;
	struct MC      *next;
}               MCstruct;


/*-------------------------------------------------------------------------*/
/*-------------------------------FUNCTION PROTOTYPES-----------------------*/
/*-------------------------------------------------------------------------*/

short           SEQ2CH(short);
void            print_usage_bernoulli(void);
void            print_usage_sankoff(void);
/*
 * void   p_error (char *msg); void   p_usage_error (char *msg);
 */
int             findMaxNumMotif(IPtype IP);
int             findMaxMotifLen(IPtype IP);
char            complement(char ch);
void            sRandom(Model B, long x);
long            Random(void);
double          drand(void);	/* BT 5/9/97 */
int             RandomInterval(int a, int b);	/* BT 5/9/97 */
void            rm_file(char *fname);
int             FileExists(char *fname);
void            print_EnterArgs(void);
void            BeginTime(SimTime * S);
void            EndTime(Model B, SimTime * S, FILE * fpt);
void            FreeData(Model B);
int             FileSize(char *filename);

#endif
