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
/* $Id: common.c,v 1.46 2009/08/12 13:38:24 Bill Exp $                    */
/* */
/* Author:   Eric C. Rouchka, 1996.6                                      */
/* Jun Zhu, 1996.8                                              */
/* Bill Thompson  1997.1                                        */
/* */
/* Description :  common subroutines                                      */
/**************************************************************************/


#include <limits.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include "common.h"
#ifdef _GUI_
#include "GUI_common.h"
#endif
#include "rgamma.h"

#define r_off   12
#define DIM(A) (sizeof(A)/sizeof((A)[0]))
/**----------------------------------------------------------------------**/
/**                        RANDOM NUMBER GENERATOR DECLARATIONS          **/
/**----------------------------------------------------------------------**/

static long     state[33] = {
	(long) 0xd53f1852, (long) 0xdfc78b83, (long) 0x4f256096, (long) 0xe643df7,
	(long) 0x82c359bf, (long) 0xc7794dfa, (long) 0xd5e9ffaa, (long) 0x2c8cb64a,
	(long) 0x2f07b334, (long) 0xad5a7eb5, (long) 0x96dc0cde, (long) 0x6fc24589,
	(long) 0xa5853646, (long) 0xe71576e2, (long) 0xdae30df, (long) 0xb09ce711,
	(long) 0x5e56ef87, (long) 0x4b4b0082, (long) 0x6f4f340e, (long) 0xc5bb17e8,
	(long) 0xd788d765, (long) 0x67498087, (long) 0x9d7aba26, (long) 0x261351d4,
	(long) 0x411ee7ea, (long) 0x393a263, (long) 0x2c5a5835, (long) 0xc115fcd8,
	(long) 0x25e9132c, (long) 0xd0c6e906, (long) 0xc2bc5b2d, (long) 0x6c065c98,
(long) 0x6e37bd55};

static long    *rJ = &state[r_off], *rK = &state[
						 DIM(state) - 1];

/* Global variables for FREEP and FREEPP loop counters - BT 1/21/97 */

	short           GL_GUI = 1;	/* Set to 1 if run from GUI */

	int             GL_I__;
	int             GL_J__;
	int             GL_Finished;

	int             glVerMajor = 3;	/* BT 3/28/97 */
	int             glVerMinor = 10;
	int             glRev = 1;
	char            glRevDate[12] = __DATE__;

	int             nComp[4] = {1, 0, 3, 2};

	int             GUI_argc;
	char          **GUI_argv;

/**************************** SEQ2CH ************************************/
/* For convenience in calculation, we enumerate amino acid to 1 to 20   */
/* corresponding to a to t.  Some special symbols we need to take care  */
/* are V,W,Y.                                                           */
/************************************************************************/

	short           SEQ2CH(short x)
{
	if (x == 1) {		/* symbol for V  */
		return 21;
	} else if (x == 9) {	/* symbol for W */
		return 22;
	} else if (x == 14) {	/* symbol for Y */
		return 24;
	} else			/* other symbols */
		return x;
}


void 
print_usage_bernoulli(void)
{
#ifdef _DEBUG_
	printf("\nGibbs %d.%02d.%03d  %s Debugging Version\n", glVerMajor, glVerMinor, glRev, glRevDate);
#elif defined _MPI_
	printf("\nGibbs %d.%02d.%03d  %s MPI Version\n", glVerMajor, glVerMinor, glRev, glRevDate);
#else
	printf("\nGibbs %d.%02d.%03d  %s\n", glVerMajor, glVerMinor, glRev, glRevDate);
#endif
	printf("USAGE (site sampler) :: Gibbs file lengths {flags}\n");
	printf("USAGE (motif sampler):: Gibbs file lengths expect {flags}\n");
	printf("USAGE (recursive sampler):: Gibbs file lengths expect -E max_sites {flags}\n");
	printf("USAGE (centroid sampler):: Gibbs file lengths expect -bayes -E max_sites {flags}\n\n");
	printf("lengths = <int>[,<int>] : width of motif to be found\n");
	printf("expect  = <int>[,<int>] : expect number of motif elements\n");
	printf("max_sites  = <int> : max. sites/seq\n");
	printf("\npossible flags:\n\n");

	printf("-A                         init sample from prior\n");
	printf("-B <bkgnd_filename>        Background Composition Model\n");
	printf("-C <cutoff_value>          cutoff for near optimal sampler\n");
	printf("-D <seqs[,aligned_seqs]>   Homologous sequences, seqs default = 2, aligned_segs def. = all seqs\n");
	printf("-E max_sites               Set max sites/seq, use recursive sampler\n");
	printf("-F                         Do not use fragmentation\n");
#ifndef _PUBLIC2_
	printf("-G                         Group Sampler\n");
	printf("-H <weight_filename>       Sequence Weight File\n");
#endif				/* _PUBLIC2_ */
	printf("-I <mnum, beg, end>*       direct repeat model between beg and end\n");
	printf("-J                         Fragment sites in center only\n");
	printf("-K <map>                   Alt. method of sampling sites/seq, min map to start (optional) (-E option only)\n");
	printf("-L                         Motif sample before recursive (-E only)\n");
	printf("-M <mnum, width>*          maximum widths for fragmentation\n");
	printf("-N <scan_filename>         output data for Scan\n");	/* BT 9/17/97 */
	printf("-O <prout_filename>        output informative priors\n");
	printf("-P <prior_filename>        file of informative priors\n");
	printf("-Q <sample_filename>       save sample counts in file\n");	/* BT 9/17/97 */
	printf("-R <mnum, beg, end>*       palindromic model between beg and end\n");
	printf("-S <num_Seeds>             number of seeds to try\n");
	printf("-U <spacing_filename>      Spacing Model\n");
	printf("-V <no. of seqs>           Verify Mode\n");
	printf("-W <pseudo_wt>             pseudosite weight (between 0 and 1)\n");
	printf("-X <min><,max<,step,<t>>>  Parallel Tempering (MPI version only)\n");
	printf("-Y                         Calculate default pseudocount weight\n");
	printf("-Z                         Don't write progress info\n");
	printf("-a <mnum, beg, end, pal>*  Concentrate between beg and end\n");
	printf("-align_centroid            Align the centroid sites\n");
	printf("-b                         Sample from background\n");
	printf("-bayes <burnin,<samples>>  Bayesian sampling\n");
	printf("-bayes_counts              Print counts from Bayesian sampling in output\n");
	printf("-c <mnum, beg, end>*       collapse alphabet between beg and end\n");
	printf("-d <mnum, min, max>*       Allow width to vary\n");	/* BT 6/27/97 */
	/*
	 * printf("-e                         Run Expectation/maximization
	 * algorithm\n");
	 *//* BT 09/22/04 */
	printf("-f <cutoff_factor>         cutoff factor for recursive sampler\n");
	printf("-frag                      alternate fragmentation sampling\n");
	printf("-freq                      Print Frequency Solution (default for optimal models)\n");
	printf("-g                         Sample along length of site\n");	/* BT 9/10/97 */
	printf("-h                         this message\n");
	printf("-hierarchical_model <iter> hierarchical model for sites/seq\n");
	printf("-i <iteration_num>         number of iterations to try\n");
	printf("-j <period>                Frag/Shift period\n");	/* BT 9/11/98 */
	printf("-k <iteration_num>         number of iterations to sample after plateau\n");
	printf("-l <control_filename>      Wilcoxon Signed Rank test\n");	/* BT 7/7/97 */
	printf("-m                         Do not maximize after near optimal sampling\n");
	printf("-n                         Use nucleic acid alphabet\n");
	printf("-nopt                      Don't print Near Optimal output\n");
	printf("-nopt_disp <value>         Min probability for Nearopt display\n");
	printf("-no_cred                   Don't calculate credibility limits for centroid.\n");
	printf("-o <out_filename>          file where results will be written\n");
	printf("-opt                       Print Near Optimal output (default for optimal models)\n");
	printf("-p <plateau_per>           number of periods a maximum value hasn't changed\n");
	printf("-pred_update               use preditive update to estimate motif model\n");
	printf("-q                         Sample width counts\n");	/* BT 3/25/98 */
	printf("-r                         turn off reverse complements with DNA\n");
	printf("-s <seedval>               random number generator seed \n");
	printf("-sample_model              sample motif model from dirichlet\n");
	printf("-t                         Display sites used in near optimal sampling\n");
	printf("-u                         Display output from suboptimal sampler\n");	/* BT 3/19/97 */
	printf("-v <overlap_value>         %% to allow overlap at ends of sequence\n");	/* BT 4/23/97 */
	printf("-w <pseduo_cnt_wt>         pseduocount weight\n");
	printf("-wilcox                    Wilcoxon sequences included in fasta file.\n");
	printf("-x                         Do not remove low complexity regions\n");
	printf("-y                         Don't print frequency solution\n");
	printf("\n");
}

void 
print_usage_sankoff(void)
{
	printf("Gibbs -PSankoff -INFILE=inputfile\n");
	printf(" [-MAX_BLOCK=max_blocks]\n");
	printf(" [-MATRIX_NUM=matrix_num]\n");
	printf(" [-OUTFILE=outputfile]\n");
	printf(" [-BEST_ALIGN]\n");
	printf(" [-OUT_SEQ]\n");
	printf(" [-BACK_SAMPLE]\n");
	printf(" [-BACK_S_BLOCKS=back_sample_blocks]\n");
	printf(" [-BACK_S_NUMBER=back_sample_number]\n");
	printf(" [-BACK_S_OUTSEQ] \n");
	printf(" [-BACK_CUTOFF]\n");
	printf(" [-BACK_C_VALUE=back_cutoff_value]\n");
}


/*
 * void print_EnterArgs(void) { printf("\n\n\n\tARGUMENTS SHOULD BE ENTERED
 * IN THE FOLLOWING USAGE :: \n\n"); printf("(site sampler)  :: file lengths
 * {flags}\n"); printf("(motif sampler) :: file lengths expect {flags}\n");
 * print_flags(); printf("\n\nEnter in the arguments followed by <return> :
 * \n"); }
 */


/************************** findMaxMotifLen **********************/
/* */
/* Find maximum motif length                                   */
/*****************************************************************/

int 
findMaxMotifLen(IPtype IP)
{
	int             i;
	int             max = 0;
	for (i = 0; i < IP->nNumMotifTypes; i++) {
		if (IP->nMaxMotifLen[i] > max)
			max = IP->nMaxMotifLen[i];
	}
	return max;
}

/************************** findMaxNumMotif  *********************/
/* */
/* Find maximum number of sites in a single type motif         */
/*****************************************************************/

int 
findMaxNumMotif(IPtype IP)
{
	int             i;
	int             max = 0;
	int             value;

	for (i = 0; i < IP->nNumMotifTypes; i++) {
		value = IP->nNumMotifs[i][FORWARD];
		if (IP->RevComplement)
			value += IP->nNumMotifs[i][REVERSE];
		if (value > max)
			max = value;
	}
	return max;		/* BT 1/2/2001 */
}

/******************* complement *****************************/
/* */
/* find the complementary nucleotide                    */
/************************************************************/

char 
complement(char ch)
{
	static char    *dnaAlpha = "abcdnxX";
	static char    *compDna = "badcnee";
	char           *p;

	if ((p = strchr(dnaAlpha, ch)) != NULL) {
		return compDna[p - dnaAlpha];
	} else {
		printf("CH = %c", ch);
		p_internal_error("Invalid Nucleotide");
	}

	return '\0';
}


char 
xcomplement(char ch)
{
	char            retval;

	switch (ch) {
	case 'a':
		retval = 'b';
		break;
	case 'b':
		retval = 'a';
		break;
	case 'c':
		retval = 'd';
		break;
	case 'd':
		retval = 'c';
		break;
	case 'n':
		/* retval = 'e';  *//* represent X here */
		retval = 'n';	/* represent X here *//* BT 5/12/97 */
		break;
	case 'x':
		retval = 'e';	/* represent X here */
		break;
	case 'X':
		retval = 'e';	/* represent X here */
		break;
	default:
		retval = '\0';
		printf("CH = %c", ch);
		p_internal_error("Invalid Nucleotide");
		break;
	}
	return retval;
}


/*
        Additive random number generator

        Modelled after "Algorithm A" in
        Knuth, D. E. (1981). The art of computer programming, volume 2, page 27
.

        7/26/90 WRG
*/

void 
sRandom(Model B, long x)
{
	register long   i;

	state[0] = x;
	/* linear congruential initializer */
	for (i = 1; i < DIM(state); ++i) {
		state[i] = 1103515245 * state[i - 1] + 12345;
	}

	rJ = &state[r_off];
	rK = &state[DIM(state) - 1];

	for (i = 0; i < 10 * DIM(state); ++i)
		Random();

	zigset(x);
}


/*
        Random --  return value in the range 0 <= x < 2**31 - 1
*/
long 
Random(void)
{
	register long   r;

	r = *rK;
	r += *rJ--;
	*rK-- = r;

	if (rK < state)
		rK = &state[DIM(state) - 1];
	else if (rJ < state)
		rJ = &state[DIM(state) - 1];
	return (r >> 1) & 0x7fffffff;	/* discard the least-random bit */
}

/* return a random double between 0 (inclusive) and 1.0 (exclusive) */
double 
drand(void)
{				/* BT 5/9/97 */
	return ((double) Random() / ((double) 0x7fffffff + 1.0));	/* BT 11/10/97 */
}


/* Return a random integer between a and b inclusive (a, b > 0) */
int 
RandomInterval(int a, int b)
{				/* BT 5/9/97 */
	int             nRes;

	nRes = floor(drand() * (b + 1 - a) + a);
	return nRes;
}


void 
rm_file(char *fname)
{
	remove(fname);
}


int 
FileExists(char *fname)
{
	struct stat     st;

	if (stat(fname, &st) == 0)
		return TRUE;
	else
		return FALSE;
}


/*******************  BeginTime  *****************/

void 
BeginTime(SimTime * S)
{
	struct timeval  tp;

	gettimeofday(&tp, NULL);
	S->begin = tp.tv_sec;
}


/*******************  EndTime   ******************/
void 
EndTime(Model B, SimTime * S, FILE * fpt)
{
	struct timeval  tp;

	gettimeofday(&tp, NULL);
	S->end = tp.tv_sec;
	S->total = S->end - S->begin;
	if (!B->IP->is_defined[cl_Z])
		fprintf(fpt, "Total Time %ld sec (%f min)\n", S->total, (double) S->total / 60.0);
}


int 
FileSize(char *filename)
{
	struct stat     st;

	if (!stat(filename, &st))
		return (st.st_size);
	else
		return 0;
}
