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
/* $Id: parameters.c,v 1.29 2009/08/12 13:38:24 Bill Exp $         */
/* */
/* Author:   Eric C. Rouchka, 1996,6.                                     */
/* Jun Zhu, 1996.8                                              */
/* Bill Thompson, 1/24/96                                       */
/* Description :  This file contains the functions that are necessary for */
/* stripping the command line arguments and setting the    */
/* appropriate variables                                   */
/**************************************************************************/

#include "parameters.h"

/*------------------------------------------------------------------------*/
/*------------------------local function prototypes-----------------------*/
/*------------------------------------------------------------------------*/
void            strip_palandromes(int argc, char **argv, IPtype IP, int *pos);
void            strip_repeats(int argc, char **argv, IPtype IP, int *pos);
void            strip_collapsed(int argc, char **argv, IPtype IP, int *pos);
void            strip_cutoff(int argc, char **argv, IPtype IP, int *i);
void            strip_plateau(int argc, char **argv, IPtype IP, int *i);
void            strip_adjust(int argc, char **argv, IPtype IP, int *i);
void            strip_seedval(int argc, char **argv, IPtype IP, int *i);
void            strip_seednum(int argc, char **argv, IPtype IP, int *i);
void            strip_sitewt(int argc, char **argv, IPtype IP, int *i);
void            strip_countwt(int argc, char **argv, IPtype IP, int *i);
void            strip_iterations(int argc, char **argv, IPtype IP, int *i);
short           strip_outfile(int argc, char **argv, IPtype IP, int *i, short bGui);	/* BT 4/25/97 */
short           strip_priorfile(int argc, char **argv, IPtype IP, int *i, short bGui);	/* BT 4/25/97 */
short           strip_scanfile(int argc, char **argv, IPtype IP, int *i, short bGui);	/* BT 9/17/97 */
short           strip_proutfile(int argc, char **argv, IPtype IP, int *i, short bGui);	/* BT 9/26/97 */
void            strip_concen(int argc, char **argv, IPtype IP, int *pos);
void            strip_overlap(int argc, char **argv, IPtype IP, int *i);	/* BT 4/23/97 */
void            strip_width(int argc, char **argv, IPtype IP, int *pos);	/* BT 7/18/97 */
short           strip_wilcoxon(int argc, char **argv, IPtype IP, int *i, short bGui);	/* BT 1/28/98 */
void            strip_maxblocks(int argc, char **argv, IPtype IP, int *i);	/* 2/11/98 */
void            strip_rcutoff(int argc, char **argv, IPtype IP, int *i);	/* 2/11/98 */
void            strip_samplecount(int argc, char **argv, IPtype IP, int *i);	/* BT 2/20/98 */
short           strip_samplefile(int argc, char **argv, IPtype IP, int *i, short bGui);
short           strip_bkgndfile(int argc, char **argv, IPtype IP, int *i, short bGui);
short           strip_spacingfile(int argc, char **argv, IPtype IP, int *i, short bGui);
short           strip_weightfile(int argc, char **argv, IPtype IP, int *i, short bGui);
void            strip_sample_iterations(int argc, char **argv, IPtype IP, int *i);
void            strip_frag_width(int argc, char **argv, IPtype IP, int *i);
void            strip_exchange(int argc, char **argv, Model B, int *i);
void            strip_ksamplemax(int argc, char **argv, Model B, int *i);
void            strip_homolog_seqs(int argc, char **argv, Model B, int *i);
void            strip_verify(int argc, char **argv, IPtype IP, int *i);
void            strip_nearopt_cutoff(int argc, char **argv, IPtype IP, int *i);
void            strip_hierarchical(int argc, char **argv, Model B, int *i);
void            strip_bayes(int argc, char **argv, Model B, int *i);

/*------------------------------------------------------------------------*/



/* =================================================================== */
/* FUNCTION NAME : strip_palandromes                                 */
/* */
/* DESCRIPTION : this function strips the palandromic information    */
/* =================================================================== */

void 
strip_palandromes(int argc, char **argv, IPtype IP, int *pos)
{
	int             temp;
	int            *nval;
	int             t, beg, end, max_index;
	int             i, j;

	NEW(nval, 300, int);	/* BT 12/03/99 */
	if (argc < (*pos) + 2)
		p_error("Missing palindromic range arguments");
	temp = ParseIntegers(argv[(*pos) + 1], nval, "Getting palindromic ranges");
	if ((temp % 3) != 0)
		p_error("Mismatched Palindrome begin and end ranges");
	for (i = 0; i < temp / 3; i++) {
		t = nval[i * 3 + 1] - 1;
		beg = nval[i * 3 + 2];
		end = nval[i * 3 + 3];
		max_index = IP->nMotifLen[t] / 2 + IP->nMotifLen[t] % 2;

		if ((beg < 1) || (end < 1))
			p_error("Palandrome subscript underflow");
		if ((beg > max_index) || (end > max_index))
			p_error("Palandrome subscript overflow");
		if (beg > end)
			p_error("Palandrome begin pos is greater than end pos");
		beg--;
		end--;
		for (j = beg; j <= end; j++) {
			if (IP->AltModel->Repeat[t][j] || IP->AltModel->Repeat[t][IP->nMotifLen[t] - j - 1])
				p_error("Positions are defined as both palindromes and repeats.");
			IP->AltModel->Palandromic[t][j] = TRUE;
			IP->AltModel->Palandromic[t][IP->nMotifLen[t] - j - 1] = TRUE;
		}
	}
	(*pos)++;
	IP->is_defined[cl_R] = TRUE;
	free(nval);
}

/* =================================================================== */
/* FUNCTION NAME : strip_repeats                                     */
/* */
/* DESCRIPTION : this function strips the repeat information         */
/* =================================================================== */

void 
strip_repeats(int argc, char **argv, IPtype IP, int *pos)
{
	int             temp;
	int            *nval;
	int             t, beg, end, max_index;
	int             i, j;

	NEW(nval, 300, int);	/* BT 12/03/99 */
	if (argc < (*pos) + 2)
		p_error("Missing repeat range arguments");
	temp = ParseIntegers(argv[(*pos) + 1], nval, "Getting repeat ranges");
	if ((temp % 3) != 0)
		p_error("Mismatched repeat begin and end ranges");
	for (i = 0; i < temp / 3; i++) {
		t = nval[i * 3 + 1] - 1;
		beg = nval[i * 3 + 2];
		end = nval[i * 3 + 3];
		max_index = IP->nMotifLen[t] / 2 + IP->nMotifLen[t] % 2;

		if ((beg < 1) || (end < 1))
			p_error("Repeat subscript underflow");
		if ((beg > max_index) || (end > max_index))
			p_error("Repeat subscript overflow");
		if (beg > end)
			p_error("Repeat begin pos is greater than end pos");
		beg--;
		end--;
		for (j = beg; j <= end; j++) {
			if (IP->AltModel->Palandromic[t][j] || IP->AltModel->Palandromic[t][IP->nMotifLen[t] - j - 1])
				p_error("Positions are defined as both palindromes and repeats.");
			IP->AltModel->Repeat[t][j] = TRUE;
			IP->AltModel->Repeat[t][IP->nMotifLen[t] - j - 1] = TRUE;
		}
	}
	(*pos)++;
	IP->is_defined[cl_I] = TRUE;
	IP->is_defined[cl_R] = TRUE;
	free(nval);
}


void 
strip_concen(int argc, char **argv, IPtype IP, int *pos)
/* =================================================================== */
/* FUNCTION NAME : strip_concen                                      */
/* */
/* DESCRIPTION : this function strips the concentrated information   */
/* =================================================================== */

{
	int             temp;
	int            *nval;
	int             t, beg, end, max_index;
	int             palindrome;
	int             i, j;

	/** NEED TO CHANGE ALL OF THIS **/
	NEW(nval, 300, int);
	if (argc < (*pos) + 2)
		p_error("Missing concentrated region arguments");
	temp = ParseIntegers(argv[(*pos) + 1], nval, "Getting concentrated region");
	if ((temp % 4) != 0)
		p_error("Mismatched concentrated region arguments");
	for (i = 0; i < temp / 4; i++) {
		t = nval[i * 4 + 1] - 1;
		beg = nval[i * 4 + 2];
		end = nval[i * 4 + 3];
		IP->AltModel->Concen->NumConcen[t]++;
		palindrome = nval[i * 4 + 4];
		if (palindrome)
			max_index = IP->nMotifLen[t] / 2 + IP->nMotifLen[t] % 2;
		else
			max_index = IP->nMotifLen[t];

		if ((beg < 1) || (end < 1))
			p_error("Concentrated Region subscript underflow");
		if ((beg > max_index) || (end > max_index))
			p_error("Concentrated Region subscript overflow");
		if (beg > end)
			p_error("Concentrated Region begin pos is greater than end pos");
		beg--;
		end--;
		for (j = beg; j <= end; j++) {
			IP->AltModel->Concen->Concentrated[t][j] =
				IP->AltModel->Concen->NumConcen[t];
			if (palindrome)
				IP->AltModel->Concen->Concentrated[t][IP->nMotifLen[t] - j - 1] =
					IP->AltModel->Concen->NumConcen[t];
		}
	}
	(*pos)++;
	IP->is_defined[cl_a] = TRUE;
	free(nval);
}




void 
strip_collapsed(int argc, char **argv, IPtype IP, int *pos)
{
	/*
	 * ===================================================================
	 * ==
	 */
	/* FUNCTION NAME : strip_collapsed                                     */
	/* */
	/*
	 * DESCRIPTION : This function strips out the parameters that fix
	 * where
	 */
	/* the collapsed alphabet will occur within a motif model */
	/*
	 * ===================================================================
	 * ==
	 */

	int            *nval;
	int             temp, j, k;
	int             t, beg, end, max_index;
	NEW(nval, 300, int);

	if (argc < (*pos) + 2)
		p_error("Missing collapsed range arguments");
	temp = ParseIntegers(argv[(*pos) + 1], nval, "Getting collapsed range");
	if ((temp % 3) != 0)
		p_error("Mismatched collapsed region begin and end ranges");
	for (j = 0; j < temp / 3; j++) {
		t = nval[j * 3 + 1] - 1;
		beg = nval[j * 3 + 2];
		end = nval[j * 3 + 3];
		max_index = IP->nMotifLen[t];
		if ((beg < 1) || (end < 1))
			p_error("Collapsed subscript underflow");
		if ((beg > max_index) || (end > max_index))
			p_error("Collapsed subscript overflow");
		if (beg > end)
			p_error("Collapsed begin pos is greater than end pos");
		beg--;
		end--;
		for (k = beg; k <= end; k++)
			IP->AltModel->Collapsed[t][k] = TRUE;
	}
	(*pos)++;
	IP->is_defined[cl_c] = TRUE;
	free(nval);
}

void 
strip_cutoff(int argc, char **argv, IPtype IP, int *i)
/* =========================================================== */
/* FUNCTION NAME : strip_cutoff                              */
/* */
/* DESCRIPTION : strips off the value for the cutoff used    */
/* in near optimal sampling                    */
/* =========================================================== */

{
	double          dtemp;

	if (argc < (*i) + 2)
		p_usage_error("-C: MISSING ARGUMENTS");
	if (!ParseReals(argv[(*i) + 1], &dtemp))
		p_error("Error in forming value after -C");
	IP->dCutoff = dtemp;
	if ((IP->dCutoff >= 1.0) || (IP->dCutoff < 0))
		p_error("Near Opt cutoff must be between 0 and 1");
	IP->is_defined[cl_C] = TRUE;
	(*i)++;
}


void 
strip_plateau(int argc, char **argv, IPtype IP, int *i)
/* ============================================================ */
/* FUNCTION NAME : strip_plateau                              */
/* */
/* DESCRIPTION : strips off the value of the maximum number of */
/* plateau periods are calculated between       */
/* local maximas                                */
/* ============================================================ */

{
	int             temp;

	if (argc < (*i) + 2)
		p_usage_error("-p: MISSING ARGUMENTS");
	temp = atoi(argv[(*i) + 1]);
	if (temp < 5)
		p_error("Plateau periods must be >= 5");
	IP->nPlateauPeriods = temp;
	IP->is_defined[cl_p] = TRUE;
	(*i)++;
}


void 
strip_adjust(int argc, char **argv, IPtype IP, int *i)
/* ============================================================ */
/* FUNCTION NAME : strip_adjust                               */
/* */
/* DESCRIPTION : strips off the value of the maximum number of */
/* iterations between fragmentation or shif     */
/* ============================================================ */

{
	int             temp;

	if (argc < (*i) + 2)
		p_usage_error("-j: MISSING ARGUMENTS");
	temp = atoi(argv[(*i) + 1]);
	if (temp < 0)
		p_error("Adjust period must be greater than 0");
	IP->nAdjustPeriod = temp;
	IP->is_defined[cl_j] = TRUE;
	(*i)++;
}

void 
strip_seedval(int argc, char **argv, IPtype IP, int *i)
/* =============================================================== */
/* FUNCTION NAME : strip_seedval                                 */
/* */
/* DESCRIPTION : strips of the initial seed value to be used     */
/* =============================================================== */

{
	if (argc < (*i) + 2)
		p_usage_error("-s: MISSING ARGUMENTS");
	IP->lSeedVal = atol(argv[(*i) + 1]);
	IP->is_defined[cl_s] = TRUE;
	(*i)++;
}

void 
strip_seednum(int argc, char **argv, IPtype IP, int *i)
/* =============================================================== */
/* FUNCTION NAME : strip_seednum                                 */
/* */
/* DESCRIPTION : strips of the number of random seeds to be used */
/* =============================================================== */

{
	int             temp;

	if (argc < (*i) + 2)
		p_usage_error("-S: Missing Arguments");
	temp = atoi(argv[(*i) + 1]);
	if (temp < 1)
		p_error("Must try at least one seed");
	IP->nSeeds = temp;
	IP->is_defined[cl_S] = TRUE;
	(*i)++;
}

void 
strip_sitewt(int argc, char **argv, IPtype IP, int *i)
/* ============================================================= */
/* FUNCTION NAME : strip_sitewt                                */
/* */
/* DESCRIPTION : strips the pseudo site weight to be used      */
/* ============================================================= */

{
	double          dtemp;

	if (argc < (*i) + 2)
		p_usage_error("-W: MISSING ARGUMENTS");
	if (!ParseReals(argv[(*i) + 1], &dtemp))
		p_error("Error in forming value after -W");
	IP->dPseudoSiteWt = dtemp;
	if ((IP->dPseudoSiteWt >= 1.0) || (IP->dPseudoSiteWt < 0))
		p_error("Pseudosite weight must be between 0 and 1");
	IP->is_defined[cl_W] = TRUE;
	(*i)++;
}

void 
strip_countwt(int argc, char **argv, IPtype IP, int *i)
/* ============================================================= */
/* FUNCTION NAME : strip_countwt                               */
/* */
/* DESCRIPTION : strips the pseudo count weight to be used     */
/* ============================================================= */
{
	double          dtemp;

	if (argc < (*i) + 2)
		p_usage_error("-w: MISSING ARGUMENTS");
	if (!ParseReals(argv[(*i) + 1], &dtemp))
		p_error("Error in forming value after -w");	/* BT 1/31/97 */
	IP->dPseudoCntWt = dtemp;
	if ((IP->dPseudoCntWt > 1.0) || (IP->dPseudoCntWt < 0))
		p_error("Pseudocount weight must be defined between 0 and 1");
	IP->is_defined[cl_w] = TRUE;
	(*i)++;
}

void 
strip_iterations(int argc, char **argv, IPtype IP, int *i)
/* ============================================================= */
/* FUNCTION NAME : strip_iterations                            */
/* */
/* DESCRIPTION : strips the maximum number of iterations used  */
/* in finding the near optimal solution          */
/* ============================================================= */
{
	int             temp;

	if (argc < (*i) + 2)
		p_usage_error("-i: MISSING ARGUMENTS");
	temp = atoi(argv[(*i) + 1]);
	if (temp < 1)
		p_error("Iteration count must be at least 1");
	IP->nMaxIterations = temp;
	IP->is_defined[cl_i] = TRUE;
	(*i)++;
}

short 
strip_outfile(int argc, char **argv, IPtype IP, int *i, short bGui)
/* ============================================================== */
/* FUNCTION NAME : strip_outfile                                */
/* */
/* DESCRIPTION : strips the name of the output file where the   */
/* results will be printed                        */
/****************************************************************/
{
	char            tmpstr[256];

	if (argc < (*i) + 2)
		p_usage_error("-o: Missing Argument");
	if (IP->nRank == 0)
		IP->Datafiles->out_fpt = fopen(argv[(*i) + 1], "w");
	else
		IP->Datafiles->out_fpt = fopen(argv[(*i) + 1], "a");

	if (IP->Datafiles->out_fpt == NULL) {	/* BT 4/25/97 */
		sprintf(tmpstr, "Cannot open output file %s", argv[(*i) + 1]);
		fprintf(stderr, "In strip_outfile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	IP->is_defined[cl_o] = TRUE;
	NEW(IP->Datafiles->output_filename, strlen(argv[(*i) + 1]) + 1, char);
	strcpy(IP->Datafiles->output_filename, argv[(*i) + 1]);
	(*i)++;

	return TRUE;
}

short 
strip_priorfile(int argc, char **argv, IPtype IP, int *i, short bGui)
/* =============================================================== */
/* FUNCTION NAME : strip_priorfile                               */
/* */
/* DESCRIPTION : strips the name of the file where informative   */
/* priors will be read from                        */
/* =============================================================== */

{
	char            tmpstr[255];

	if (argc < (*i) + 2)
		p_usage_error("-P: Minning argument");
	IP->Datafiles->prior_fpt = fopen(argv[(*i) + 1], "r");
	if (IP->Datafiles->prior_fpt == NULL) {
		sprintf(tmpstr, "Cannot open prior file %s", argv[(*i) + 1]);
		fprintf(stderr, "In strip_priorfile args\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);	/* BT 3/21/97 */
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	IP->is_defined[cl_P] = TRUE;
	NEW(IP->Datafiles->prior_filename, strlen(argv[(*i) + 1]) + 1, char);
	strcpy(IP->Datafiles->prior_filename, argv[(*i) + 1]);
	(*i)++;
	return TRUE;
}


short 
strip_scanfile(int argc, char **argv, IPtype IP, int *i, short bGui)
{				/* BT 9/17/97 */
	/* ============================================================== */
	/* FUNCTION NAME : strip_scanfile                               */
	/* */
	/* DESCRIPTION : strips the name of the file where scan the     */
	/* results will be printed                        */
	/****************************************************************/
	char            tmpstr[256];

	if (argc < (*i) + 2)
		p_usage_error("-N: MISSING ARGUMENTS");
	IP->Datafiles->sn_fpt = fopen(argv[(*i) + 1], "w");

	if (IP->Datafiles->sn_fpt == NULL) {	/* BT 4/25/97 */
		sprintf(tmpstr, "Cannot open scan output file %s", argv[(*i) + 1]);
		fprintf(stderr, "In strip_scanfile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);	/* BT 3/21/97 */
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	IP->is_defined[cl_N] = TRUE;
	NEW(IP->Datafiles->ScanFilename, strlen(argv[(*i) + 1]) + 1, char);	/* BT 1/27/97 */
	strcpy(IP->Datafiles->ScanFilename, argv[(*i) + 1]);
	(*i)++;

	return TRUE;
}


short 
strip_proutfile(int argc, char **argv, IPtype IP, int *i, short bGui)
{				/* BT 9/26/97 */
	/* ============================================================== */
	/* FUNCTION NAME : strip_proutfile                              */
	/* */
	/* DESCRIPTION : strips the name of the file where prior info   */
	/* results will be printed                        */
	/****************************************************************/
	char            tmpstr[256];

	if (argc < (*i) + 2)
		p_usage_error("-O: MISSING ARGUMENTS");
	IP->Datafiles->prout_fpt = fopen(argv[(*i) + 1], "w");

	if (IP->Datafiles->prout_fpt == NULL) {	/* BT 4/25/97 */
		sprintf(tmpstr, "Cannot open prior output file %s", argv[(*i) + 1]);
		fprintf(stderr, "In strip_proutfile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);	/* BT 3/21/97 */
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	IP->is_defined[cl_O] = TRUE;
	NEW(IP->Datafiles->PriorOutFilename, strlen(argv[(*i) + 1]) + 1, char);	/* BT 1/27/97 */
	strcpy(IP->Datafiles->PriorOutFilename, argv[(*i) + 1]);
	(*i)++;

	return TRUE;
}


void 
strip_overlap(int argc, char **argv, IPtype IP, int *i)
{				/* BT 4/23/97 */
	/* =========================================================== */
	/* FUNCTION NAME : strip_overlap                             */
	/* */
	/* DESCRIPTION : strips off the value for the overlap %      */
	/* =========================================================== */
	double          dtemp;

	if (argc < (*i) + 2)
		p_usage_error("-v: MISSING ARGUMENTS");
	if (!ParseReals(argv[(*i) + 1], &dtemp))
		p_error("Error in forming value after -v");
	IP->glOverlapParam = dtemp;
	if ((IP->glOverlapParam >= 1.0) || (IP->glOverlapParam < 0))
		p_error("Overlap value must be between 0 and 1");
	IP->is_defined[cl_v] = TRUE;
	(*i)++;
}


void 
strip_width(int argc, char **argv, IPtype IP, int *pos)
{
	/*
	 * ===================================================================
	 * ==
	 */
	/* FUNCTION NAME : strip_width                                         */
	/* */
	/*
	 * DESCRIPTION : This function strips out the parameters that fix
	 * where
	 */
	/* the min and max width values.                         */
	/*
	 * ===================================================================
	 * ==
	 */

	int            *nval;
	int             temp, i;
	int             minWidth, maxWidth;
	int             t;
	char            cval[128];

	NEW(nval, 300, int);

	if (argc < (*pos) + 2)
		p_error("Missing width arguments");
	temp = ParseIntegers(argv[(*pos) + 1], nval, "Getting width");
	if ((temp % 3) != 0)
		p_error("Mismatched width ranges");
	for (i = 0; i < temp / 3; i++) {
		t = nval[i * 3 + 1] - 1;
		minWidth = nval[i * 3 + 2];
		if (minWidth < MIN_MOTIF_LEN) {
			sprintf(cval, "Min width must be greater than %d.", MIN_MOTIF_LEN);
			p_error(cval);
		}
		maxWidth = nval[i * 3 + 3];
		if (maxWidth < minWidth)
			p_error("Min width greater than max width");
		if (IP->nMotifLen[t] < minWidth) {
			sprintf(cval, "Min width must be at least %d.", IP->nMotifLen[t]);
			p_error(cval);
		}
		if (IP->nMotifLen[t] > maxWidth) {
			sprintf(cval, "Max width must be greater than %d.", IP->nMotifLen[t]);
			p_error(cval);
		}
		IP->nMinMotifLen[t] = minWidth;
		IP->nMaxMotifLen[t] = maxWidth;
	}

	IP->is_defined[cl_d] = TRUE;
	IP->is_defined[cl_F] = TRUE;
	(*pos)++;
	free(nval);
}


/* BT 1/28/98 */
short 
strip_wilcoxon(int argc, char **argv, IPtype IP, int *i, short bGui)
/* =============================================================== */
/* FUNCTION NAME : strip_wilcoxon                                */
/* */
/* DESCRIPTION : strips the name of the Wilcoxon test control    */
/* file. May be blank                               */
/* =============================================================== */

{
	char            tmpstr[255];

	if ((argc > (*i) + 1) && (argv[(*i) + 1][0] != '-')) {
		IP->Datafiles->control_fpt = fopen(argv[(*i) + 1], "r");
		if (IP->Datafiles->control_fpt == NULL) {
			sprintf(tmpstr, "Cannot open Wilcoxon Control file %s", argv[(*i) + 1]);
			fprintf(stderr, "In strip_priorfile args\n");
#ifdef _GUI_
			if (bGui) {
				wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
				fprintf(stderr, tmpstr);
			} else
				p_error(tmpstr);
#else
			p_error(tmpstr);
#endif
			(*i)++;
			return FALSE;
		} else {
			NEW(IP->Datafiles->ControlFileName, strlen(argv[(*i) + 1]), char);
			strcpy(IP->Datafiles->ControlFileName, argv[(*i) + 1]);
			(*i)++;
		}
	}
	IP->is_defined[cl_l] = TRUE;
	IP->is_defined[cl_wilcox] = TRUE;
	return TRUE;
}



void 
strip_maxblocks(int argc, char **argv, IPtype IP, int *i)
/* =============================================================== */
/* FUNCTION NAME : strip_maxblocks                               */
/* */
/* DESCRIPTION : strips off the number of motifs per sequence    */
/* for the recursive sampler                       */
/* =============================================================== */

{
	int             temp;
	int             nval[3];

	if (argc < (*i) + 2)
		p_usage_error("-E: MISSING ARGUMENTS");

	temp = ParseIntegers(argv[(*i) + 1], nval, "Getting max sites");
	if (temp < 1)
		p_error("strip_maxblocks: -E must have at least one motif per sequence");

	if (temp > 1) {
		IP->nMinBlocks = nval[1];
		IP->nMaxBlocks = nval[2];
		if (IP->nMinBlocks > IP->nMaxBlocks) {
			p_error("strip_maxblocks: -E min sites/seq must be <= max.");
		}
	} else {
		IP->nMinBlocks = 0;
		IP->nMaxBlocks = nval[1];
	}

	(*i)++;
	IP->is_defined[cl_E] = TRUE;
	IP->site_samp = FALSE;
}


void 
strip_rcutoff(int argc, char **argv, IPtype IP, int *i)
/* =========================================================== */
/* FUNCTION NAME : strip_rcutoff                              */
/* */
/* DESCRIPTION : strips off the value for the cutoff used    */
/* in recursive sampling                       */
/* =========================================================== */

{
	double          dtemp;

	if (argc < (*i) + 2)
		p_usage_error("-f: MISSING ARGUMENTS");
	if (!ParseReals(argv[(*i) + 1], &dtemp))
		p_error("Error in forming value after -f");
	IP->dRCutoff = dtemp;
	if ((IP->dRCutoff >= 1.0) || (IP->dRCutoff < 0))
		p_error("Recursion cutoff must be between 0 and 1");
	IP->is_defined[cl_f] = TRUE;
	(*i)++;
}


void 
strip_samplecount(int argc, char **argv, IPtype IP, int *i)
/* =============================================================== */
/* FUNCTION NAME : strip_samplecount                             */
/* */
/* DESCRIPTION : strips off the value of samples for post. width */
/* distribution                                    */
/* =============================================================== */

{
	if (argc < (*i) + 2)
		p_usage_error("-q: MISSING ARGUMENTS");
	IP->nSampleCnt = atol(argv[(*i) + 1]);
	IP->is_defined[cl_q] = TRUE;
	(*i)++;
}


short 
strip_samplefile(int argc, char **argv, IPtype IP, int *i, short bGui)
/* BT 3/25/98 */
/* ============================================================== */
/* FUNCTION NAME : strip_samplefile                             */
/* */
/* DESCRIPTION : strips the name of the file where sample inf   */
/* for alignment will be saved                    */
/****************************************************************/
{
	char            tmpstr[FILENAME_MAX + 256];
	char            tName[FILENAME_MAX + 256];

	if (argc < (*i) + 2)
		p_usage_error("-Q: MISSING ARGUMENTS");

	sprintf(tName, "%s.%d", argv[(*i) + 1], IP->nRank);
	IP->Datafiles->sample_fpt = fopen(tName, "w");

	if (IP->Datafiles->sample_fpt == NULL) {
		sprintf(tmpstr, "Cannot open sample output file %s", argv[(*i) + 1]);
		fprintf(stderr, "In strip_samplefile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	NEW(IP->Datafiles->SampleFileName, strlen(tName) + 1, char);
	strcpy(IP->Datafiles->SampleFileName, tName);

	IP->is_defined[cl_Q] = TRUE;

	sprintf(tName, "%s.occur.%d", argv[(*i) + 1], IP->nRank);

	IP->Datafiles->occur_fpt = fopen(tName, "w");

	if (IP->Datafiles->occur_fpt == NULL) {
		sprintf(tmpstr, "Cannot open occurence output file %s", tName);
		fprintf(stderr, "In strip_samplefile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	NEW(IP->Datafiles->OccurFileName, strlen(tName) + 1, char);
	strcpy(IP->Datafiles->OccurFileName, tName);

	sprintf(tName, "%s.occur.near.%d", argv[(*i) + 1],
		IP->nRank);

	IP->Datafiles->near_occur_fpt = fopen(tName, "w");

	if (IP->Datafiles->near_occur_fpt == NULL) {
		sprintf(tmpstr, "Cannot open near optimal occurence output file %s",
			tName);
		fprintf(stderr, "In strip_samplefile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif
		return FALSE;
	}
	NEW(IP->Datafiles->NearOccurFileName, strlen(tName) + 1, char);
	strcpy(IP->Datafiles->NearOccurFileName, tName);

	sprintf(tName, "%s.sites.%d", argv[(*i) + 1], IP->nRank);

	IP->Datafiles->sites_fpt = fopen(tName, "w");

	if (IP->Datafiles->sites_fpt == NULL) {
		sprintf(tmpstr, "Cannot open sites output file %s", tName);
		fprintf(stderr, "In strip_samplefile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif
		return FALSE;
	}
	fprintf(IP->Datafiles->sites_fpt, "seed iter motif seq pos seqLen width temp\n");

	NEW(IP->Datafiles->SitesFileName, strlen(tName) + 1, char);
	strcpy(IP->Datafiles->SitesFileName, tName);

	sprintf(tName, "%s.sitespr.%d", argv[(*i) + 1], IP->nRank);

	IP->Datafiles->sitespr_fpt = fopen(tName, "w");

	if (IP->Datafiles->sitespr_fpt == NULL) {
		sprintf(tmpstr, "Cannot open sites prior output file %s", tName);
		fprintf(stderr, "In strip_samplefile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif
		return FALSE;
	}
	fprintf(IP->Datafiles->sitespr_fpt,
		"seq seed iter probs (0...) selected\n");

	NEW(IP->Datafiles->SitesPrFileName, strlen(tName) + 1, char);
	strcpy(IP->Datafiles->SitesPrFileName, tName);

	sprintf(tName, "%s.prob.%d", argv[(*i) + 1], IP->nRank);

	IP->Datafiles->prob_fpt = fopen(tName, "w");

	if (IP->Datafiles->prob_fpt == NULL) {
		sprintf(tmpstr, "Cannot open probability output file %s", tName);
		fprintf(stderr, "In strip_samplefile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	fprintf(IP->Datafiles->prob_fpt,
		"seed iter motif pos pA pT pC pG temp\n");	/* BT 7/26/05 */

	NEW(IP->Datafiles->ProbFileName, strlen(tName) + 1, char);
	strcpy(IP->Datafiles->ProbFileName, tName);

	sprintf(tName, "%s.dist.%d", argv[(*i) + 1], IP->nRank);

	IP->Datafiles->dist_fpt = fopen(tName, "w");

	if (IP->Datafiles->dist_fpt == NULL) {
		sprintf(tmpstr, "Cannot open p(k|R) output file %s", tName);
		fprintf(stderr, "In strip_samplefile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif
		return FALSE;
	}
	fprintf(IP->Datafiles->dist_fpt,
		"seed iter seq k P(k|R)\n");

	NEW(IP->Datafiles->DistFileName, strlen(tName) + 1, char);
	strcpy(IP->Datafiles->DistFileName, tName);

	if (!IP->is_defined[cl_k])
		IP->nPostPlateauIter = POST_PLATEAU_SAMPLES;
	(*i)++;

	return TRUE;
}


short 
strip_bkgndfile(int argc, char **argv, IPtype IP, int *i, short bGui)
/* ============================================================== */
/* FUNCTION NAME : strip_bkgnd                                  */
/* */
/* DESCRIPTION : strips the name of the filecontaining          */
/* background composition data                    */
/****************************************************************/
{
	char            tmpstr[256];

	if (argc < (*i) + 2)
		p_usage_error("-B: Missing arguments");
	IP->Datafiles->bkgnd_fpt = fopen(argv[(*i) + 1], "r");

	if (IP->Datafiles->bkgnd_fpt == NULL) {
		sprintf(tmpstr, "Cannot open background composition file %s", argv[(*i) + 1]);
		fprintf(stderr, "In strip_bkgndfile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	IP->is_defined[cl_B] = TRUE;
	NEW(IP->Datafiles->BkgndFileName, strlen(argv[(*i) + 1]) + 1, char);	/* BT 1/27/97 */
	strcpy(IP->Datafiles->BkgndFileName, argv[(*i) + 1]);
	(*i)++;

	return TRUE;
}


short 
strip_spacingfile(int argc, char **argv, IPtype IP, int *i, short bGui)
/* ============================================================== */
/* FUNCTION NAME : strip_spacingfile                            */
/* */
/* DESCRIPTION : strips the name of the file containing the     */
/* spacing model                                  */
/****************************************************************/
{
	char            tmpstr[256];

	if (argc < (*i) + 2)
		p_usage_error("-U: MISSING ARGUMENTS");
	IP->Datafiles->spacing_fpt = fopen(argv[(*i) + 1], "r");

	if (IP->Datafiles->spacing_fpt == NULL) {
		sprintf(tmpstr, "Cannot open spacing model file %s", argv[(*i) + 1]);
		fprintf(stderr, "In strip_spacingfile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	IP->is_defined[cl_U] = TRUE;
	NEW(IP->Datafiles->SpacingFileName, strlen(argv[(*i) + 1]) + 1, char);	/* BT 1/27/97 */
	strcpy(IP->Datafiles->SpacingFileName, argv[(*i) + 1]);
	(*i)++;

	return TRUE;
}


short 
strip_weightfile(int argc, char **argv, IPtype IP, int *i, short bGui)
/* ============================================================== */
/* FUNCTION NAME : strip_weightfile                             */
/* */
/* DESCRIPTION : strips the name of the file containing         */
/* sequence weight data                           */
/****************************************************************/
{
	char            tmpstr[256];

	if (argc < (*i) + 2)
		p_usage_error("-H: MISSING ARGUMENTS");
	IP->Datafiles->weight_fpt = fopen(argv[(*i) + 1], "r");

	if (IP->Datafiles->weight_fpt == NULL) {
		sprintf(tmpstr, "Cannot open weight file %s", argv[(*i) + 1]);
		fprintf(stderr, "In strip_weightfile\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif

		return FALSE;
	}
	NEW(IP->Datafiles->WeightFileName, strlen(argv[(*i) + 1]) + 1, char);	/* BT 1/27/97 */
	strcpy(IP->Datafiles->WeightFileName, argv[(*i) + 1]);
	(*i)++;

	IP->is_defined[cl_H] = TRUE;

	return TRUE;
}


void 
strip_sample_iterations(int argc, char **argv, IPtype IP, int *i)
/* ============================================================= */
/* FUNCTION NAME : strip_sample_iterations                     */
/* */
/* DESCRIPTION : strips the maximum number of iterations used  */
/* for sampling after plateau                    */
/* ============================================================= */
{
	int             temp;

	if (argc < (*i) + 2)
		p_usage_error("-k: MISSING ARGUMENTS");
	temp = atoi(argv[(*i) + 1]);
	if (temp < 1)
		p_error("Iteration count must be at least 1");
	IP->nPostPlateauIter = temp;
	IP->is_defined[cl_k] = TRUE;
	(*i)++;
}


void 
strip_exchange(int argc, char **argv, Model B, int *pos)
/*
 *=============================================================
 * FUNCTION NAME : strip_exchange
 *
 * DESCRIPTION : strips information for exchange monte carlo
 *
 * MODIFICATIONS HISTORY:
 *
 * 6/30/00: New Temperature schedule,taking temp as random var
 *          currTemp=dMinTemp + drand()*dMaxTemp
 *=============================================================
 */
{
	int             temp;
	double         *nval;

	if (!B->IP->nMPI)
		p_error("The -X option requires the MPI version of Gibbs.");

	NEW(B->AN, 1, AnealStruct);

	B->AN->nExchangeIterations = EXCHANGE_ITERATIONS;
	B->AN->dMaxTemp = MAXTEMP;
	B->AN->dMinTemp = MINTEMP;
	B->AN->exchangePeriod = 0;
	B->AN->bExchange = TRUE;
	B->IP->is_defined[cl_X] = TRUE;

	if ((argc > (*pos) + 1) && isdigit((int) argv[(*pos) + 1][0])) {
		NEW(nval, 300, double);
		temp = ParseReals(argv[(*pos) + 1], nval);
		(*pos)++;

		B->AN->dMinTemp = nval[0];
		if (B->AN->dMinTemp <= 0)
			p_error("The minimum temperature must be > 0.");

		if (nval[1] != 0) {
			B->AN->dMaxTemp = nval[1];
			if (B->AN->dMaxTemp <= 0)
				p_error("The maximum temperature must be positive.");
		}
		if (nval[2] != 0 && nval[3] == 0) {
			B->AN->nExchangeIterations = nval[2];
			if (B->AN->nExchangeIterations <= 0)
				p_error("The number of steps must be >= 1.");
		}
		if (nval[3] != 0) {
			B->AN->exchangePeriod = nval[3];
			if (B->AN->exchangePeriod <= 0)
				p_error("The minimum number of milliseconds between exchanges must be > 0.");
		}
		free(nval);
	}
	printf("min temp = %f max temp = %f iters = %d period = %d\n",
	       B->AN->dMinTemp, B->AN->dMaxTemp, B->AN->nExchangeIterations, (int) B->AN->exchangePeriod);

	B->AN->currTemp = B->AN->dMaxTemp;
}


void 
strip_frag_width(int argc, char **argv, IPtype IP, int *pos)
/* ============================================================ */
/* FUNCTION NAME : strip_frag_width                           */
/* */
/* DESCRIPTION : strips off the value of the width multiplier */
/* for fragmentation                            */
/* ============================================================ */
{
	int             temp;
	int            *nval;
	int             width;
	int             t;
	int             i;

	NEW(nval, 300, int);
	if (argc < (*pos) + 2)
		p_error("-M: Missing fragmentation width arguments");
	temp = ParseIntegers(argv[(*pos) + 1], nval, "Getting fragmentation widths");
	if ((temp % 2) != 0)
		p_error("-M: Mismatched fragmentation begin and end ranges");
	for (i = 0; i < temp / 2; i++) {
		t = nval[i * 2 + 1] - 1;
		width = nval[i * 2 + 2];

		if (width < 0)
			p_error("-M: Negative fragmentation width");
		if (width < IP->nMotifLen[t])
			p_error("-M: Fragmentation widths must be greater than or equal to motif lengths");

		IP->nMaxFragWidth[t] = width;
	}

	(*pos)++;
	IP->is_defined[cl_M] = TRUE;
	free(nval);
}


void 
strip_ksamplemax(int argc, char **argv, Model B, int *i)
/* ============================================================ */
/* FUNCTION NAME : strip_ksmaplemax                           */
/* */
/* DESCRIPTION : strips off the value of the map value at     */
/* which to start sampling from map for k       */
/* ============================================================ */

{
	double          temp;
	char           *endptr;

	if (argc > (*i) + 1) {
		temp = strtod(argv[(*i) + 1], &endptr);
		if (strcmp(endptr, "") == 0) {
			B->IP->dKSampleMap = temp;
			(*i)++;
		}
	}
	B->IP->is_defined[cl_K] = TRUE;
}


void 
strip_homolog_seqs(int argc, char **argv, Model B, int *i)
/* ============================================================ */
/* FUNCTION NAME : strip_homolog_seqs                         */
/* */
/* DESCRIPTION : strips off the number of homologous seqs     */
/* ============================================================ */

{
	double          temp;
	int            *nval;

	NEW(B->Phylo->phyloSpecies, 1, int);

	B->Phylo->phyloSpecies[0] = 2;
	B->Phylo->treeCount = 0;

	if (argc > (*i) + 1 && isdigit((int) argv[(*i) + 1][0])) {
		NEW(nval, 300, int);
		temp = ParseIntegers(argv[(*i) + 1], nval, "Getting phylogenetic sequences");
		(*i)++;
		B->Phylo->phyloSpecies[0] = nval[1];
		if (B->Phylo->phyloSpecies[0] < 1)
			p_error("-D: Number of species must be >= 1.");
		if (temp > 1)
			B->Phylo->maxPhyloSeq = nval[2] - 1;
		free(nval);
	}
	B->IP->is_defined[cl_D] = TRUE;

	/*
	 * for(t = 0; t < B->IP->nNumMotifTypes; t++) { for( j = 0; j <
	 * B->IP->nNumMotifs[t][FORWARD] % B->Phylo->phyloSpecies[0]; j++ )
	 * B->IP->nNumMotifs[t][FORWARD]++; for( j = 0; j <
	 * B->IP->nNumMotifs[t][REVERSE] % B->Phylo->phyloSpecies[0]; j++ )
	 * B->IP->nNumMotifs[t][REVERSE]++; }
	 */
}


void 
strip_verify(int argc, char **argv, IPtype IP, int *i)
/* ============================================================ */
/* FUNCTION NAME : strip_verify                               */
/* */
/* DESCRIPTION : strips off the value of the maximum number of */
/* sequences to be used for verification        */
/* ============================================================ */

{
	int             temp;

	if (argc < (*i) + 2)
		p_usage_error("-V: Missing arguments");
	temp = atoi(argv[(*i) + 1]);
	if (temp <= 0)
		p_error("-V: Sequences to verify must be >= 0");
	IP->nVerifySeq = temp;
	IP->is_defined[cl_V] = TRUE;
	(*i)++;
}


void 
strip_nearopt_cutoff(int argc, char **argv, IPtype IP, int *i)
/* =========================================================== */
/* FUNCTION NAME : strip_nearopt_cutoff                      */
/* */
/* DESCRIPTION : strips off the value for the display cutoff */
/* used in near optimal sampling               */
/* =========================================================== */

{
	double          dtemp;

	if (argc < (*i) + 2)
		p_usage_error("-nopt_disp: MISSING ARGUMENTS");
	if (!ParseReals(argv[(*i) + 1], &dtemp))
		p_error("Error in forming value after -nopt_disp");
	IP->dNearOptCutoff = dtemp;
	if ((IP->dCutoff >= 1.0) || (IP->dCutoff < 0))
		p_error("Near Opt display cutoff must be between 0 and 1");
	IP->is_defined[cl_nopt_disp] = TRUE;
	(*i)++;
}


void 
strip_hierarchical(int argc, char **argv, Model B, int *i)
/* ============================================================ */
/* FUNCTION NAME : strip_hierarchical                         */
/* */
/* DESCRIPTION : send information to master about sites/seq   */
/* update pseudocounts                          */
/* ============================================================ */

{
	int             temp;
	int            *nval;
	char            host[256];

	B->IP->hier_iter = 10;
	if (argc > (*i) + 1 && isdigit((int) argv[(*i) + 1][0])) {
		NEW(nval, 300, int);
		temp = ParseIntegers(argv[(*i) + 1], nval,
			 "-hm: Getting iterations for hierarchical model.");
		(*i)++;
		B->IP->hier_iter = nval[1];
		if (B->IP->hier_iter < 1)
			p_error("-hm: The number of iterations must be >= 1.");
		free(nval);
	}
	gethostname(host, 256);
	NEW(B->IP->Datafiles->hierTempFileName, FILENAME_MAX, char);
	sprintf(B->IP->Datafiles->hierTempFileName, "hier.%s.%d.out", host, (int) getpid());

	B->IP->is_defined[cl_hm] = TRUE;
}


void 
strip_bayes(int argc, char **argv, Model B, int *i)
/* ============================================================ */
/* FUNCTION NAME : strip_bayes                                */
/* */
/* DESCRIPTION : set up Bayesian sampling                     */
/* ============================================================ */

{
	int             temp;
	int            *nval;

	B->IP->burnInPeriod = DEF_BURN_IN;
	B->IP->bayesSamplePeriod = DEF_BAYES_SAMPLE;

	if (argc > (*i) + 1 && isdigit((int) argv[(*i) + 1][0])) {
		NEW(nval, 300, int);
		temp = ParseIntegers(argv[(*i) + 1], nval,
		    "-bayes: Getting sample period for Bayesian sampling.");
		(*i)++;
		B->IP->burnInPeriod = nval[1];
		if (B->IP->burnInPeriod < 1)
			p_error("-bayes: The number of burn-in iterations must be >= 1.");

		if (temp > 1) {
			B->IP->bayesSamplePeriod = nval[2];
			if (B->IP->bayesSamplePeriod < 1)
				p_error("-bayes: The number of sample iterations must be >= 1.");
		}
		free(nval);
	}
	B->IP->is_defined[cl_bayes] = TRUE;
	B->IP->is_defined[cl_sample_model] = TRUE;
	B->IP->is_defined[cl_nopt] = TRUE;
	B->IP->is_defined[cl_opt] = FALSE;
	B->IP->is_defined[cl_y] = TRUE;
	B->IP->is_defined[cl_freq] = FALSE;
	B->IP->is_defined[cl_m] = TRUE;
	B->IP->is_defined[cl_max] = FALSE;
	B->IP->is_defined[cl_pred_update] = FALSE;
}


/************************* stripargs ***************************************/
/* */
/* This function strips out all of the command line arguments              */
/* that set up the program as desired                        */
/***************************************************************************/

short 
stripargs(int argc, char **argv, Model B, short bGui)
/* BT 3/21/97 */
{
	int             i, t, temp, *nval;
	char           *tmpstr;
	int             beg_val = 4;
	IPtype          IP;

	if (strcmp(argv[1], "-h") == 0) {
		print_usage_bernoulli();
		exit(0);
	}
	IP = B->IP;

	IP->argc = argc;
	IP->argv = argv;

	NEW(tmpstr, 80, char);
	NEW(nval, 300, int);	/* BT 12/3/99 */
	NEW(IP->Datafiles, 1, files);
	if (argc < 3)
		/* filename, motif length and */
		/* expected number are needed */
		p_usage_error("MISSING ARGUMENTS");

	for (i = 0; i < sizeof(IP->is_defined) / sizeof(short); i++)
		IP->is_defined[i] = FALSE;
	/* First argument is the input sequence file name */
	IP->Datafiles->fpt = fopen(argv[1], "r");
	if (IP->Datafiles->fpt == NULL) {
		sprintf(tmpstr, "Cannot open data file %s", argv[1]);
		fprintf(stderr, "In strip args\n");
#ifdef _GUI_
		if (bGui) {
			wid_error((Widget) toplevel, (Widget) 0, tmpstr, 0, 0);
			fprintf(stderr, tmpstr);
		} else
			p_error(tmpstr);
#else
		p_error(tmpstr);
#endif
		return FALSE;
	}
	IP->Datafiles->filename = (char *) malloc(strlen(argv[1]) * sizeof(char) + 1);	/* BT 1/24/97 */
	strcpy(IP->Datafiles->filename, argv[1]);

	IP->Datafiles->InputFilename = (char *) malloc(strlen(argv[1]) * sizeof(char) + 1);	/* BT 7/2/97 */
	strcpy(IP->Datafiles->InputFilename, argv[1]);

	/* second arg is the motif widths */
	IP->nNumMotifTypes =
		ParseIntegers(argv[2], nval, "Getting Motif Lengths");
	NEW(IP->dposterior_prob, IP->nNumMotifTypes, double);
	NEW(IP->nMotifLen, IP->nNumMotifTypes, int);
	NEW(IP->DOF, IP->nNumMotifTypes, int);
	for (t = 0; t < IP->nNumMotifTypes; t++) {
		if (nval[t + 1] <= 0)	/* BT 10/24/03 */
			p_error("stripargs - Getting Motif Lengths: motif lengths must be > 0");
		IP->nMotifLen[t] = nval[t + 1];
	}
	NEWP(IP->nNumMotifs, IP->nNumMotifTypes, int);
	/* each motif can be forward and backward */
	for (t = 0; t < IP->nNumMotifTypes; t++)
		NEW(IP->nNumMotifs[t], 2, int);

	NEW(IP->col_shift, IP->nNumMotifTypes, int);	/* BT 6/6/2001 */

	NEW(IP->nMaxFragWidth, IP->nNumMotifTypes, int);
	for (t = 0; t < IP->nNumMotifTypes; t++)
		IP->nMaxFragWidth[t] = ((double) IP->nMotifLen[t]) * MAX_WIDTH_MULT;

	NEW(IP->nMinMotifLen, IP->nNumMotifTypes, int);	/* BT 5/9/97 */
	NEW(IP->nMaxMotifLen, IP->nNumMotifTypes, int);
	NEW(IP->nInputMotifLen, IP->nNumMotifTypes, int);
	for (t = 0; t < IP->nNumMotifTypes; t++) {
		IP->nMinMotifLen[t] = MIN_MOTIF_LEN;
		IP->nMaxMotifLen[t] = 5 * IP->nMotifLen[t];
		IP->nInputMotifLen[t] = IP->nMotifLen[t];
	}

	NEWPP(IP->nWidthCnts, IP->nNumMotifTypes, double);	/* BT 2/18/98 */
	for (t = 0; t < IP->nNumMotifTypes; t++) {
		NEWP(IP->nWidthCnts[t], 2, double);
		NEW(IP->nWidthCnts[t][0], IP->nMaxMotifLen[t] + 1, double);
		NEW(IP->nWidthCnts[t][1], IP->nMaxMotifLen[t] + 1, double);
	}

	if ((argc == 3) || (argv[3][0] == '-')) {
		/* site sampler */
		beg_val = 3;
		IP->site_samp = TRUE;
	} else {
		/* motif sampler */
		beg_val = 4;
		IP->site_samp = FALSE;
		temp =
			ParseIntegers(argv[3], nval, "Getting Initial Motif Guesses");
		for (t = 0; t < IP->nNumMotifTypes; t++) {
			/* Initialize all motifs in forward direction */
			if (nval[t + 1] < 0)	/* BT 10/24/03 */
				p_error("stripargs - Getting InitialMotif Counts: motif counts must be >= 0");
			IP->nNumMotifs[t][FORWARD] = nval[t + 1];
			IP->nNumMotifs[t][REVERSE] = 0;
		}
	}
	/* default values for amino acids */
	IP->nAlphaLen = 20;
	IP->RevComplement = FALSE;

	/* allocate space for palandromic and collapsed alphabet in case  */
	/* of nucleoid acids                                              */
	NEW(IP->AltModel, 1, AltModelStruct);
	NEWP(IP->AltModel->Collapsed, IP->nNumMotifTypes, int);
	NEWP(IP->AltModel->Palandromic, IP->nNumMotifTypes, int);
	NEWP(IP->AltModel->Repeat, IP->nNumMotifTypes, int);
	NEW(IP->AltModel->Concen, 1, Concenstruct);
	NEWP(IP->AltModel->Concen->Concentrated, IP->nNumMotifTypes, short);
	NEW(IP->AltModel->Concen->NumConcen, IP->nNumMotifTypes, short);

	IP->nMaxBlocks = IP->nNumMotifTypes;

	for (t = 0; t < IP->nNumMotifTypes; t++) {
		NEW(IP->AltModel->Collapsed[t], IP->nMaxMotifLen[t], int);
		NEW(IP->AltModel->Palandromic[t], IP->nMaxMotifLen[t], int);
		NEW(IP->AltModel->Repeat[t], IP->nMaxMotifLen[t], int);
		NEW(IP->AltModel->Concen->Concentrated[t], IP->nMaxMotifLen[t], short);
		for (i = 0; i < IP->nMaxMotifLen[t]; i++) {
			IP->AltModel->Collapsed[t][i] = FALSE;
			IP->AltModel->Palandromic[t][i] = FALSE;
			IP->AltModel->Repeat[t][i] = FALSE;
		}
	}

	/* print defaults */
	IP->is_defined[cl_opt] = TRUE;
	IP->is_defined[cl_freq] = TRUE;
	IP->is_defined[cl_max] = TRUE;
	IP->is_defined[cl_pred_update] = TRUE;

	IP->ticks = sysconf(_SC_CLK_TCK);

	/* begin to find optional arguments, the first fews are required */
	for (i = beg_val; i < argc; i++)
		if (strcmp(argv[i], "-r") == 0) {
			/* turn off reverse complement */
			IP->RevComplement = FALSE;
			IP->is_defined[cl_r] = TRUE;
		} else if (strcmp(argv[i], "-R") == 0)
			/* Palandromic Sequences */
			strip_palandromes(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-I") == 0)
			/* Repeat Sequences */
			strip_repeats(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-x") == 0)
			/* DON'T remove low complexity regions */
			IP->is_defined[cl_x] = TRUE;
		else if (strcmp(argv[i], "-c") == 0)
			/* collapse alphabet */
			strip_collapsed(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-a") == 0)
			/* Concentrated regions */
			strip_concen(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-C") == 0)
			/* Near opt. cutoff */
			strip_cutoff(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-n") == 0) {
			/* default values for Nucleic Acid Alphabet */
			IP->nAlphaLen = 4;
			/*
			 * make sure -r option take effect anywhere it
			 * appears
			 */
			if (IP->is_defined[cl_r])
				IP->RevComplement = FALSE;
			else
				IP->RevComplement = TRUE;
		} else if (strcmp(argv[i], "-p") == 0)
			/* PLATEAU PERIODS */
			strip_plateau(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-j") == 0)
			/* ADJUST PERIODS */
			strip_adjust(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-s") == 0)
			/* initial SEED VALUE */
			strip_seedval(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-S") == 0)
			/* NUMBER OF SEEDS */
			strip_seednum(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-W") == 0)
			/* PSEUDO SITE WEIGHT */
			strip_sitewt(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-w") == 0)
			/* PSEUDO COUNT WEIGHT */
			strip_countwt(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-i") == 0)
			/* MAX ITERATIONS */
			strip_iterations(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-o") == 0) {
			/* OUTPUT FILE DEFINED */
			if (!strip_outfile(argc, argv, IP, &i, bGui))	/* BT 4/25/97 */
				return FALSE;
		} else if (strcmp(argv[i], "-N") == 0) {
			/* SCAN FILE DEFINED */
			if (!strip_scanfile(argc, argv, IP, &i, bGui))	/* BT 9/17/97 */
				return FALSE;
		} else if (strcmp(argv[i], "-t") == 0)
			/* display sites for nearopt */
			IP->is_defined[cl_t] = TRUE;
		else if (strcmp(argv[i], "-m") == 0) {
			/* turn off maximization sampling */
			IP->is_defined[cl_m] = TRUE;
			IP->is_defined[cl_max] = FALSE;
		} else if (strcmp(argv[i], "-e") == 0)
			/* turn on EM sampling */
			/*
			 * IP->is_defined[cl_e] = TRUE; *//* BT 09/22/04
			 * remove em option
			 */
			p_error("-e option is deprecated.");
		else if (strcmp(argv[i], "-F") == 0)
			/* Turn off fragmentation */
			IP->is_defined[cl_F] = TRUE;
		else if (strcmp(argv[i], "-d") == 0) {
			/* Vary Width */
			strip_width(argc, argv, IP, &i);
			for (t = 0; t < IP->nNumMotifTypes; t++) {	/* realloc alphabet
									 * arrays */
				int             oldLen;
				oldLen = IP->nMaxMotifLen[t];
				IP->AltModel->Collapsed[t] =
					(int *) realloc(IP->AltModel->Collapsed[t],
					 IP->nMaxMotifLen[t] * sizeof(int));
				IP->AltModel->Palandromic[t] =
					(int *) realloc(IP->AltModel->Palandromic[t],
					 IP->nMaxMotifLen[t] * sizeof(int));
				IP->AltModel->Repeat[t] =
					(int *) realloc(IP->AltModel->Repeat[t],
					 IP->nMaxMotifLen[t] * sizeof(int));
				IP->AltModel->Concen->Concentrated[t] =
					(short *) realloc(IP->AltModel->Concen->Concentrated[t],
				       IP->nMaxMotifLen[t] * sizeof(short));
				if (IP->nMaxMotifLen[t] > oldLen) {
					for (i = oldLen + 1; i < IP->nMaxMotifLen[t]; i++) {
						IP->AltModel->Collapsed[t][i] = FALSE;
						IP->AltModel->Palandromic[t][i] = FALSE;
						IP->AltModel->Repeat[t][i] = FALSE;
					}
				}
			}
		} else if (strcmp(argv[i], "-u") == 0)
			/* Suboptimal sampler output */
			IP->is_defined[cl_u] = TRUE;	/* BT 3/19/97 */
		else if (strcmp(argv[i], "-v") == 0)	/* BT 4/23/97 */
			/* Strip overlap paramater */
			strip_overlap(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-l") == 0) {	/* BT 4/23/97 */
			/* Wilcoxon test */
			if (!strip_wilcoxon(argc, argv, IP, &i, bGui))	/* BT 1/28/98 */
				return FALSE;
		}
	/* sample along length */
		else if (strcmp(argv[i], "-g") == 0)	/* BT 9/10/97 */
			IP->is_defined[cl_g] = TRUE;
		else if (strcmp(argv[i], "-h") == 0) {
			/* HELP MODE */
			print_usage_bernoulli();
			exit(0);
		} else if (strcmp(argv[i], "-P") == 0) {
			/* PRIORS FILE DEFINED */
			if (!strip_priorfile(argc, argv, IP, &i, bGui))	/* BT 4/25/97 */
				return FALSE;
		} else if (strcmp(argv[i], "-O") == 0) {
			/* PRIORS OUTPUT FILE DEFINED */
			if (!strip_proutfile(argc, argv, IP, &i, bGui))	/* BT 4/25/97 */
				return FALSE;
		}
#ifndef _PUBLIC_
	/* Recursuive sampling */
		else if (strcmp(argv[i], "-E") == 0)	/* BT 2/11/97 */
			strip_maxblocks(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-f") == 0)	/* BT 2/11/97 */
			strip_rcutoff(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-z") == 0)	/* BT 6/10/98 */
			IP->is_defined[cl_z] = TRUE;
		else if (strcmp(argv[i], "-Q") == 0) {
			if (!strip_samplefile(argc, argv, IP, &i, bGui))
				return FALSE;
		} else if (strcmp(argv[i], "-k") == 0)
			strip_sample_iterations(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-b") == 0)
			IP->is_defined[cl_b] = TRUE;
		else if (strcmp(argv[i], "-B") == 0) {
			if (!strip_bkgndfile(argc, argv, IP, &i, bGui))
				return FALSE;
		} else if (strcmp(argv[i], "-U") == 0) {
			if (!strip_spacingfile(argc, argv, IP, &i, bGui))
				return FALSE;
		} else if (strcmp(argv[i], "-A") == 0)
			IP->is_defined[cl_A] = TRUE;
#ifndef _PUBLIC2_
		else if (strcmp(argv[i], "-H") == 0) {
			if (!strip_weightfile(argc, argv, IP, &i, bGui))
				return FALSE;
		} else if (strcmp(argv[i], "-G") == 0)
			IP->is_defined[cl_G] = TRUE;
#endif				/* _PUBLIC2_ */
		else if (strcmp(argv[i], "-D") == 0)
			strip_homolog_seqs(argc, argv, B, &i);
		else if (strcmp(argv[i], "-X") == 0)	/* BT 3/30/99 */
			strip_exchange(argc, argv, B, &i);
		else if (strcmp(argv[i], "-V") == 0)
			strip_verify(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-K") == 0)
			strip_ksamplemax(argc, argv, B, &i);
		else if (strcmp(argv[i], "-y") == 0)
			IP->is_defined[cl_y] = TRUE;
		else if (strcmp(argv[i], "-Y") == 0)
			IP->is_defined[cl_Y] = TRUE;
		else if (strcmp(argv[i], "-q") == 0)	/* BT 2/20/97 */
			strip_samplecount(argc, argv, IP, &i);
#endif				/* _PUBLIC_ */
		else if (strcmp(argv[i], "-M") == 0)
			strip_frag_width(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-Z") == 0)
			IP->is_defined[cl_Z] = TRUE;
		else if (strcmp(argv[i], "-J") == 0) {	/* BT 9/7/2001 */
			if (IP->is_defined[cl_F])
				p_error("Remove -F option when using -J option.");
			IP->is_defined[cl_J] = TRUE;
		} else if (strcmp(argv[i], "-L") == 0)
			IP->is_defined[cl_L] = TRUE;
		else if (strcmp(argv[i], "-frag") == 0)
			IP->is_defined[cl_frag] = TRUE;
		else if (strcmp(argv[i], "-nopt") == 0)
			IP->is_defined[cl_nopt] = TRUE;
		else if (strcmp(argv[i], "-nopt_disp") == 0)
			strip_nearopt_cutoff(argc, argv, IP, &i);
		else if (strcmp(argv[i], "-sample_model") == 0)
			IP->is_defined[cl_sample_model] = TRUE;
		else if (strcmp(argv[i], "-hierarchical_model") == 0 ||
			 strcmp(argv[i], "-hm") == 0)
			strip_hierarchical(argc, argv, B, &i);
		else if (strcmp(argv[i], "-wilcox") == 0)
			IP->is_defined[cl_wilcox] = TRUE;
		else if (strcmp(argv[i], "-bayes") == 0)
			strip_bayes(argc, argv, B, &i);
		else if (strcmp(argv[i], "-opt") == 0) {
			IP->is_defined[cl_opt] = TRUE;
			IP->is_defined[cl_nopt] = FALSE;
		} else if (strcmp(argv[i], "-freq") == 0) {
			IP->is_defined[cl_freq] = TRUE;
			IP->is_defined[cl_y] = FALSE;
		} else if (strcmp(argv[i], "-max") == 0) {
			IP->is_defined[cl_max] = TRUE;
			IP->is_defined[cl_m] = FALSE;
		} else if (strcmp(argv[i], "-bayes_counts") == 0)
			IP->is_defined[cl_bayes_counts] = TRUE;
		else if (strcmp(argv[i], "-pred_update") == 0) {
			B->IP->is_defined[cl_pred_update] = TRUE;
			B->IP->is_defined[cl_sample_model] = FALSE;
		} else if (strcmp(argv[i], "-align_centroid") == 0)
			IP->is_defined[cl_align_centroid] = TRUE;
		else if (strcmp(argv[i], "-no_cred") == 0)
			IP->is_defined[cl_no_cred] = TRUE;
		else {
			sprintf(tmpstr, "unrecognized command argument: %s\n\n", argv[i]);
			printargs(argc, argv, IP, stdout);	/* BT 1/9/98 */
			printf("%s", tmpstr);
			p_usage_error(tmpstr);
		}
	free(tmpstr);
	free(nval);

	if (IP->is_defined[cl_X])
		NEW(B->AN->results, B->IP->nSeeds + 1, MaxResults);

#ifndef _MPI_
	printargs(argc, argv, IP, stdout);	/* BT 1/9/98 */
#endif

	return TRUE;		/* BT 3/21/97 */
}


void 
printargs(int argc, char **argv, IPtype IP, FILE * tmpFpt)
{				/* BT 1/9/98 */
	FILE           *fpt;
	int             i;

	fpt = IP->Datafiles->out_fpt;

	if (!IP->is_defined[cl_Z])
		fprintf(tmpFpt, "%s ", argv[0]);
	if (fpt != NULL)
		fprintf(fpt, "%s ", argv[0]);

	for (i = 1; i < argc; i++) {
		if (!IP->is_defined[cl_Z])
			fprintf(tmpFpt, "%s ", argv[i]);
		if (fpt != NULL)
			fprintf(fpt, "%s ", argv[i]);
	}

	if (!IP->is_defined[cl_Z])
		fprintf(tmpFpt, "\n");
	if (fpt != NULL)
		fprintf(fpt, "\n\n");
}


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


void 
print_options(Model B)
{				/* BT 1/27/97 */
	FILE           *fpt;
	int             t;
	IPtype          IP;
	char            currDir[PATH_MAX + 1];

	IP = B->IP;

	fpt = IP->Datafiles->out_fpt;

#ifdef _DEBUG_
	fprintf(fpt, "\nGibbs %d.%02d.%03d  %s Debugging Version\n", glVerMajor, glVerMinor, glRev, glRevDate);
#elif defined _MPI_
	fprintf(fpt, "\nGibbs %d.%02d.%03d  %s MPI Version\n", glVerMajor, glVerMinor, glRev, glRevDate);
#else
	fprintf(fpt, "\nGibbs %d.%02d.%03d  %s\n", glVerMajor, glVerMinor, glRev, glRevDate);
#endif

	if (IP->comment != NULL) {
		fprintf(fpt, "%s", IP->comment);
		fprintf(fpt, "\n");
	}
	fprintf(fpt, "Data file: %s\n", IP->Datafiles->InputFilename);
	if (IP->Datafiles->output_filename != NULL)
		fprintf(fpt, "Output file: %s\n", IP->Datafiles->output_filename);
	if (IP->Datafiles->prior_filename != NULL)
		fprintf(fpt, "Priors file: %s\n", IP->Datafiles->prior_filename);
	if (IP->Datafiles->ScanFilename != NULL)
		fprintf(fpt, "Scan file: %s\n", IP->Datafiles->ScanFilename);
	if (IP->Datafiles->PriorOutFilename != NULL)
		fprintf(fpt, "Prior output file: %s\n", IP->Datafiles->PriorOutFilename);
	if (IP->Datafiles->ControlFileName != NULL)
		fprintf(fpt, "Wilcoxon control file: %s\n", IP->Datafiles->ControlFileName);
	if (IP->Datafiles->SampleFileName != NULL)
		fprintf(fpt, "Sample counts file: %s\n", IP->Datafiles->SampleFileName);
	if (IP->Datafiles->BkgndFileName != NULL)
		fprintf(fpt, "Background Composition Model file: %s\n", IP->Datafiles->BkgndFileName);
	if (IP->Datafiles->SpacingFileName != NULL)
		fprintf(fpt, "Spacing Model file: %s\n", IP->Datafiles->SpacingFileName);
	if (IP->Datafiles->WeightFileName != NULL)
		fprintf(fpt, "SequenceWeight file: %s\n", IP->Datafiles->WeightFileName);

	getcwd(currDir, PATH_MAX * sizeof(char) + 1);
	fprintf(fpt, "Current directory: %s\n", currDir);

	fprintf(fpt, "The following options are set:\n");

	fprintf(fpt, "Concentrated Region          %5s    Sequence type        %5s\n",
		IP->is_defined[cl_a] ? "True " : "False",
		IP->is_defined[cl_t] ? "True " : "False");
	fprintf(fpt, "Collapsed Alphabet           %5s    Pseudocount weight   %5s\n",
		IP->is_defined[cl_c] ? "True " : "False",
		IP->is_defined[cl_w] ? "True " : "False");
	fprintf(fpt, "Use Expectation/Maximization %5s    Don't Xnu sequence   %5s\n",
		IP->is_defined[cl_e] ? "True " : "False",
		IP->is_defined[cl_x] ? "True " : "False");
	fprintf(fpt, "Help flag                    %5s    Near optimal cutoff  %5s\n",
		IP->is_defined[cl_h] ? "True " : "False",
		IP->is_defined[cl_C] ? "True " : "False");
	fprintf(fpt, "Number of iterations         %5s    Don't fragment       %5s\n",
		IP->is_defined[cl_i] ? "True " : "False",
		IP->is_defined[cl_F] ? "True " : "False");
	fprintf(fpt, "Don't use map maximization   %5s    Repeat regions       %5s\n",
		IP->is_defined[cl_m] ? "True " : "False",
		IP->is_defined[cl_I] ? "True " : "False");
	fprintf(fpt, "Output file                  %5s    Informed priors file %5s\n",
		IP->is_defined[cl_o] ? "True " : "False",
		IP->is_defined[cl_P] ? "True " : "False");
	fprintf(fpt, "Plateau periods              %5s    palindromic sequence %5s\n",
		IP->is_defined[cl_p] ? "True " : "False",
		IP->is_defined[cl_R] ? "True " : "False");
	fprintf(fpt, "Don't Reverse complement     %5s    Number of seeds      %5s\n",
		IP->is_defined[cl_r] ? "True " : "False",
		IP->is_defined[cl_S] ? "True " : "False");
	fprintf(fpt, "Seed Value                   %5s    Pseudosite weight    %5s\n",
		IP->is_defined[cl_s] ? "True " : "False",
		IP->is_defined[cl_W] ? "True " : "False");
	fprintf(fpt, "Suboptimal sampler output    %5s    Overlap              %5s\n",
		IP->is_defined[cl_u] ? "True " : "False",
		IP->is_defined[cl_v] ? "True " : "False");	/* BT 4/23/97 */
	fprintf(fpt, "Allow width to vary          %5s    Wilcoxon signed rank %5s\n",
		IP->is_defined[cl_d] ? "True " : "False",
		IP->is_defined[cl_l] ? "True " : "False");	/* BT 4/23/97 */
	fprintf(fpt, "Sample along length          %5s    Output Scan File     %5s\n",
		IP->is_defined[cl_g] ? "True " : "False",
		IP->is_defined[cl_N] ? "True " : "False");	/* BT 9/10/97 */
	fprintf(fpt, "Output prior file            %5s    Modular Sampler      %5s\n",
		IP->is_defined[cl_O] ? "True " : "False",
		IP->is_defined[cl_E] ? "True " : "False");
	fprintf(fpt, "Ignore Spacing Model         %5s    Sample Background    %5s\n",
		IP->is_defined[cl_z] ? "True " : "False",
		IP->is_defined[cl_b] ? "True " : "False");
	fprintf(fpt, "Bkgnd Comp Model             %5s    Init from prior      %5s\n",
		IP->is_defined[cl_B] ? "True " : "False",
		IP->is_defined[cl_A] ? "True " : "False");
	fprintf(fpt, "Homologous Seq pairs         %5s    Parallel Tempering   %5s\n",
		IP->is_defined[cl_D] ? "True " : "False",
		IP->is_defined[cl_X] ? "True " : "False");
	fprintf(fpt, "Group Sampler                %5s    No progress info     %5s\n",
		IP->is_defined[cl_G] ? "True " : "False",
		IP->is_defined[cl_Z] ? "True " : "False");
	fprintf(fpt, "Fragment from middle         %5s    Verify Mode          %5s\n",
		IP->is_defined[cl_J] ? "True " : "False",
		IP->is_defined[cl_V] ? "True " : "False");
	fprintf(fpt, "Alternate sample on k        %5s    No freq. soln.       %5s\n",
		IP->is_defined[cl_K] ? "True " : "False",
		IP->is_defined[cl_y] ? "True " : "False");
	fprintf(fpt, "Calc. def. pseudo wt.        %5s    Motif/Recur smpl     %5s\n",
		IP->is_defined[cl_Y] ? "True " : "False",
		IP->is_defined[cl_L] ? "True " : "False");
	fprintf(fpt, "Phylogenetic Sampling        %5s    Supress Near Opt.    %5s\n",
		B->Phylo->phyloTree ? "True " : "False",
		IP->is_defined[cl_nopt] ? "True " : "False");
	fprintf(fpt, "Nearopt display cutoff       %5s    Sample model         %5s\n",
		IP->is_defined[cl_nopt_disp] ? "True " : "False",
		IP->is_defined[cl_sample_model] ? "True " : "False");
	fprintf(fpt, "Hierarchical Model           %5s    Centroid model       %5s\n",
		IP->is_defined[cl_hm] ? "True " : "False",
		IP->is_defined[cl_bayes] ? "True " : "False");
	fprintf(fpt, "Print Bayesian Counts        %5s    Align Centroid       %5s\n",
		IP->is_defined[cl_bayes_counts] ? "True " : "False",
		IP->is_defined[cl_align_centroid] ? "True " : "False");
	fprintf(fpt, "Calculate Credibility Limits %5s    Frequency Bkgnd.     %5s\n",
		IP->is_defined[cl_bayes] && !IP->is_defined[cl_no_cred] ? "True " : "False",
		IP->is_defined[cl_freq_background] ? "True " : "False");
	fprintf(fpt, "\n");

	fprintf(fpt, "site_samp            = %12d\n", IP->site_samp);
	fprintf(fpt, "nMotifLen            = ");
	for (t = 0; t < IP->nNumMotifTypes; t++) {	/* 4/25/97 */
		if (t > 0) {
			fprintf(fpt, ", ");
			fprintf(fpt, "%d", IP->nMotifLen[t]);
		} else
			fprintf(fpt, "%12d", IP->nMotifLen[t]);
	}
	fprintf(fpt, "\n");
	fprintf(fpt, "nAlphaLen            = %12d\n", IP->nAlphaLen);
	fprintf(fpt, "nNumMotifs           = ");
	for (t = 0; t < IP->nNumMotifTypes; t++) {	/* 4/25/97 */
		if (t > 0) {
			fprintf(fpt, " ,");
			fprintf(fpt, "%d", IP->nNumMotifs[t][FORWARD] + IP->nNumMotifs[t][REVERSE]);
		} else
			fprintf(fpt, "%12d", IP->nNumMotifs[t][FORWARD] + IP->nNumMotifs[t][REVERSE]);

	}
	fprintf(fpt, "\n");
	fprintf(fpt, "dPseudoCntWt         = %12.6g\n", IP->dPseudoCntWt);
	fprintf(fpt, "dPseudoSiteWt        = %12.6g\n", IP->dPseudoSiteWt);
	if (!IP->is_defined[cl_bayes])
		fprintf(fpt, "nMaxIterations       = %12d\n", IP->nMaxIterations);
	fprintf(fpt, "lSeedVal             = %12ld\n", IP->lSeedVal);
#ifdef _MPI_
	fprintf(fpt, "lSeedVal0            = %12ld\n", IP->lSeedVal0);	/* BT 02/25/03 */
#endif
	fprintf(fpt, "nPlateauPeriods      = %12d\n", IP->nPlateauPeriods);
	fprintf(fpt, "nSeeds               = %12d\n", IP->nSeeds);
	fprintf(fpt, "nNumMotifTypes       = %12d\n", IP->nNumMotifTypes);
	fprintf(fpt, "dCutoff              = %12.6g\n", IP->dCutoff);
	fprintf(fpt, "dNearOptDispCutoff   = %12.6g\n", IP->dNearOptCutoff);
	fprintf(fpt, "RevComplement        = %12d\n", IP->RevComplement);
	fprintf(fpt, "glOverlapParam       = %12.6g\n", IP->glOverlapParam);	/* BT 4/21/97 */
	if (IP->is_defined[cl_q])
		fprintf(fpt, "Post. width cnt      = %12d\n", IP->nSampleCnt);	/* BT 2/20/98 */
	if (IP->is_defined[cl_d]) {
		fprintf(fpt, "Min Motif Width     = ");
		for (t = 0; t < IP->nNumMotifTypes; t++) {	/* 4/25/97 */
			if (t > 0)
				fprintf(fpt, ",%d", IP->nMinMotifLen[t]);
			else
				fprintf(fpt, "%12d", IP->nMinMotifLen[t]);
		}
		fprintf(fpt, "\n");
		fprintf(fpt, "Max Motif Width      = ");
		for (t = 0; t < IP->nNumMotifTypes; t++) {	/* 4/25/97 */
			if (t > 0)
				fprintf(fpt, ",%d", IP->nMaxMotifLen[t]);
			else
				fprintf(fpt, "%12d", IP->nMaxMotifLen[t]);
		}
		fprintf(fpt, "\n");
	}
	if (IP->is_defined[cl_E]) {
		fprintf(fpt, "Max Sites/seq        = %12d\n", IP->nMaxBlocks);
		fprintf(fpt, "Min Sites/Seq        = %12d\n", IP->nMinBlocks);
	}
	fprintf(fpt, "Rcutoff factor       = %12g\n", IP->dRCutoff);
	fprintf(fpt, "Post Plateau Samples = %12d\n", IP->nPostPlateauIter);
	if (IP->is_defined[cl_X]) {
		fprintf(fpt, "Min Temp             = %12g\n", B->AN->dMinTemp);
		fprintf(fpt, "Max Temp             = %12g\n", B->AN->dMaxTemp);
	}
	fprintf(fpt, "Frag/Shft Per.       = %12d\n", IP->nAdjustPeriod);
	if (!IP->is_defined[cl_F]) {
		fprintf(fpt, "Frag width           = ");
		for (t = 0; t < IP->nNumMotifTypes; t++) {	/* 4/25/97 */
			if (t > 0)
				fprintf(fpt, ",%d", B->F->nMaxLen[t]);
			else
				fprintf(fpt, "%12d", B->F->nMaxLen[t]);
		}
		fprintf(fpt, "\n");
	}
	if (IP->is_defined[cl_bayes]) {
		fprintf(fpt, "Burn-in              = %12d\n", IP->burnInPeriod);
		fprintf(fpt, "Sample Period        = %12d\n", IP->bayesSamplePeriod);
	}
	fprintf(fpt, "\n\n");

	for (t = 0; t < IP->nNumMotifTypes; t++) {
		if (IP->nMotifLen[t] == 1 && !IP->is_defined[cl_F])
			p_error("Motif models of length = 1 are not compatible with fragmentation.");
	}
}


/******************* get_input ***********************************/
/* */
/* DESCRIPTION : This function sets the default values           */
/* for parameters that have not been               */
/* entered in through the command line             */
/* arguments                                       */
/*****************************************************************/

void 
get_inputs(Model B)
{

	int             t, j, i;
	int             tot;
	int             nFor;

	B->IP->nNumProcessed = 0;
	if (!B->IP->is_defined[cl_W])
		B->IP->dPseudoSiteWt = PSEUDO_SITE_WT;
	if (!B->IP->is_defined[cl_w])
		B->IP->dPseudoCntWt = PSEUDO_CNT_WT;
	if (!B->IP->is_defined[cl_p])
		B->IP->nPlateauPeriods = PLATEAU_PERIODS;
	if (!B->IP->is_defined[cl_i])
		B->IP->nMaxIterations = MAX_ITERATIONS;
	if (!B->IP->is_defined[cl_S])
		B->IP->nSeeds = NUM_SEEDS;
	if (!B->IP->is_defined[cl_C])
		B->IP->dCutoff = NEAROPT_CUTOFF_PERCNT;
	if (!B->IP->is_defined[cl_nopt_disp])
		B->IP->dNearOptCutoff = NEAROPT_CUT;
	if (!B->IP->is_defined[cl_o])
		B->IP->Datafiles->out_fpt = stdout;
	if (!B->IP->is_defined[cl_s])
		B->IP->lSeedVal = (long) time(NULL);
	if (!B->IP->is_defined[cl_v])
		B->IP->glOverlapParam = 0.0;	/* BT 4/23/97 */
	if (!B->IP->is_defined[cl_f])
		B->IP->dRCutoff = RCUTOFF_FACTOR;	/* BT 2/11/97 */
	if (!B->IP->is_defined[cl_q])
		B->IP->nSampleCnt = 0;	/* BT 2/20/98 */
	if (!B->IP->is_defined[cl_j])
		B->IP->nAdjustPeriod = 5;	/* BT 9/11/98 */
	if (!B->IP->is_defined[cl_nopt])
		B->IP->is_defined[cl_opt] = TRUE;
	if (!B->IP->is_defined[cl_y])
		B->IP->is_defined[cl_freq] = TRUE;
	if (!B->IP->is_defined[cl_m])
		B->IP->is_defined[cl_max] = TRUE;

	if (B->IP->RevComplement) {
		/* Update number  of motifs for reverse comp. */
		for (t = 0; t < B->IP->nNumMotifTypes; t++) {
			tot = B->IP->nNumMotifs[t][FORWARD];
			nFor = tot / 2;
			if (B->IP->is_defined[cl_D])
				nFor = 0;
			B->IP->nNumMotifs[t][REVERSE] = nFor;
			B->IP->nNumMotifs[t][FORWARD] = tot - nFor;
		}
	}
	for (t = 0; t < B->IP->nNumMotifTypes; t++) {
		B->IP->DOF[t] = B->IP->nMotifLen[t] * (B->IP->nAlphaLen - 1);
		for (i = 0; i < B->IP->nMotifLen[t]; i++) {
			if (B->IP->AltModel->Collapsed[t][i])
				B->IP->DOF[t] -= (B->IP->nAlphaLen) / 2.0 - 1.0;
			if ((B->IP->AltModel->Palandromic[t][i]) &&
			    (i < (B->IP->nMotifLen[t] / 2 + B->IP->nMotifLen[t] % 2)))
				B->IP->DOF[t] -= B->IP->nAlphaLen - 1.0;
			if ((B->IP->AltModel->Repeat[t][i]) &&
			    (i < (B->IP->nMotifLen[t] / 2 + B->IP->nMotifLen[t] % 2)))
				B->IP->DOF[t] -= B->IP->nAlphaLen - 1.0;
		}
		if (B->IP->is_defined[cl_a]) {
			/* DO CONCENTRATED REGIONS */
		}
		/* ???????????? */
		/** NOW NEED TO SET FOR VARIOUS MODELS **/
	}

	if (!B->IP->is_defined[cl_F]) {	/* fragmentation */
		NEW(B->F, 1, FragStruct);
		NEWP(B->F->nColMask, B->IP->nNumMotifTypes, int);
		NEWP(B->F->nOldColMask, B->IP->nNumMotifTypes, int);
		NEW(B->F->nMaxLen, B->IP->nNumMotifTypes, int);
		NEW(B->F->FragWidth, B->IP->nNumMotifTypes, int);
		NEW(B->F->nOldFragWidth, B->IP->nNumMotifTypes, int);
		NEW(B->F->shift, B->IP->nNumMotifTypes, int);
		NEWP(B->F->fragPos, B->IP->nNumMotifTypes, int);
		NEWP(B->F->fragInitMask, B->IP->nNumMotifTypes, int);
		/*
		 * nvFragCnts is initialized in setFragCnt - we do this to
		 * prevent
		 */
		/* a problem on the second run */
		B->F->nvFragCnts = NULL;	/* BT 2/7/97 */

		for (t = 0; t < B->IP->nNumMotifTypes; t++) {
			B->F->nMaxLen[t] = B->IP->nMaxFragWidth[t];
			B->F->FragWidth[t] = B->IP->nMotifLen[t];
			B->F->shift[t] = 0;
			NEW(B->F->nColMask[t], B->F->nMaxLen[t], int);
			NEW(B->F->nOldColMask[t], B->F->nMaxLen[t], int);
			NEW(B->F->fragPos[t], B->IP->nMotifLen[t], int);
			NEW(B->F->fragInitMask[t], B->F->nMaxLen[t], int);
			i = (B->F->nMaxLen[t] - B->IP->nMotifLen[t]) / 2;
			for (j = 0; j < i; j++)
				B->F->nColMask[t][j] = 0;
			for (j = i; j < i + B->IP->nMotifLen[t]; j++)
				B->F->nColMask[t][j] = 1;
			for (j = 0; j < B->IP->nMotifLen[t]; j++)
				B->F->fragPos[t][j] = 0;

			i = (B->F->nMaxLen[t] - B->IP->nMotifLen[t]) / 2;	/* BT 07/26/2001 */
			for (j = 0; j < i; j++)
				B->F->fragInitMask[t][j] = 0;
			for (j = i; j < i + B->IP->nMotifLen[t]; j++)
				B->F->fragInitMask[t][j] = 1;
			for (j = i + B->IP->nMotifLen[t]; j < B->F->nMaxLen[t]; j++)
				B->F->fragInitMask[t][j] = 0;
		}
	} else			/* do not fragment */
		B->F = NULL;

#ifndef _MPI_
	/* print_options( B ); *//* BT 1/27/97 */
#endif
}


short 
ParseReals(char *str, double *values)
{
	long            n;
	double          k;

	if (!isdigit((int) str[0]))
		return FALSE;
	for (n = 1; str[0] != 0;) {
		if (str[0] == ',') {
			n++;
			str++;
		} else if (isdigit((int) str[0])) {
			if (sscanf(str, "%lf", &k) != 1) {
				return FALSE;
			} else {
				values[n - 1] = k;
				while (isdigit((int) str[0]) || str[0] == '.')
					str++;
			}
		} else
			return FALSE;
	}
	return TRUE;
}


int 
ParseIntegers(char *str, int *values, char *msg)
{
	int             n, k;

	if (!isdigit((int) str[0]))
		p_error("Motif lengths must be integers > 0");
	for (n = 1; str[0] != 0;) {
		if (str[0] == ',') {
			n++;
			str++;
		} else if (isdigit((int) str[0])) {
			if (sscanf(str, "%d", &k) != 1) {
				p_error(msg);
			} else {
				values[n] = k;
				while (isdigit((int) str[0]))
					str++;
			}
		} else
			p_error(msg);
	}
	return n;
}


void 
PrintTempOut(FILE * fpt, char *fmt,...)
{
	va_list         argp;

	va_start(argp, fmt);
	vfprintf(fpt, fmt, argp);
	va_end(argp);
	fflush(fpt);
#ifdef _MPI_
	fsync(fileno(fpt));
#endif
}


void 
PrintSeqDescriptions(Model M)
{
	int             i = 1;
	int             newlen;

	PrintTempOut(M->IP->Datafiles->out_fpt, "Sequences to be Searched:\n");
	PrintTempOut(M->IP->Datafiles->out_fpt, "_________________________\n");

	for (i = 0; i < M->IP->nNumSequences; i++) {
		PrintTempOut(M->IP->Datafiles->out_fpt, "#%-3d ", i + 1);
		PrintTempOut(M->IP->Datafiles->out_fpt, "%s\n", M->IP->fastaHeader[i]);
	}

	if (M->IP->nAlphaLen == 20)
		newlen = M->IP->nSeqLen - M->IP->nNumProcessed;
	else
		newlen = M->First->nTotBack;
	PrintTempOut(M->IP->Datafiles->out_fpt,
		"Processed Sequence Length: %d Total sequence length: %d\n",
		     newlen, M->IP->nSeqLen);
}
