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
/* $Id: counts.c,v 1.9 2008/04/15 16:26:42 Bill Exp $            */
/* */
/* Author :       Eric C. Rouchka, 1996.7                                 */
/* Jun Zhu, 1996.8                                         */
/* */
/* Description :  This file contains the functions that are used in       */
/* dealing with the observed and pseudocounts.  The        */
/* following functions are included here:                  */
/* set_pseudo_counts                                    */
/* set_counts                                           */
/* adjust_counts                                        */
/* reset_counts                                         */
/* copy_counts                                          */
/* get_informative_priors                               */
/**************************************************************************/

/*------------------------------------------------------------------------*/
/* INCLUDE FILES                           */
/*------------------------------------------------------------------------*/
#include "counts.h"

/* ----------------------------------------------------------------------- */
/* Local routines                             */
/* ----------------------------------------------------------------------- */


/****************** set_pseudo_counts ************************/
/* */
/* DESCRIPTION : This function sets the pseudocounts by      */
/* setting them equal to a fraction of the     */
/* observed counts and then reading in any     */
/* informed priors if provided                 */
/* =========================================================== */


void 
set_pseudo_counts(Model B, IPtype IP, Ctype C)
{
	int             i, j, t;
	int             newlen;
	int             nNumMotifs;
	double          fraction;

	if (IP->nAlphaLen == 20)
		newlen = IP->nSeqLen - IP->nNumProcessed;
	else
		newlen = C->nTotBack;

	if (newlen == 0)
		p_error("set_pseudo_counts: invalid count. Possibly you did not use -n with nucleotide data.");

	/* Initialize Pseudocounts to zero  */
	C->dSumBGPseudo = 0.0;
	C->dTotSumMotifPseudo = 0.0;

	for (t = 0; t < IP->nNumMotifTypes; t++) {
		nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]);
		for (i = 0; i < IP->nAlphaLen; i++) {
			fraction = (double) C->fCounts[t][BG][i] / (double) newlen;
			C->dPseudoCounts[t][BG][i] =
				fraction * IP->dPseudoCntWt * (double) nNumMotifs;
			if (C->dPseudoCounts[t][BG][i] < MIN_PSEUDOCNT)
				C->dPseudoCounts[t][BG][i] = MIN_PSEUDOCNT;
			if (t == 0)
				C->dSumBGPseudo += C->dPseudoCounts[t][BG][i];
		}
	}

	/* Set the motif pseudocounts to the background pseudocounts */
	/*
	   for(t = 0; t < IP->nNumMotifTypes; t++) {
	     for(j = 1; j <=IP->nMotifLen[t]; j++) {
	       for(i = 0; i < IP->nAlphaLen; i++) {
	            C->dPseudoCounts[t][j][i] =
	                C->dPseudoCounts[t][BG][i];
		    C->dSumPseudo[t][MOTIF][i] = C->dSumBGPseudo;
	         }
	      }
	   }
	*/

	/* BT 3/05/97 */

	for (t = 0; t < IP->nNumMotifTypes; t++) {
		for (i = 0; i < IP->nAlphaLen; i++) {
			C->dSumPseudo[t][MOTIF][i] = C->dSumBGPseudo;
			for (j = 1; j <= IP->nMaxMotifLen[t]; j++) {
				C->dPseudoCounts[t][j][i] = C->dPseudoCounts[t][BG][i];
				if (t == 0)
					C->dTotSumMotifPseudo += C->dPseudoCounts[t][j][i];
			}
		}
	}

	if (IP->is_defined[cl_P]) {
		ReadPriorFile(B, C);	/* BT 7/3/98 */
	}
}


void 
get_informative_priors(Model B, IPtype IP, Ctype C)
{
	/* =========================================================== */
	/* FUNCTION NAME : get_informative_priors                    */
	/* */
	/* DESCRIPTION : This function reads in the informative      */
	/* motif alignment prior data from a file      */
	/* and sets the motif pseudocounts accordingly */
	/* =========================================================== */

	int             i, j, index;
	double          pseudocnt;
	double          dsum = 0.0;
	char            aminoacids[20] = "acdefghiklmnpqrstbjo";
	double        **bg_pseudo_count;	/* for quick access */
	int             t;
	int             nNumMotifs;
	int             ch;
	char            line[1024];
	int             n;
	double          weight[21];
	char           *ptr;
	int             nEnd;

	bg_pseudo_count = C->dPseudoCounts[0];
	C->dTotSumMotifPseudo = 0.0;

	for (t = 0; t < IP->nNumMotifTypes; t++) {	/* BT 4/25/97 */
		nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]);	/* BT 5/7/97 */
		bg_pseudo_count = C->dPseudoCounts[t];
		ch = fgetc(IP->Datafiles->prior_fpt);
		while ((ch != EOF) && (ch == '\n')) {	/* skip blank lines */
			ch = fgetc(IP->Datafiles->prior_fpt);
		}

		if (IP->is_defined[cl_d])	/* BT 12/19/97 */
			nEnd = IP->nMaxMotifLen[t];
		else
			nEnd = IP->nMotifLen[t];

		for (j = 1; j <= nEnd; j++) {
			n = 0;
			if (j > 1)
				ch = fgetc(IP->Datafiles->prior_fpt);
			if (ch == '\n')
				break;	/* BT 12/24/97 */
			while ((ch != EOF) && (ch != '\n')) {
				line[n++] = ch;
				ch = fgetc(IP->Datafiles->prior_fpt);
			}
			line[n] = '\0';
			if (n != 0) {	/* BT 5/30/97 */
				ptr = line;
				for (i = 0; i < IP->nAlphaLen; i++)
					weight[i] = strtod(ptr, &ptr);
				weight[IP->nAlphaLen] = strtod(ptr, &ptr);
				if (weight[IP->nAlphaLen] == 0.0)
					weight[IP->nAlphaLen] = IP->dPseudoCntWt;
				for (i = 0; i < IP->nAlphaLen; i++) {
					pseudocnt = weight[i] / 100.0;
					if (IP->nAlphaLen == 20)
						index = (int) aminoacids[i] - 97;
					else
						index = i;
					/*
					 * bg_pseudo_count[j][index] *=
					 * pseudocnt;
					 */
					bg_pseudo_count[j][index] = pseudocnt * weight[IP->nAlphaLen] * (double) nNumMotifs;	/* BT 5/7/97 */
					if (bg_pseudo_count[j][index] < MIN_PSEUDOCNT)
						bg_pseudo_count[j][index] = MIN_PSEUDOCNT;
					dsum += bg_pseudo_count[j][index];
				}
			}
		}
	}
	C->dTotSumMotifPseudo = dsum;
}


/********************  copy_counts **************************************/
/* */
/* DESCRIPTION : Copies the original counts -- no motifs have been      */
/* chosen at this point.                                  */
/************************************************************************/

void 
copy_counts(Model B)
{
	int             i, j, t;
	Ctype           C, First;	/* for quick access the count data */
	RPType          RP;

	C = B->C;
	First = B->First;
	RP = B->RP;

	C->nTotBack = First->nTotBack;
	C->dBackCnt = First->dBackCnt;
	C->dSumBGPseudo = First->dSumBGPseudo;
	C->dTotSumMotifPseudo = First->dTotSumMotifPseudo;

	for (t = 0; t < B->IP->nNumMotifTypes; t++) {
		for (i = 0; i < B->IP->nAlphaLen; i++) {
			C->dSumPseudo[t][MOTIF][i] = First->dSumPseudo[t][MOTIF][i];	/* BT 5/16/97 */
			for (j = 0; j <= B->IP->nMaxMotifLen[t]; j++) {
				C->fCounts[t][j][i] = 0;
				C->wCounts[t][j][i] = 0;
				C->dPseudoCounts[t][j][i] = First->dPseudoCounts[t][j][i];
			}
		}
		/* we need to remember the backround counts */
		for (i = 0; i < B->IP->nAlphaLen; i++) {
			C->fCounts[t][BG][i] = First->fCounts[t][BG][i];
			C->wCounts[t][BG][i] = First->wCounts[t][BG][i];
		}
	}

	for (t = 0; t < B->IP->nNumMotifTypes; t++) {
		for (i = 0; i < B->IP->RevComplement + 1; i++) {
			C->dmodel_sites[t][i] = First->dmodel_sites[t][i];
			C->dmodel_pseudo[t][i] = First->dmodel_pseudo[t][i];
		}
		C->dtot_sites[t] = First->dtot_sites[t];
		C->dtot_pseudo[t] = First->dtot_pseudo[t];
		C->dTot[t] = First->dTot[t];
		C->dbg_pseudo[t] = First->dbg_pseudo[t];
	}

	for (i = 0; i < B->IP->nNumSequences; i++) {
		for (j = 0; j < SequenceLength(B, i); j++) {
			for (t = 0; t < B->IP->nNumMotifTypes; t++) {
				RP->sitePos[i][j][t].nMotifStart = 0;
				RP->sitePos[i][j][t].nRevComp = 0;
			}
		}
	}

	for (i = 0; i < B->IP->nSeqLen; i++)
		B->RP->nInMotif[i] = 0;
}

/********************* set_counts ***********************************/
/* */
/* */
/* DESCRIPTION : Sets the counts of the number of times each        */
/* individual residue occurs throughout the sequence  */
/********************************************************************/

void 
set_counts(Model B)
{

	int             i, t, index;
	char           *Residual;	/* for quick access sequence */
	IPtype          IP;	/* for quick access */
	Ctype           First;	/* for quick access */
	int             nSeq;
	double          wt;

	Residual = (*B->Seq->R);
	IP = B->IP;
	First = B->First;
	B->C->nMaskCount = 0;	/* BT 4/16/97 */
	First->dBackCnt = 0;

	IP->nMotifSampleCnt = MOTIF_SAMPLE_CNT;	/* BT 5/7/97 */

	for (t = 0; t < IP->nNumMotifTypes; t++) {
		for (i = 0; i <= IP->nAlphaLen; i++) {
			First->fCounts[t][BG][i] = 0;
			First->wCounts[t][BG][i] = 0;
		}
	}

	for (i = 0; i < IP->nSeqLen; i++) {
		/* Go through the sequence                         */
		/* and set the background  counts according to the */
		/* residues found                                 */

		if (B->WT == NULL)
			wt = 1.0;
		else {
			nSeq = SequenceFromPosition(B, i);
			wt = GetSeqWeight(B, nSeq, t);
		}

		if (Residual[i] == 'x' ||
		    (Residual[i] == 'n' && IP->nAlphaLen != 20)) {
			First->nMaskCount += wt;
			index = B->IP->nAlphaLen;
		} else
			index = CH2INT(i, B);


		First->fCounts[0][BG][index] += 1.0;
		First->wCounts[0][BG][index] += wt;
		First->dBackCnt += wt;

	}

	for (t = 1; t < IP->nNumMotifTypes; t++) {
		for (i = 0; i <= IP->nAlphaLen; i++) {
			First->fCounts[t][BG][i] = First->fCounts[0][BG][i];
			First->wCounts[t][BG][i] = First->wCounts[0][BG][i];
		}
	}

	for (t = 0; t < IP->nNumMotifTypes; t++)
		SetPossibleSites(B, t);	/* BT 4/16/97 */

	B->C->nMaskCount = First->nMaskCount;	/* BT 4/16/97 */
	if (B->WT == NULL) {
		First->nTotBack = IP->nSeqLen - B->C->nMaskCount;
		First->dBackCnt = First->nTotBack;
	} else {
		for (First->dBackCnt = 0, First->nTotBack = 0, i = 0; i < IP->nAlphaLen; i++) {
			First->nTotBack += First->wCounts[BG][BG][i];
			First->dBackCnt += First->wCounts[BG][BG][i];
		}
	}
}


/************* SetPossibleSites           BT 4/16/97 ****************/
/* SetPossibleSites	recount the possible sites, making sure we  */
/* don't allow overlap of low complexity       */
/* regions.				    */
/* */
/* B	the current mode				    */
/* t	the motif number				    */
/********************************************************************/

void 
SetPossibleSites(Model B, int t)
{
	int             i;
	int             n = 0;
	char           *residual;
	int             nMotifLen;

	nMotifLen = MotifWidth(B, t);

	residual = (*B->Seq->R);
	B->IP->nPossSites[t] = 0;

	for (i = 0; i < B->IP->nNumSequences; i++) {
		while (n < (*B->Seq->nvEndLocs)[i]) {	/* BT 6/20/97 */
			if (!((B->IP->nAlphaLen == 4 && residual[n] == 'n') ||
			      (residual[n] == 'x'))) {
				if (n < (*B->Seq->nvEndLocs)[i] - nMotifLen + 1) {	/* BT 6/20/97 */
					if (!((B->IP->nAlphaLen == 4 &&
					       residual[n + nMotifLen - 1] == 'n') ||
					      (residual[n + nMotifLen - 1] == 'x')))
						B->IP->nPossSites[t]++;
					else
						n += nMotifLen - 1;	/* BT 6/20/97 */
				}
			}
			n++;
		}
	}
}


/*************   reset_counts  *************************/
/* */
/* DESCRIPTION : This function goes through and        */
/* for each of the initial motif         */
/* starting positions adds the counts    */
/* of the residues to the motif model    */
/*******************************************************/

void 
reset_counts(Model B, int **startPos, PoSition ** Pos)
{
	int             i;
	int             t;
	int             startloc;
	int            *startPos_t;	/* startPos_t=startPos[t] */
	PoSition       *Pos_t;
	IPtype          IP;

	IP = B->IP;

	for (t = 0; t < IP->nNumMotifTypes; t++) {
		startPos_t = startPos[t];
		Pos_t = Pos[t];
		for (i = 0; i < IP->nNumMotifs[t][FORWARD] + IP->nNumMotifs[t][REVERSE]; i++) {
			startloc = startPos_t[i];
			adjust_counts(B, ADD, startloc, t, Pos_t[startloc].RevComp);
		}
	}
}


/*********************  adjust_counts  ***********************/
/* */
/* DESCRIPTION : This function updates the background and    */
/* motif residue counts by either adding or    */
/* deleting a motif (based on the flag)        */
/* starting at location begin_loc.             */
/* =========================================================== */

void 
adjust_counts(Model B, int flag, int begin_loc,
	      int t, short RevComp)
{
	int             i, j, index, mid;
	double          add_val;/* BT 4/18/2000 */
	int             n, J, compindex, orig_index;
	int             last = -1, curr, first = -1;
	int             curr_fwd, curr_back, last_fwd = -1, last_back;
	char            processed;
	char           *Res_ptr = B->Seq->R[0];	/* for quick access sequence
						 * residue */
	int             compindex2;
	int             nSeq;
	int             nPos;
	int             frag_width = 0;
	double          wt;
	double          addCnt;

	nSeq = SequenceFromPosition(B, begin_loc);
	nPos = begin_loc - SequenceStartPos(B, nSeq);
	B->RP->sitePos[nSeq][nPos][t].nMotifStart = flag;
	B->RP->sitePos[nSeq][nPos][t].nRevComp = RevComp;

	if (!B->IP->is_defined[cl_F]) {
		first = FirstCol(B->F, t);
		last = last_fwd = -1;
		last_back = B->F->nMaxLen[t];
		frag_width = LastCol(B->F, t) - first + 1;
	}
	if (flag == ADD) {
		add_val = 1;
		if (B->WT != NULL)
			wt = GetSeqWeight(B, nSeq, t);
		else
			wt = 1;

		B->RP->sitePos[nSeq][nPos][t].addVal = wt;
	} else {
		add_val = -1;
		if (B->WT != NULL)
			wt = -B->RP->sitePos[nSeq][nPos][t].addVal;
		else
			wt = -1;
		B->RP->sitePos[nSeq][nPos][t].addVal = 0;
	}
	addCnt = wt;

	/*
	 * J = B->IP->nMotifLen[t];  mid = J / 2 + J % 2 - 1;
	 *//* BT 4/18/97 */

	J = 0;
	for (i = 0; i < B->IP->nNumSequences; i++) {	/* BT 3/24/97 - allow
							 * overlap */
		if (begin_loc < (*B->Seq->nvEndLocs)[i]) {
			J = min(B->IP->nMotifLen[t], (*B->Seq->nvEndLocs)[i] - begin_loc);
			break;
		}
	}

	if (J == 0)
		return;		/* BT 4/18/97 */

	mid = J / 2 + J % 2 - 1;

	/*--------------------------*/
	/* DEALING WITH AMINO ACIDS */
	/*--------------------------*/

	if (B->IP->nAlphaLen == 20) {
		for (j = 0; j < J; j++) {
			if (!B->IP->is_defined[cl_F]) {	/* Fragmentation is used */
				curr = NextCol(B->F, t, last);
				if (curr == -1)
					p_internal_error("No next column");
				last = curr;
				curr -= first;
				if (curr < 0)
					p_internal_error("Invalid location");
			} else
				curr = j;

			B->RP->nInMotif[begin_loc + curr] = flag;

			processed = Res_ptr[begin_loc + curr];
			if ((processed != 'x') && (processed != 'X') && (processed != 'u')) {
				index = (int) (processed) - 97;
				B->C->fCounts[t][j + 1][index] += add_val;
				B->C->wCounts[t][j + 1][index] += wt;
				for (i = 0; i < B->IP->nNumMotifTypes; i++) {
					B->C->fCounts[i][BG][index] -= add_val;
					B->C->wCounts[i][BG][index] -= wt;
				}
			}
		}
	}
	/*--------------------------*/
	/* DEALING WITH NUCLEOTIDES */
	/*--------------------------*/

	else {
		for (j = 0; j < J; j++) {
			if (!B->IP->is_defined[cl_F]) {	/* Fragmentation */
				curr_fwd = NextCol(B->F, t, last_fwd);
				/*
				 * curr_back = PrevCol(B->F, t, last_back); *//* 9
				 * /8/99
				 */
				last_fwd = curr_fwd;
				/* last_back = curr_back; *//* 9/8/99 */
				curr_fwd -= first;
				/* curr_back -= first; *//* 9/8/99 */

				curr_back = frag_width - curr_fwd - 1;	/* BT 08/11/99 */

				if ((curr_fwd < 0) || (curr_fwd > B->F->nMaxLen[t])) {
					printf("first %d curr_fwd = %d\n", first, curr_fwd);
					p_internal_error("Invalid forward column");
				}
				if ((curr_back < 0) || (curr_back > B->F->nMaxLen[t]))
					p_internal_error("Invalid backward column");
			} else {/* No Fragmentation */
				curr_fwd = j;
				curr_back = J - j - 1;
			}
			if (!RevComp) {
				curr = curr_fwd;
			} else {
				curr = curr_back;
			}


			processed = Res_ptr[begin_loc + curr_fwd];
			orig_index = (int) processed - 97;
			if ((processed == 'X') || (processed == 'n') ||
			    (processed == 'x') || (processed == 'N')) {	/* BT 4/9/97 */
				orig_index = 4;
			}
			if (!RevComp) {	/* FORWARD MOTIF */
				n = begin_loc + curr_fwd;
				index = orig_index;
				compindex = CH2INTCOMP(n, B);
				compindex2 = CH2INTCOMP(begin_loc + curr_back, B);
			} else {/* REVERSE COMPLEMENT */
				n = begin_loc + curr_back;
				index = CH2INTCOMP(n, B);
				compindex = orig_index;
				compindex2 = CH2INTCOMP(begin_loc + curr_fwd, B);
			}

			B->RP->nInMotif[begin_loc + curr] = flag;

			for (i = 0; i < B->IP->nNumMotifTypes; i++) {
				B->C->fCounts[i][BG][orig_index] -= add_val;
				B->C->wCounts[i][BG][orig_index] -= wt;
			}
			B->C->nTotBack -= wt;
			B->C->dBackCnt -= addCnt;
			B->C->fCounts[t][j + 1][index] += add_val;
			B->C->wCounts[t][j + 1][index] += wt;
		}
	}

#ifdef _DEBUG_
	for (j = 0; j < J; j++) {
		for (i = 0; i < B->IP->nAlphaLen; i++) {
			if (B->C->fCounts[t][j + 1][i] < -EPS)
				p_internal_error("Negative frequency count.");
		}
	}
#endif
}
