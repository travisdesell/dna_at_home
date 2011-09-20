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
/* $Id: column.c,v 1.8 2009/04/23 18:43:53 Bill Exp $            */
/* */
/* Author :       Eric C. Rouchka May 29, 1996                            */
/* Jun Zhu, October 20, 1996                               */
/* */
/* Description :  Contains the function ColumnShift which takes the       */
/* alignment, shifts it left and right, and samples one of */
/* the alignments based on the map probability of each of  */
/* the alignments                                          */
/**************************************************************************/


#include "column.h"

/********************   ColumnShift   *****************************/
/* */
/* */
/* DESCRIPTION : This function takes the current motif alignment, */
/* shifts it left and right, and samples one of     */
/* those alignments based on the map probability of */
/* each of the alignments                           */
/******************************************************************/

void 
ColumnShift(Model B, Mlist M, PoSition ** Pos, ProbStruct * P)
{
	int             t;
	int             i;
	int             maxleft;
	int             maxright;
	MlistEl        *curr;
	double         *ProbArray;
	double          dProb;
	int             shift;
	double          dMotifProb;
	double          dBGProb;
	double          totalProb;
	double          r;
	int             index;
	int             pos;
	int            *motifSave;

	if (M == NULL)
		return;

	for (t = 0; t < B->IP->nNumMotifTypes; t++) {
		maxleft = maxright = B->IP->nMotifLen[t] / 2;

		curr = M[t]->Motifs;
		while (curr != NULL) {	/* Determine the max */
			not_in_motif(Pos, curr->pos, B, t);

			for (i = 1; i <= maxleft; i++) {
				if (!curr->RevComp)
					pos = curr->pos - i;
				else
					pos = curr->pos + i;
				if (!PossStartPos(Pos, pos, t, B)) {
					maxleft = i - 1;
					break;
				}
			}

			for (i = 1; i <= maxright; i++) {
				if (!curr->RevComp)
					pos = curr->pos + i;
				else
					pos = curr->pos - i;
				if (!PossStartPos(Pos, pos, t, B)) {
					maxright = i - 1;
					break;
				}
			}

			set_in_motif(Pos, curr->pos, B, t, curr->RevComp);
			curr = curr->next;
		}

		if (maxleft > 0 || maxright > 0) {
			curr = M[t]->Motifs;
			while (curr != NULL) {
				not_in_motif(Pos, curr->pos, B, t);
				curr = curr->next;
			}

			NEW(ProbArray, maxleft + maxright + 1, double);
			NEW(motifSave, NUMMOTIFS(B->IP->nNumMotifs[t]), int);

			curr = M[t]->Motifs;
			i = 0;
			while (curr != NULL) {
				motifSave[i] = curr->pos;
				i++;
				curr = curr->next;
			}

			totalProb = 0;
			for (index = 0, shift = -maxleft; shift <= maxright; shift++) {
				curr = M[t]->Motifs;
				while (curr != NULL) {
					adjust_counts(B, DELETE, curr->pos, t, curr->RevComp);
					if (!curr->RevComp) {
						if (shift != -maxleft)
							curr->pos -= (shift - 1);
						curr->pos += shift;
					} else {
						if (shift != -maxleft)
							curr->pos += (shift - 1);
						curr->pos -= shift;
					}
					adjust_counts(B, ADD, curr->pos, t, curr->RevComp);
#ifdef _DEBUG_
					CheckCounts(B);
#endif
					curr = curr->next;
				}

				if (!B->IP->is_defined[cl_sample_model])
					update_prob(P, B, (!B->IP->is_defined[cl_b]));

				dProb = 1.0;
				curr = M[t]->Motifs;
				while (curr != NULL) {
					if (!IsPhyloSeq(B, curr->seq_num) ||
					    curr->seq_num == B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[curr->seq_num]])
						dProb *= in_motif_prob(B, *P, curr->pos, t, curr->RevComp,
						TRUE, &dMotifProb, &dBGProb);	/* BT 2/9/98 */

					curr = curr->next;
				}
				ProbArray[index] = dProb;
				totalProb += dProb;
				index++;
			}

			curr = M[t]->Motifs;
			i = 0;
			while (curr != NULL) {
				adjust_counts(B, DELETE, curr->pos, t, curr->RevComp);
				curr->pos = motifSave[i];
				adjust_counts(B, ADD, curr->pos, t, curr->RevComp);
				i++;
				curr = curr->next;
			}

			for (i = 0; i < index; i++)
				ProbArray[i] /= totalProb;

			r = drand();
			dProb = 0;
			shift = 0;
			for (i = 0; i < index; i++) {
				dProb += ProbArray[i];
				if (dProb > r) {
					shift = -maxleft + i;
					break;
				}
			}

			B->IP->col_shift[t] = shift;

			curr = M[t]->Motifs;
			while (curr != NULL) {
				adjust_counts(B, DELETE, curr->pos, t, curr->RevComp);
				if (!curr->RevComp)
					curr->pos += shift;
				else
					curr->pos -= shift;
				adjust_counts(B, ADD, curr->pos, t, curr->RevComp);
				set_in_motif(Pos, curr->pos, B, t, curr->RevComp);
#ifdef _DEBUG_
				CheckCounts(B);
#endif
				curr = curr->next;
			}

			free(ProbArray);
			free(motifSave);
			if (!B->IP->is_defined[cl_sample_model])
				update_prob(P, B, (!B->IP->is_defined[cl_b]));
		}
	}
}
