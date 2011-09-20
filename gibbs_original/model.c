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
#include "model.h"

void            MetropolisHastingsMotifModel(Model B, ProbStruct * P, ProbStruct * P1, int t);
void            MetropolisHastingsBkgndModel(Model B, ProbStruct * P, ProbStruct * P1);


/* sample probability models from the current counts and pseudocounts */

void 
SampleModel(Model B, ProbStruct * P)
{
	Ctype           C;
	IPtype          IP;
	int             t;
	int             n;
	double         *counts;
	ProbStruct      P1;

	C = B->C;
	IP = B->IP;

	if (!IP->is_defined[cl_E])
		p_error("SampleModel: -E must be users with -sample_model.");

	init_prob(&P1, IP);
	NEW(counts, IP->nAlphaLen, double);

	if (!IP->is_defined[cl_B] && !IP->is_defined[cl_freq_background]) {
		for (n = 0; n < IP->nAlphaLen; n++)
			counts[n] = C->wCounts[0][BG][n] + C->dPseudoCounts[0][BG][n];
		rdirichlet(B, IP->nAlphaLen, counts, P1.dvInBGProb[0]);
		MetropolisHastingsBkgndModel(B, P, &P1);

		for (t = 1; t < IP->nNumMotifTypes; t++) {
			for (n = 1; n < IP->nAlphaLen; n++)
				P->dvInBGProb[t][n] = P->dvInBGProb[0][n];
		}
	}
	for (t = 0; t < IP->nNumMotifTypes; t++) {
		SampleAModel(B, t, P1.dvInMotifProb[t], C->wCounts[t]);
		MetropolisHastingsMotifModel(B, P, &P1, t);
	}

	free(counts);
	free_prob(&P1, IP);
}


void 
MetropolisHastingsMotifModel(Model B, ProbStruct * P, ProbStruct * P1, int t)
{
	RPType          RP;
	Ctype           C;
	IPtype          IP;
	int             n;
	int             seq;
	int             nLen;
	int             nPos;
	int             j;
	double          motifProb;
	double          bkgProb;
	double          fels0 = 0;
	double          fels1 = 0;
	double          pTheta0 = 0;
	double          pTheta1 = 0;
	double          ratio;
	double          r;

	/* FILE       *fpt = NULL; */
	/* char       *dumpFile = NULL; */

	C = B->C;
	IP = B->IP;
	RP = B->RP;

	for (seq = 0; seq < IP->nNumSequences; seq += SeqInc(B, seq)) {
		nLen = SequenceLength(B, seq);
		for (n = 0; n < nLen; n++) {
			if (RP->sitePos[seq][n][t].nMotifStart) {
				nPos = n + SequenceStartPos(B, seq);

				in_motif_prob(B, *P, nPos, t, RP->sitePos[seq][n][t].nRevComp, FALSE,
					      &motifProb, &bkgProb);
				fels0 += log(motifProb);

				in_motif_prob(B, *P1, nPos, t, RP->sitePos[seq][n][t].nRevComp, FALSE,
					      &motifProb, &bkgProb);
				fels1 += log(motifProb);
			}
		}
	}

	for (j = 0; j < IP->nMotifLen[t]; j++) {
		for (n = 0; n < IP->nAlphaLen; n++) {
			pTheta0 += C->wCounts[t][j + 1][n] * log(P->dvInMotifProb[t][j][n]);
			pTheta1 += C->wCounts[t][j + 1][n] * log(P1->dvInMotifProb[t][j][n]);
		}
	}

	r = drand();
	ratio = (fels1 - pTheta1) - (fels0 - pTheta0);

	/* TEST - 09/05/07 */

	if (IP->bayesSampleIter >= 0) {
		RP->bayesSampleCount++;
		/* if( IP->Datafiles->output_filename ) */
		/* { */
		/* NEW( dumpFile, FILENAME_MAX, char ); */
		/*
		 * sprintf( dumpFile, "%s.%d.model",
		 * IP->Datafiles->output_filename, IP->nRank );
		 */
		/* fpt = fopen( dumpFile, "a" ); */
		/* fprintf( fpt, "%d %d %g %g %g %g %g  ",  */
		/*
		 * IP->nSeedRun, IP->bayesSampleIter, ratio, fels1, pTheta1,
		 * fels0, pTheta0 );
		 */
		/* } */
	}
	if (ratio > 0 || r <= exp(ratio)) {
		if (IP->bayesSampleIter >= 0) {
			/* if( IP->Datafiles->output_filename ) */
			/* fprintf( fpt, "accepted\n" );       */
			RP->acceptCount++;
		}
		for (j = 0; j < IP->nMotifLen[t]; j++) {
			for (n = 0; n < IP->nAlphaLen; n++) {
				P->dvInMotifProb[t][j][n] = P1->dvInMotifProb[t][j][n];
			}
		}
	}
	/* else */
	/* { */
	/* if( IP->bayesSampleIter >= 0 && IP->Datafiles->output_filename ) */
	/* fprintf( fpt, "not_accepted\n" ); */
	/* } */

	/* if( IP->bayesSampleIter >= 0 && IP->Datafiles->output_filename ) */
	/* { */
	/* fclose( fpt );	    */
	/* free( dumpFile ); */
	/* } */
}


void 
MetropolisHastingsBkgndModel(Model B, ProbStruct * P, ProbStruct * P1)
{
	RPType          RP;
	Ctype           C;
	IPtype          IP;
	PhyloType       PH;
	int             n;
	int             seq;
	int             nLen;
	int             nPos;
	int             t;
	double          fels0 = 0;
	double          fels1 = 0;
	double          pTheta0 = 0;
	double          pTheta1 = 0;
	double          ratio;
	double          r;
	int             inBkgnd;
	PhyloTree       tree;
	char           *R;

	R = (*(B)->Seq->R);
	C = B->C;
	IP = B->IP;
	RP = B->RP;
	PH = B->Phylo;

	for (seq = 0; seq < IP->nNumSequences; seq += SpeciesInc(B, seq)) {
		nLen = SequenceLength(B, seq);

		n = 0;
		while (n < nLen) {
			inBkgnd = TRUE;
			for (t = 0; t < IP->nNumMotifTypes; t++) {
				if (RP->sitePos[seq][n][t].nMotifStart) {
					inBkgnd = FALSE;
					n += MotifWidth(B, t);
					break;
				}
			}

			if (inBkgnd) {
				nPos = n + SequenceStartPos(B, seq);
				if (IsPhyloSeq(B, seq) && PH->treeCount > 0) {
					if (PossibleAlignedBkgndPos(B, nPos, nLen, seq)) {
						tree = PH->phyloTree[PH->phyloIndex[seq]];
						CalcSubsMatrix(B, tree, P->dvInBGProb[0]);
						fels0 += log(CalcFelsProb(B, tree, P->dvInBGProb[0], nPos, FALSE));

						CalcSubsMatrix(B, tree, P1->dvInBGProb[0]);
						fels1 += log(CalcFelsProb(B, tree, P1->dvInBGProb[0], nPos, FALSE));
					}
				} else {
					fels0 += log(P->dvInBGProb[0][RCH2INT(nPos, R)]);
					fels1 += log(P1->dvInBGProb[0][RCH2INT(nPos, R)]);
				}
				n++;
			}
		}
	}

	for (n = 0; n < IP->nAlphaLen; n++) {
		pTheta0 += C->wCounts[0][BG][n] * log(P->dvInBGProb[0][n]);
		pTheta1 += C->wCounts[0][BG][n] * log(P1->dvInBGProb[0][n]);
	}

	r = drand();
	ratio = (fels1 - pTheta1) - (fels0 - pTheta0);
	if (ratio > 0 || r <= exp(ratio)) {
		for (n = 0; n < IP->nAlphaLen; n++) {
			P->dvInBGProb[0][n] = P1->dvInBGProb[0][n];
		}
	}
}
