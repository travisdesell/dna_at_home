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


#include "centroid.h"


#define FLANK       5
#define CENT_START  1
#define CENT_POS    2
#define FLANK_POS   3
#define LINE_SIZE   FILENAME_MAX+2048
#define TOKEN_SIZE  1024

typedef struct {
	char           *header;
	char           *seq;
	int             pos;
	int             isPhyloSeq;
}               SeqStruct;

typedef struct {
	int             seq;
	int             pos;
	char           *line;
}               LineStruct;


int             CalcCentroid(Model B, int seq, int *centroidPos, FILE * fastaFpt, FILE * wilcoxFpt);
int 
PrintCentroid(Model B, int seq, int *results, int *sum, int *centroidPos,
	      FILE * fastaFpt, FILE * wilcoxFptr);
void            AlignCentroidModel(Model B, int *centroidPos, int totalSites, char *fastaFile);
void            WritePriorFile(Model B, char *priorFile, char *fastaFile, char *origPriorFile);
void            SortFastaFile(Model B, char *fastaFile, FILE * priorFpt);
void            PrintTree(Model B, int seq, FILE * priorFpt);
void            PrintTree2(PhyloTree tree, FILE * priorFpt);
void            GetProgramName(Model B, char *name);
void            ConstructCommandLine(Model B, char *fastaFile, char *command, char *origPriorFile);
void            ParseOutputFile(Model B, char *outFile, FILE * fpt);
void            ParseLine(Model B, char *line, int *prevSeq, int *seqCnt, LineStruct * lineStruct);
void            GetLeftFlank(Model B, int seq, int pos, int rev, int siteWidth, char *flank);
void            GetRightFlank(Model B, int seq, int pos, int rev, int siteWidth, char *flank);
void            WriteCentroidWilcoxon(Model B, char *wilcoxFile);
void            Credibility(Model B, int *centroidPos);
int             SampleDist(SampleStruct * centroid, int centCnt, SampleStruct * sample, int sampleCnt);
int             OverlapDist(int start1, int width1, int start2, int width2);
int             compare_ints(const void *a, const void *b);
int             CompLineStruct(const void *p1, const void *p2);
void            DumpSamples(Model B);
void            DumpDistance(Model B, int seq, int *dist);


void 
Centroid(Model B)
{
	IPtype          IP;
	RPType          RP;
	FILE           *fpt;
	int             seq;
	int             totalSites = 0;
	int            *centroidPos;
	char           *fastaFile = NULL;
	FILE           *fastaFpt = NULL;
	char           *wilcoxFile = NULL;
	FILE           *wilcoxFpt = NULL;

	IP = B->IP;
	RP = B->RP;

	NEW(centroidPos, IP->nSeqLen, int);

	if (IP->is_defined[cl_align_centroid]) {
		fastaFile = GibbsTempnam(NULL, "cen");
		fastaFpt = fopen(fastaFile, "w");
		if (fastaFpt == NULL)
			p_error("Centroid: Unable to open FASTA output file.");

		chmod(fastaFile, S_IRUSR | S_IWUSR);
	}
	if (IP->is_defined[cl_wilcox]) {
		wilcoxFile = GibbsTempnam(NULL, "wil");
		wilcoxFpt = fopen(wilcoxFile, "w");
		if (wilcoxFpt == NULL)
			p_error("Centroid: Unable to open Wilcoxon output file.");

		chmod(wilcoxFile, S_IRUSR | S_IWUSR);
	}
	for (seq = 0; seq < IP->nNumSequences; seq++)
		totalSites += CalcCentroid(B, seq, centroidPos, fastaFpt, wilcoxFpt);

	fpt = B->IP->Datafiles->out_fpt;
	fprintf(fpt, "Num Sites: %d\n", totalSites);

	if (totalSites) {
		fprintf(fpt, "\nColumn 1 :  Sequence Number, Site Number\n");
		fprintf(fpt, "Column 2 :  Left End Location\n");
		fprintf(fpt, "Column 4 :  Motif Element\n");
		fprintf(fpt, "Column 6 :  Right End Location\n");
		fprintf(fpt, "Column 7 :  Probability of Element\n");
		fprintf(fpt, "Column 8 :  Sequence Description from FastA input\n\n");
	}
	if (IP->is_defined[cl_wilcox]) {
		fclose(wilcoxFpt);
		if (totalSites > 0)
			WriteCentroidWilcoxon(B, wilcoxFile);
		remove(wilcoxFile);
		free(wilcoxFile);
	}
	if (!IP->is_defined[cl_no_cred])
		Credibility(B, centroidPos);

	if (IP->is_defined[cl_sample_model] && B->Phylo->phyloTree) {
		fprintf(fpt, "M-H proposal accepted %d times out of %d sampled (%5.2f%%)\n\n",
			B->RP->acceptCount, RP->bayesSampleCount,
			(((double) B->RP->acceptCount) / (RP->bayesSampleCount)) * 100);
	}
	if (IP->is_defined[cl_align_centroid]) {
		fclose(fastaFpt);
		if (totalSites > 0)
			AlignCentroidModel(B, centroidPos, totalSites, fastaFile);
		remove(fastaFile);
		free(fastaFile);
	}
	free(centroidPos);
}


int 
CalcCentroid(Model B, int seq, int *centroidPos, FILE * fastaFpt, FILE * wilcoxFpt)
{
	RPType          RP;
	IPtype          IP;
	int           **cumStarts;
	int             nLen;
	int             totalSites;
	int             nWidths;
	int             pos;
	int            *best;
	int            *perfect;
	int            *site;
	int            *results;
	int            *sum;
	int             minWidth;
	int             w;
	int             w2;
	int             iWidth;
	int             jWidth;
	int             sufficientOverlaps;
	int             requiredOverlap;
	int             perfectOverlaps;
	int             alternate;
	int             samples;
	int             siteCount;

	IP = B->IP;
	RP = B->RP;

	nWidths = NumWidths(B);
	nLen = SequenceLength(B, seq);
	minWidth = MinPossibleWidth(B);
	samples = IP->bayesSamplePeriod * IP->nSeeds;

	NEWP(cumStarts, nLen + 1, int);
	for (pos = 0; pos <= nLen; pos++)
		NEW(cumStarts[pos], nWidths, int);

	NEW(best, nLen + 1, int);
	NEW(perfect, nLen + 1, int);
	NEW(site, nLen + 1, int);
	NEW(results, nLen + 1, int);
	NEW(sum, nLen + 1, int);

	totalSites = 0;
	for (w = 0; w < nWidths; w++) {
		for (pos = 1; pos <= nLen; pos++)
			cumStarts[pos][w] = cumStarts[pos - 1][w] + RP->totalPosCounts[seq][pos - 1][w];
		totalSites += cumStarts[nLen][w];
	}

	best[0] = totalSites;
	perfect[0] = 0;
	site[0] = 0;
	for (pos = 1; pos <= nLen; pos++) {
		best[pos] = best[pos - 1];
		perfect[pos] = perfect[pos - 1];
		site[pos] = 0;

		for (w = 0; w < nWidths; w++) {
			iWidth = minWidth + w;
			if (iWidth > pos)
				break;

			sufficientOverlaps = 0;
			for (w2 = 0; w2 < nWidths; w2++) {
				jWidth = minWidth + w2;
				requiredOverlap = 1 + (iWidth > jWidth ? iWidth : jWidth) / 2;
				if (requiredOverlap > iWidth || requiredOverlap > jWidth)
					continue;

				sufficientOverlaps += cumStarts[pos - requiredOverlap + 1][w2]
					- (pos + requiredOverlap >= iWidth + jWidth ?
					   cumStarts[pos - iWidth - jWidth + requiredOverlap][w2] : 0);
			}

			alternate = best[pos - iWidth] + samples - 2 * sufficientOverlaps;
			perfectOverlaps = perfect[pos - iWidth] + RP->totalPosCounts[seq][pos - iWidth][w];

			if (best[pos] > alternate ||
			    (best[pos] == alternate && perfect[pos] < perfectOverlaps)) {
				best[pos] = alternate;
				perfect[pos] = perfectOverlaps;
				site[pos] = iWidth;
				sum[pos] = sufficientOverlaps;
			}
		}
	}

	for (pos = nLen; pos > 0; /* nothing */ ) {
		if (site[pos] != 0) {
			jWidth = site[pos];
			for (iWidth = site[pos]; iWidth > 0; --iWidth)
				results[--pos] = jWidth;
		} else
			results[--pos] = 0;
	}

	siteCount = PrintCentroid(B, seq, results, sum, centroidPos, fastaFpt, wilcoxFpt);

	FREEP(cumStarts, nLen + 1);
	free(best);
	free(perfect);
	free(site);
	free(results);
	free(sum);

	return (siteCount);
}


int 
PrintCentroid(Model B, int seq, int *results, int *sum, int *centroidPos,
	      FILE * fastaFpt, FILE * wilcoxFpt)
{
	IPtype          IP;
	FILE           *fpt;
	int             nLen;
	int             pos;
	int             i;
	int             j;
	char           *letters;
	char           *desc;
	double          prob;
	int             start;
	int             maxWidth;
	int             width;
	int             siteCount = 0;

	IP = B->IP;

	maxWidth = MaxPossibleWidth(B);

	NEW(desc, max(maxWidth + 1, 51), char);
	strncpy(desc, B->IP->fastaHeader[seq], 50);

	nLen = SequenceLength(B, seq);
	start = SequenceStartPos(B, seq);
	fpt = B->IP->Datafiles->out_fpt;
	letters = &((*B->Seq->R)[start]);

	for (pos = 0; pos < nLen; /* nothing */ ) {
		if (results[pos]) {
			width = results[pos];
			fprintf(fpt, "%5d, %d %7d ", seq + 1, siteCount + 1, pos + 1);
			for (j = 0, i = max(0, pos - FLANK); i < pos; i++, j++) {
				desc[j] = curr_ch(letters[i], B->IP->nAlphaLen, FALSE);
				centroidPos[start + i] = FLANK_POS;
			}
			desc[j] = '\0';
			fprintf(fpt, "%5s", desc);
			fprintf(fpt, " ");

			if (IP->is_defined[cl_align_centroid]) {
				fprintf(fastaFpt, ">%d %d\n", seq, max(0, pos - FLANK));
				fprintf(fastaFpt, "%s", desc);
			}
			memset(desc, ' ', maxWidth);
			desc[maxWidth] = '\0';
			for (j = 0; j < width; j++, pos++) {
				desc[j] = curr_ch(letters[pos], B->IP->nAlphaLen, TRUE);
				if (j == 0)
					centroidPos[start + pos] = CENT_START;
				else
					centroidPos[start + pos] = CENT_POS;
			}
			fprintf(fpt, "%s", desc);
			fprintf(fpt, " ");

			if (IP->is_defined[cl_align_centroid]) {
				for (j = 0; j < width; j++)
					fprintf(fastaFpt, "%c", desc[j]);
			}
			for (j = 0, i = pos; i < min(pos + FLANK, nLen); i++, j++) {
				desc[j] = curr_ch(letters[i], B->IP->nAlphaLen, FALSE);
				centroidPos[start + i] = FLANK_POS;
			}
			desc[j] = '\0';
			fprintf(fpt, "%5s", desc);

			if (IP->is_defined[cl_align_centroid])
				fprintf(fastaFpt, "%s\n", desc);

			fprintf(fpt, " %7d", pos);
			prob = ((double) sum[pos]) / (IP->bayesSamplePeriod * IP->nSeeds);
			fprintf(fpt, " %5.2f ", prob);

			if (IP->is_defined[cl_wilcox]) {
				if (strstr(IP->fastaHeader[seq], "Wilcoxon") != IP->fastaHeader[seq])
					fprintf(wilcoxFpt, "0.0 %.4f\n", prob);
				else
					fprintf(wilcoxFpt, "%.4f 0.0\n", prob);
			}
			strncpy(desc, IP->fastaHeader[seq], 50);
			fprintf(fpt, "%s\n", desc);
			siteCount++;
		} else
			pos++;
	}

	free(desc);

	return (siteCount);
}


void 
AlignCentroidModel(Model B, int *centroidPos, int totalSites, char *fastaFile)
{
	IPtype          IP;
	RPType          RP;
	FILE           *fpt;
	char           *outFile;
	char           *tempFile;
	char           *priorFile;
	char           *command;
	struct tm       tm_buf;
	struct tm      *localTime;
	time_t          tloc;
	char           *origPriorFile;
	int             ret;

	IP = B->IP;
	RP = B->RP;

	fpt = B->IP->Datafiles->out_fpt;

	fprintf(fpt, "\n\n\n\n\n======================================================================\n");
	fprintf(fpt,
		"======================== Aligned Centroid Sites ======================\n");
	fprintf(fpt,
		"======================================================================\n\n");

	NEW(command, LINE_SIZE, char);
	outFile = GibbsTempnam(NULL, "cen");
	tempFile = GibbsTempnam(NULL, "cen");
	priorFile = GibbsTempnam(NULL, "cen");
	NEW(origPriorFile, LINE_SIZE, char);

	ConstructCommandLine(B, fastaFile, command, origPriorFile);

	WritePriorFile(B, priorFile, fastaFile, origPriorFile);

	strcat(command, "-o ");
	strcat(command, outFile);
	strcat(command, " -P ");
	strcat(command, priorFile);
	strcat(command, " >");
	strcat(command, tempFile);

	if (!IP->is_defined[cl_Z]) {
		time(&tloc);
		localTime = localtime_r(&tloc, &tm_buf);
		fprintf(stdout, "Launching %s at %d:%d:%d\n", command,
		  localTime->tm_hour, localTime->tm_min, localTime->tm_sec);
	}
	ret = system(command);
	if (!IP->is_defined[cl_Z]) {
		time(&tloc);
		localTime = localtime_r(&tloc, &tm_buf);
		fprintf(stdout, "Finished %s at %d %d %d return = %d\n", command,
			localTime->tm_hour, localTime->tm_min, localTime->tm_sec, ret);
	}
	ParseOutputFile(B, outFile, fpt);
	print_col_desc(fpt);

	free(command);
	remove(outFile);
	remove(tempFile);
	remove(priorFile);
	free(tempFile);
	free(outFile);
	free(priorFile);
	free(origPriorFile);
}


void 
WritePriorFile(Model B, char *priorFile, char *fastaFile, char *origPriorFile)
{
	IPtype          IP;
	RPType          RP;
	FILE           *priorFpt;
	FILE           *origFpt;
	char           *line;

	IP = B->IP;
	RP = B->RP;

	priorFpt = fopen(priorFile, "w");

	fprintf(priorFpt, ">BLOCKS 1 1\n");
	fprintf(priorFpt, "0 1\n");
	fprintf(priorFpt, ">\n\n");

	if (strlen(origPriorFile) > 0) {
		NEW(line, LINE_SIZE + 1, char);
		origFpt = fopen(origPriorFile, "r");
		while (GetALine(origFpt, line, LINE_SIZE, FALSE)) {
			if (strstr(line, ">PSEUDO") != NULL) {
				fprintf(priorFpt, "%s\n", line);
				while (GetALine(origFpt, line, LINE_SIZE, FALSE)) {
					fprintf(priorFpt, "%s\n", line);
					if (strstr(line, ">") != NULL) {
						fprintf(priorFpt, "\n");
						break;
					}
				}
			}
		}
		fclose(origFpt);
		free(line);
	}
	if (B->Phylo->treeCount > 0 || B->WT)
		SortFastaFile(B, fastaFile, priorFpt);

	fclose(priorFpt);
}


void 
SortFastaFile(Model B, char *fastaFile, FILE * priorFpt)
{
	IPtype          IP;
	RPType          RP;
	FILE           *fastaFpt;
	char           *line;
	char           *token;
	char           *ptr;
	SeqStruct     **seqs;
	int             i;
	int             j;
	int             siteCnt = 0;
	int             seqCnt = 0;
	int             seqCnt2 = 0;
	int             tokenCnt;
	int             sq = 0;
	int             prevSq = -1;
	int            *treeIndex;
	int             treeStart;
	int             treeEnd;
	int             origTree;
	int             pos;
	int             origPos;
	int            *posIndex;
	int            *phyloSpeciesSample;

	IP = B->IP;
	RP = B->RP;

	NEWP(seqs, IP->nNumSequences, SeqStruct);
	for (i = 0; i < IP->nNumSequences; i++) {
		NEW(seqs[i], RP->nMaxBlocks, SeqStruct);
	}

	NEW(line, LINE_SIZE + 1, char);
	NEW(token, LINE_SIZE + 1, char);

	fastaFpt = fopen(fastaFile, "r");
	while (GetALine(fastaFpt, line, LINE_SIZE, FALSE)) {
		if (strstr(line, ">") != NULL) {
			tokenCnt = 0;
			ptr = line;
			ptr = GetAToken(ptr, token);
			sq = atoi(&token[1]);
			if (sq != prevSq)
				siteCnt = 0;
			prevSq = sq;
			NEW(seqs[sq][siteCnt].header, strlen(line) + 1, char);
			strcpy(seqs[sq][siteCnt].header, line);
			seqs[sq][siteCnt].isPhyloSeq = IsPhyloSeq(B, sq);
			ptr = GetAToken(ptr, token);
			pos = atoi(token);
			seqs[sq][siteCnt].pos = pos;
		} else if (strlen(line) != 0) {
			NEW(seqs[sq][siteCnt].seq, strlen(line) + 1, char);
			strcpy(seqs[sq][siteCnt].seq, line);
			siteCnt++;
			seqCnt++;
		}
	}
	fclose(fastaFpt);

	NEW(posIndex, seqCnt, int);
	NEW(treeIndex, seqCnt + 1, int);
	NEW(phyloSpeciesSample, seqCnt + 1, int);
	treeIndex[seqCnt] = -1;

	fastaFpt = fopen(fastaFile, "w");

	if (B->WT)
		fprintf(priorFpt, ">WEIGHT\n");

	for (seqCnt = 0, j = 0; j < RP->nMaxBlocks; j++) {
		for (i = 0; i < IP->nNumSequences; i++) {
			if (seqs[i][j].header && seqs[i][j].isPhyloSeq) {
				fprintf(fastaFpt, "%s\n", seqs[i][j].header);
				fprintf(fastaFpt, "%s\n", seqs[i][j].seq);

				if (B->WT)
					fprintf(priorFpt, "%d %9.4f\n", seqCnt + 1, B->WT->weight[i]);

				if (B->Phylo->phyloIndex) {
					treeIndex[seqCnt] = B->Phylo->phyloIndex[i];
					phyloSpeciesSample[seqCnt] = B->Phylo->phyloSpeciesSample[B->Phylo->phyloIndex[i]];
				}
				posIndex[seqCnt] = seqs[i][j].pos;
				seqCnt++;
			}
		}
	}

	for (seqCnt2 = seqCnt, j = 0; j < RP->nMaxBlocks; j++) {
		for (i = 0; i < IP->nNumSequences; i++) {
			if (seqs[i][j].header && !seqs[i][j].isPhyloSeq) {
				fprintf(fastaFpt, "%s\n", seqs[i][j].header);
				fprintf(fastaFpt, "%s\n", seqs[i][j].seq);

				if (B->WT)
					fprintf(priorFpt, "%d %9.4f\n", seqCnt2 + 1, B->WT->weight[i]);

				seqCnt2++;
			}
		}
	}

	fclose(fastaFpt);
	if (B->WT)
		fprintf(priorFpt, ">\n\n");

	if (B->Phylo->phyloIndex) {
		i = 0;
		while (i < seqCnt) {
			if (treeIndex[i] != -1) {
				treeStart = i;
				treeEnd = i;
				origTree = treeIndex[i];
				origPos = posIndex[i];
				while (treeIndex[i] == origTree && posIndex[i] == origPos) {
					treeEnd = i;
					i++;
				}
				fprintf(priorFpt, ">TREE %d %d", treeStart + 1, treeEnd + 1);
				if (phyloSpeciesSample[i])
					fprintf(priorFpt, " 1");
				fprintf(priorFpt, "\n");
				PrintTree(B, origTree, priorFpt);
				fprintf(priorFpt, ">\n\n");
			} else
				i++;
		}
	}
	for (i = 0; i < IP->nNumSequences; i++) {
		for (j = 0; j < RP->nMaxBlocks; j++) {
			if (seqs[i][j].seq) {
				free(seqs[i][j].header);
				free(seqs[i][j].seq);
			}
		}
		free(seqs[i]);
	}
	free(seqs);
	free(line);
	free(token);
	free(treeIndex);
	free(posIndex);
	free(phyloSpeciesSample);
}


void 
PrintTree(Model B, int seq, FILE * priorFpt)
{
	PhyloType       PH;
	PhyloTree       tree;

	PH = B->Phylo;

	tree = PH->phyloTree[seq];
	fprintf(priorFpt, "(");
	PrintTree2(tree, priorFpt);
	fprintf(priorFpt, ");\n");
}


void 
PrintTree2(PhyloTree tree, FILE * priorFpt)
{
	if (tree) {
		if (tree->species)
			fprintf(priorFpt, "%s:%g", tree->species, tree->length);
		else {
			fprintf(priorFpt, "(");
			PrintTree2(tree->leftChild, priorFpt);
			fprintf(priorFpt, ",");
			PrintTree2(tree->rightChild, priorFpt);
			fprintf(priorFpt, "):%g", tree->length);
		}
	}
}



void 
GetProgramName(Model B, char *name)
{
	char           *msg;

	if (B->IP->programName) {
		if (!FileExists(B->IP->programName)) {
			NEW(msg, LINE_SIZE, char);
			sprintf(msg, "%s could not be found. -align_centroid will not be run.", B->IP->programName);
			p_warning(msg);
			free(msg);
		}
		strncpy(name, B->IP->programName, FILENAME_MAX);
	} else
#if _MPI_
		strncpy(name, B->IP->argv[0], strlen(B->IP->argv[0]) - 4);	/* assumes MPI name ends
										 * in .MPI */
#else
		strcpy(name, B->IP->argv[0]);
#endif
}


void 
ConstructCommandLine(Model B, char *fastaFile, char *command, char *origPriorFile)
{
	IPtype          IP;
	int             i;

	IP = B->IP;

	strcpy(command, "'");
	GetProgramName(B, command + 1);
	strcat(command, "'");
	strcat(command, " ");
	strcat(command, fastaFile);
	strcat(command, " ");
	strcat(command, IP->argv[2]);
	strcat(command, " ");
	strcat(command, "-E 1");
	strcat(command, " ");
	strcat(command, "-S 10");
	strcat(command, " ");
	strcat(command, "-p 50");
	strcat(command, " ");
	strcat(command, "-Z");
	strcat(command, " ");
	if (B->Phylo->treeCount > 0)
		strcat(command, "-m -y ");

	for (i = 3; i < IP->argc; i++) {
		if (strcmp(IP->argv[i], "-F") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-I") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
			strcat(command, IP->argv[i + 1]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-J") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
			strcat(command, IP->argv[i + 1]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-M") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
			strcat(command, IP->argv[i + 1]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-R") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
			strcat(command, IP->argv[i + 1]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-a") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
			strcat(command, IP->argv[i + 1]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-c") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
			strcat(command, IP->argv[i + 1]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-n") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-r") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-x") == 0) {
			strcat(command, IP->argv[i]);
			strcat(command, " ");
		} else if (strcmp(IP->argv[i], "-P") == 0) {
			strcpy(origPriorFile, IP->argv[i + 1]);
		}
	}
}


void 
ParseOutputFile(Model B, char *outFile, FILE * fpt)
{
	FILE           *outFpt;
	char           *line;
	int             foundIt = FALSE;
	int             numMotifs;
	int             t;
	int             seq;
	int             seqCnt;
	LineStruct    **lineStruct;
	int             cnt;
	char           *token;
	char           *ptr;
	int             i;

	NEW(line, LINE_SIZE + 1, char);

	outFpt = fopen(outFile, "r");

	while (!foundIt && GetALine(outFpt, line, LINE_SIZE, FALSE)) {
		if (strstr(line, "MAP MAXIMIZATION RESULTS") != NULL ||
		    (B->Phylo->treeCount > 0 && strstr(line, "NEAR OPTIMAL RESULTS") != NULL))
			foundIt = TRUE;
		if (strstr(line, "ERROR") != NULL) {
			fprintf(fpt, "%s\n", line);
			return;
		}
	}
	GetALine(outFpt, line, LINE_SIZE, FALSE);

	foundIt = FALSE;
	while (!foundIt && GetALine(outFpt, line, LINE_SIZE, FALSE)) {
		if (strstr(line, "---------------------") != NULL) {
			fprintf(fpt, "%s\n", line);
			foundIt = TRUE;
		}
	}

	for (t = 0; t < B->IP->nNumMotifTypes; t++) {
		numMotifs = 0;
		while (GetALine(outFpt, line, LINE_SIZE, FALSE)) {
			fprintf(fpt, "%s\n", line);
			if (strstr(line, "Num Motifs:") != NULL) {
				NEW(token, TOKEN_SIZE, char);
				ptr = line;
				ptr = GetAToken(ptr, token);
				ptr = GetAToken(ptr, token);
				ptr = GetAToken(ptr, token);
				numMotifs = atoi(token);
				free(token);
				break;
			} else if (strstr(line, "No Motifs Detected") != NULL)
				break;
		}

		if (numMotifs) {
			NEWP(lineStruct, numMotifs, LineStruct);

			cnt = 0;
			seq = -1;
			seqCnt = 1;
			while (GetALine(outFpt, line, LINE_SIZE, FALSE)) {
				if (strstr(line, "*") != NULL)
					break;
				else {
					NEW(lineStruct[cnt], 1, LineStruct);
					NEW(lineStruct[cnt]->line, LINE_SIZE, char);
					ParseLine(B, line, &seq, &seqCnt, lineStruct[cnt]);
					cnt++;
				}
			}

			qsort(lineStruct, numMotifs, sizeof(LineStruct *), CompLineStruct);
			for (i = 0; i < numMotifs; i++)
				fprintf(fpt, "%s", lineStruct[i]->line);
			fprintf(fpt, "%s\n", line);

			for (i = 0; i < numMotifs; i++) {
				free(lineStruct[i]->line);
				free(lineStruct[i]);
			}
			free(lineStruct);
		}
	}

	fclose(outFpt);
	free(line);
}


void 
ParseLine(Model B, char *line, int *prevSeq, int *seqCnt, LineStruct * lineStruct)
{
	char           *token;
	char           *ptr;
	char          **items;
	int             i;
	int             tokenCnt;
	int             pos;
	int             start;
	int             end;
	int             seq;
	int             isRev;
	double          prob;
	char           *ptr2;
	char           *leftFlank;
	char           *rightFlank;
	char           *site;
	char           *rev;
	char           *desc;
	char           *linePtr;

	NEWP(items, 12, char);
	NEW(token, TOKEN_SIZE, char);
	NEW(desc, 51, char);
	NEW(leftFlank, FLANK + 1, char);
	NEW(rightFlank, FLANK + 1, char);

	tokenCnt = 0;
	ptr = line;
	while ((ptr = GetAToken(ptr, token))) {
		NEW(items[tokenCnt], strlen(token) + 1, char);
		strcpy(items[tokenCnt], token);
		tokenCnt++;
	}

	pos = atoi(items[tokenCnt - 1]);
	seq = atoi(items[tokenCnt - 2]);
	rev = items[tokenCnt - 3];
	isRev = (strcmp(items[tokenCnt - 3], "R") == 0);
	start = atoi(items[2]) + pos;
	prob = strtod(items[tokenCnt - 4], &ptr2);
	end = atoi(items[tokenCnt - 5]) + pos;

	lineStruct->seq = seq;
	lineStruct->pos = pos;
	linePtr = lineStruct->line;

	if (seq == *prevSeq)
		(*seqCnt)++;
	else
		*seqCnt = 1;

	sprintf(linePtr, "%4d,%3d ", seq + 1, *seqCnt);
	linePtr = lineStruct->line + strlen(lineStruct->line);
	*prevSeq = seq;

	sprintf(linePtr, " %6d", start);
	linePtr = lineStruct->line + strlen(lineStruct->line);

	if (items[3][0] == tolower(items[3][0]))
		site = items[4];
	else
		site = items[3];

	GetLeftFlank(B, seq, start - 1, isRev, strlen(site), leftFlank);
	GetRightFlank(B, seq, start - 1, isRev, strlen(site), rightFlank);

	sprintf(linePtr, " %-5s %s %-5s %6d %4.2f %s ", leftFlank, site, rightFlank, end, prob, rev);
	linePtr = lineStruct->line + strlen(lineStruct->line);

	strncpy(desc, B->IP->fastaHeader[seq], 50);
	sprintf(linePtr, "%s\n", desc);
	linePtr = lineStruct->line + strlen(lineStruct->line);

	free(token);
	for (i = 0; i < tokenCnt; i++)
		free(items[i]);
	free(items);
	free(desc);
	free(leftFlank);
	free(rightFlank);
}


void 
GetLeftFlank(Model B, int seq, int pos, int rev, int siteWidth, char *flank)
{
	int             start;
	int             len;
	char           *letters;
	int             i;
	int             j;

	start = SequenceStartPos(B, seq);
	len = SequenceLength(B, seq);
	letters = &((*B->Seq->R)[start]);

	if (rev) {
		for (j = 0, i = min(len - 1, pos + FLANK); i > pos; i--, j++)
			flank[j] = curr_ch(complement(letters[i]), B->IP->nAlphaLen, FALSE);
	} else {
		for (j = 0, i = max(0, pos - FLANK); i < pos; i++, j++)
			flank[j] = curr_ch(letters[i], B->IP->nAlphaLen, FALSE);
	}

	flank[j] = '\0';
}


void 
GetRightFlank(Model B, int seq, int pos, int rev, int siteWidth, char *flank)
{
	int             start;
	int             len;
	char           *letters;
	int             i;
	int             j;

	start = SequenceStartPos(B, seq);
	len = SequenceLength(B, seq);
	letters = &((*B->Seq->R)[start]);

	if (rev) {
		for (j = 0, i = pos - siteWidth; i >= max(0, pos - siteWidth - FLANK + 1); i--, j++)
			flank[j] = curr_ch(complement(letters[i]), B->IP->nAlphaLen, FALSE);
	} else {
		for (j = 0, i = pos + siteWidth; i < min(pos + siteWidth + FLANK, len); i++, j++)
			flank[j] = curr_ch(letters[i], B->IP->nAlphaLen, FALSE);
	}

	flank[j] = '\0';
}


void 
WriteCentroidWilcoxon(Model B, char *wilcoxFile)
{
	IPtype          IP;
	FILE           *fpt;
	int             i;
	char           *outfile;
	char           *comStr;
	int             ch;
	char           *msg;
#ifdef _CYGWIN_
	char           *comStr2;
#endif

	IP = B->IP;

	NEW(comStr, 3 * FILENAME_MAX + 32, char);

	outfile = GibbsTempnam(NULL, "wil");
	remove(outfile);

	for (i = strlen(IP->argv[0]) - 1; i >= 0; i--) {
		if (IP->argv[0][i] == '/')
			break;
	}
	strncpy(comStr, IP->argv[0], i + 1);
	if (strstr(IP->argv[0], ".x86") == NULL)
		strcat(comStr, "wilcox.test");
	else
		strcat(comStr, "wilcox.test.x86");

	if (!FileExists(comStr)) {
		NEW(msg, 3 * FILENAME_MAX + 32, char);
		sprintf(msg, "WriteWilcoxon: can't find %s.", comStr);
		p_warning(msg);
		free(msg);

		remove(outfile);
		free(comStr);
		free(outfile);

		return;
	}
#ifdef _CYGWIN_
	/* deal with spaces in path name for Cygwin */
	NEW(comStr2, strlen(comStr) + 32, char);
	strcpy(comStr2, "'");
	strcat(comStr2, comStr);
	strcat(comStr2, "'");
	strcpy(comStr, comStr2);
	free(comStr2);
#endif

	strcat(comStr, " ");
	strcat(comStr, wilcoxFile);
	strcat(comStr, " -P -L >");
	strcat(comStr, outfile);

	system(comStr);

	fpt = fopen(outfile, "r");
	if (fpt == NULL)
		p_error("Unable to open results file for Wilcoxon test.");

	fprintf(IP->Datafiles->out_fpt,
		"\n\n======================= Wilcoxon Signed-rank Test for Centroid Model ========================\n");

	while ((ch = getc(fpt)) != EOF)
		putc(ch, IP->Datafiles->out_fpt);
	fclose(fpt);

	remove(outfile);

	free(comStr);
	free(outfile);
}


void 
Credibility(Model B, int *centroidPos)
{
	IPtype          IP;
	RPType          RP;
	SampleStruct  **cent;
	SampleStruct   *sample;
	int             seq;
	int             i;
	int             pos;
	int             width;
	int             k;
	int            *dist;
	int             seed;
	int             samp;
	int             n;
	int             credPos;
	int             credibility;
	double          credLimit[] = {0.95, 0.90, 0.85};
	double          credAvg[] = {0.0, 0.0, 0.0, 0.0};
	char           *desc;
	double          sum;
	int             start;
	int             sampleCnt;
	int            *centCnt;

	IP = B->IP;
	RP = B->RP;

	NEW(desc, 51, char);

	NEWP(cent, IP->nNumSequences, SampleStruct);
	NEW(centCnt, IP->nNumSequences, int);
	for (seq = 0; seq < IP->nNumSequences; seq++) {
		NEW(cent[seq], SequenceLength(B, seq), SampleStruct);
		start = SequenceStartPos(B, seq);
		centCnt[seq] = 0;
		for (i = 0; i < SequenceLength(B, seq); i++) {
			if (centroidPos[i + start] == CENT_START) {
				width = 1;
				pos = i;
				i++;
				while (centroidPos[i + start] == CENT_POS) {
					width++;
					i++;
				}
				cent[seq][centCnt[seq]].pos = pos;
				cent[seq][centCnt[seq]].width = width;
				centCnt[seq]++;
			}
		}
	}

	fprintf(IP->Datafiles->out_fpt, "\nCredibility limits\n");
	fprintf(IP->Datafiles->out_fpt, "%% samples within specified number of mismatched sites\n");
	fprintf(IP->Datafiles->out_fpt, "from the sequence centroid.\n");
	fprintf(IP->Datafiles->out_fpt, "                                 Avg. dist.\n");
	fprintf(IP->Datafiles->out_fpt, "Seq.      95%%     90%%     85%%    to samples\n");
	fprintf(IP->Datafiles->out_fpt, "-------------------------------\n");

	NEW(sample, RP->nMaxBlocks, SampleStruct);
	NEW(dist, IP->nSeeds * IP->bayesSamplePeriod, int);
	for (seq = 0; seq < IP->nNumSequences; seq++) {
		bzero(dist, IP->nSeeds * IP->bayesSamplePeriod * sizeof(int));
		bzero(sample, RP->nMaxBlocks * sizeof(SampleStruct));

		for (n = 0, seed = 0; seed < IP->nSeeds; seed++) {
			for (samp = 0; samp < IP->bayesSamplePeriod; samp++) {
				for (sampleCnt = 0, k = 0; k < RP->nMaxBlocks; k++) {
					if (RP->samples[seq][seed][samp][k].width) {
						sample[sampleCnt].pos = RP->samples[seq][seed][samp][k].pos;
						sample[sampleCnt].width = RP->samples[seq][seed][samp][k].width;
						sampleCnt++;
					}
				}
				dist[n] = SampleDist(cent[seq], centCnt[seq], sample, sampleCnt);
				n++;
			}
		}

		qsort(dist, IP->nSeeds * IP->bayesSamplePeriod, sizeof(int), compare_ints);

		fprintf(IP->Datafiles->out_fpt, "%-5d", seq + 1);
		for (i = 0; i < 3; i++) {
			credPos = credLimit[i] * (double) (IP->nSeeds * IP->bayesSamplePeriod);
			credibility = dist[credPos];
			fprintf(IP->Datafiles->out_fpt, "\t%5d", credibility);
			credAvg[i] += credibility;
		}

		for (sum = 0, i = 0; i < IP->nSeeds * IP->bayesSamplePeriod; i++)
			sum += dist[i];
		sum /= (double) (IP->nSeeds * IP->bayesSamplePeriod);
		credAvg[3] += sum;

		strncpy(desc, B->IP->fastaHeader[seq], 50);
		fprintf(IP->Datafiles->out_fpt, "\t%5.2f\t%s\n", sum, desc);
	}

	fprintf(IP->Datafiles->out_fpt, "Avg.\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n",
	     credAvg[0] / IP->nNumSequences, credAvg[1] / IP->nNumSequences,
	    credAvg[2] / IP->nNumSequences, credAvg[3] / IP->nNumSequences);
	fprintf(IP->Datafiles->out_fpt, "\n");

	FREEP(cent, IP->nNumSequences);
	free(desc);
	free(dist);
	free(sample);
	free(centCnt);
}


int 
SampleDist(SampleStruct * centroid, int centCnt, SampleStruct * sample, int sampleCnt)
{
	int             dist = 0;
	int             i;
	int             j;
	int             roverlap;
	int             over;
	int             overlap;

	for (i = 0; i < centCnt; i++) {
		over = 0;
		for (j = 0; j < sampleCnt; j++) {
			roverlap = max(centroid[i].width, sample[j].width) / 2 + 1;
			overlap = OverlapDist(centroid[i].pos, centroid[i].width,
					    sample[j].pos, sample[j].width);
			if (overlap >= roverlap) {
				over = 1;
				break;
			}
		}
		dist += (1 - over);
	}

	for (i = 0; i < sampleCnt; i++) {
		over = 0;
		for (j = 0; j < centCnt; j++) {
			roverlap = max(sample[i].width, centroid[j].width) / 2 + 1;
			overlap = OverlapDist(sample[i].pos, sample[i].width,
					centroid[j].pos, centroid[j].width);
			if (overlap >= roverlap) {
				over = 1;
				break;
			}
		}
		dist += (1 - over);
	}

	return dist;
}


int 
OverlapDist(int start1, int width1, int start2, int width2)
{
	int             end1;
	int             end2;
	int             dist = 0;

	end1 = start1 + width1 - 1;
	end2 = start2 + width2 - 1;

	if (start1 <= start2 && start2 <= end1 && end1 <= end2)
		dist = end1 - start2 + 1;
	else if (start1 <= start2 && end1 >= end2)
		dist = width2;
	else if (start1 >= start2 && end1 <= end2)
		dist = width1;
	else if (start1 >= start2 && start1 <= end2 && end1 >= end2)
		dist = end2 - start1 + 1;

	return (dist);
}


int 
compare_ints(const void *a, const void *b)
{
	int            *arg1 = (int *) a;
	int            *arg2 = (int *) b;

	if (*arg1 < *arg2)
		return -1;
	else if (*arg1 == *arg2)
		return 0;
	else
		return 1;
}


int 
CompLineStruct(const void *p1, const void *p2)
{
	LineStruct     *ptr1 = *(LineStruct **) p1;
	LineStruct     *ptr2 = *(LineStruct **) p2;

	if (ptr1->seq < ptr2->seq)
		return -1;
	else if (ptr1->seq == ptr2->seq) {
		if (ptr1->pos < ptr2->pos)
			return -1;
		else if (ptr1->pos == ptr2->pos)
			return 0;
		else
			return -1;
	} else
		return 1;
}


void 
DumpSamples(Model B)
{
	IPtype          IP;
	RPType          RP;
	char           *file = "/ccmb/data/thompson/taf4b/output/samples.dump";
	FILE           *fpt;
	int             seq;
	int             seed;
	int             samp;
	int             k;

	IP = B->IP;
	RP = B->RP;

	NEW(file, FILENAME_MAX, char);
	if (IP->Datafiles->output_filename)
		sprintf(file, "%s.%d.samples", IP->Datafiles->output_filename, IP->nRank);
	else
		strcpy(file, "/ccmb/data/thompson/taf4b/output/samples.dump");

	fpt = fopen(file, "w");

	for (seq = 0; seq < IP->nNumSequences; seq++) {
		for (seed = 0; seed < IP->nSeeds; seed++) {
			for (samp = 0; samp < IP->bayesSamplePeriod; samp++) {
				for (k = 0; k < RP->nMaxBlocks; k++) {
					if (RP->samples[seq][seed][samp][k].width)
						fprintf(fpt, "%d %d %d %d %d %d\n", seq, seed, samp,
							RP->samples[seq][seed][samp][k].pos,
							RP->samples[seq][seed][samp][k].width,
							RP->samples[seq][seed][samp][k].t);
				}
			}
		}
	}
	fclose(fpt);
	free(file);
}


void 
DumpDistance(Model B, int seq, int *dist)
{
	IPtype          IP;
	char           *file = "/ccmb/data/thompson/taf4b/output/distance.dump";
	FILE           *fpt;
	int             seed;
	int             samp;
	int             n;

	IP = B->IP;

	if (seq == 0)
		fpt = fopen(file, "w");
	else
		fpt = fopen(file, "a");

	for (n = 0, seed = 0; seed < IP->nSeeds; seed++) {
		for (samp = 0; samp < IP->bayesSamplePeriod; samp++) {
			fprintf(fpt, "%d %d %d %d\n", seq, seed, samp, dist[n]);
			n++;
		}
	}

	fclose(fpt);
}


int 
GetProcInfo()
{
	int             pid;
	char            fname1[FILENAME_MAX];
	char            buffer[FILENAME_MAX];
	FILE           *procfile;
	int             itmp;
	char            ctmp;

	pid = getpid();

	sprintf(fname1, "/proc/%d/stat", pid);
	if ((procfile = fopen(fname1, "r")) > 0) {
		fprintf(stdout, "=====================================\n");

		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "Process: %d\n", itmp);
		if (fscanf(procfile, "%s ", &buffer[0]) > 0)
			fprintf(stdout, "Exec file: %s\n", buffer);
		if (fscanf(procfile, "%c ", &ctmp) > 0)
			fprintf(stdout, "State: %c\n", ctmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "Parent process: %d\n", itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "Process group: %d\n", itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "Session id: %d\n", itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "TTY: %d\n", itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "TTY owner process group: %d\n", itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Flags: 0x%x\n", itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Minor faults (no memory page): %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Minor faults, children: %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Major faults (memory page faults): %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Major faults, children: %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "utime: %d\n", itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "stime: %d\n", itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "utime, children: %d\n", itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "stime, children: %d\n", itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "jiffies remaining in current time slice: %d\n",
				itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "'nice' value: %d\n", itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "jiffies until next timeout: %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "jiffies until next SIGALRM: %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%d ", &itmp) > 0)
			fprintf(stdout, "start time (jiffies since system boot): %d\n",
				itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Virtual memory size: %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Resident set size: %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "rlim: %u\n",
				(unsigned int) itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Start of text: 0x%x\n", itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "End of text: 0x%x\n", itmp);
		if (fscanf(procfile, "%u ", &itmp) > 0)
			fprintf(stdout, "Start of stack: 0x%x\n", itmp);
		fclose(procfile);
		fprintf(stdout, "=====================================\n");

		return 0;
	} else
		return -1;
}
