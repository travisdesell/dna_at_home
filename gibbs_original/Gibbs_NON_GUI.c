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
/****************************************************************/
/* */
/* $Id: /home/junzhu1/proj/RCS/Gibbs_GUI.c,v 1.2                */
/* 1996/09/18 18:44:26 junzhuExp junzhu $                       */
/* */
/* $Author: Bill $                                            */
/* $Date: 2007/05/23 18:19:55 $                                 */
/* Description: non_graphic user interface for Gibbs,           */
/* working for batch job                           */
/****************************************************************/

#include "Gibbs_NON_GUI.h"

void 
Gibbs_NON_GUI(int argc, char **argv)
{
#ifndef _MPI_
	int             i;
	Model           M;
#endif

	/* allocate space for GUI_argv, we only allocate array of pointer,    */
	/* we do not need to locate space for each array pointed by pointers. */
	/* The pointers will accept values from other array, such as          */
	/* GUI_argv[i]=argv[i+1]                                             */
	NEWP(GUI_argv, MAX_ARG_SIZE, char);	/* BT 1/2/98 */

	if (strcmp(argv[1], "-PBernoulli") == 0) {	/* bernoulli method */
		/* construct command line */
#ifdef _MPI_
		Gibbs_MPI(argc, argv);	/* in case some uses -PBernoulli
					 * instead of -PMBernoulli */
#else
		GUI_argc = argc - 1;
		for (i = 1; i <= GUI_argc - 1; i++) {
			GUI_argv[i] = argv[i + 1];
		}

		GUI_argv[0] = argv[0];

		/* call the main function of Bernoulli */
		M = alloc_Model();
		stripargs(GUI_argc, GUI_argv, M, FALSE);	/* Get arguments from
								 * command line               */

		get_inputs(M);

		if (M->IP->is_defined[cl_l])
			Wilcoxon(M);
		else
			get_strings(M);	/* Read the sequences */

		InitRProbStruct(M);

		alloc_Counts(M);
		set_counts(M);	/* Set the initial counts     */
		set_posterior_prob(M->IP, M->First);	/* Set initial posterior
							 * prob */
		*M->Seq->ProcessedSTR = Xnu_Sequence(M);	/* preprocess the
								 * sequence    */

		set_pseudo_counts(M, M->IP, M->First);	/* Set initial pseudo
							 * counts  */

#ifndef _MPI_
		print_options(M);	/* BT 1/27/97 */
		PrintSeqDescriptions(M);
#endif

		if (M->IP->is_defined[cl_B])
			ReadBkgndComposition(M, M->First);

		if (M->IP->is_defined[cl_U])
			ReadSpacingFile(M, M->IP->Datafiles->spacing_fpt);

		if (M->IP->Datafiles->weight_fpt != NULL)
			ReadWeightFile(M);

		CalcAlignProbFromFootprint(M);

		if (!M->IP->is_defined[cl_G])
			find_best_sites(M);	/* find maximal alignment     */
		else
			GroupSample(M);

		FreeData(M);
#endif				/* not MPI */
	}
#ifdef _MPI_
	else if (strcmp(argv[1], "-PMBernoulli") == 0) {
		Gibbs_MPI(argc, argv);
	} else if (strncasecmp(argv[1], "-Hier", 5) == 0)
		Gibbs_Hierarchy(argc, argv);
#endif
	else {
#ifdef _MPI_
		Gibbs_MPI(argc, argv);	/* in case some uses -PBernoulli
					 * instead of -PMBernoulli */
#else

		GUI_argc = argc;
		for (i = 0; i <= GUI_argc - 1; i++) {
			GUI_argv[i] = argv[i];
		}

		/* call the main function of Bernoulli */
		M = alloc_Model();
		stripargs(GUI_argc, GUI_argv, M, FALSE);	/* Get arguments from
								 * command line               */
		get_inputs(M);

		if (M->IP->is_defined[cl_l])
			Wilcoxon(M);
		else
			get_strings(M);	/* Read the sequences */

		InitRProbStruct(M);

		alloc_Counts(M);
		set_counts(M);	/* Set the initial counts     */
		set_posterior_prob(M->IP, M->First);	/* Set initial posterior
							 * prob */
		*M->Seq->ProcessedSTR = Xnu_Sequence(M);	/* preprocess the
								 * sequence    */

		set_pseudo_counts(M, M->IP, M->First);	/* Set initial pseudo
							 * counts  */

#ifndef _MPI_
		print_options(M);	/* BT 1/27/97 */
		PrintSeqDescriptions(M);
#endif

		if (M->IP->is_defined[cl_B])
			ReadBkgndComposition(M, M->First);

		if (M->IP->is_defined[cl_U])
			ReadSpacingFile(M, M->IP->Datafiles->spacing_fpt);

		if (M->IP->Datafiles->weight_fpt != NULL)
			ReadWeightFile(M);

		CalcAlignProbFromFootprint(M);

		if (!M->IP->is_defined[cl_G])
			find_best_sites(M);	/* find maximal alignment     */
		else
			GroupSample(M);

		FreeData(M);
#endif
	}

	free(GUI_argv);
}
