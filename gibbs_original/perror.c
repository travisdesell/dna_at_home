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
#include "perror.h"		/* BT 4/25/97 */


void 
p_error(char *msg)
{
	char            tmp[512];

	sprintf(tmp, "\nERROR :: %s\nEXITING\n", msg);	/* BT 4/25/97 */
	printf("%s", tmp);
#ifdef _GUI_
	if (GL_GUI) {
		message(tmp);
		XFlush(idispl);
	}
#endif

#ifdef _MPI_
	MPI_Abort(MPI_COMM_WORLD, -1);
#endif
	exit(1);
}


void 
p_internal_error(char *msg)
{
	char            tmp[512];

	sprintf(tmp, "\nFATAL INTERNAL ERROR :: %s\nEXITING\n", msg);	/* BT 4/25/97 */
	printf("%s", tmp);
	printf("\nPlease contact gibbsamp@brown.edu if you cannot resolve this error.\n");
#ifdef _GUI_
	if (GL_GUI) {
		message(tmp);
		XFlush(idispl);
	}
#endif

#ifdef _MPI_
	MPI_Abort(MPI_COMM_WORLD, -1);
#endif
	exit(1);
}


void 
p_usage_error(char *msg)
{
	printf("\nERROR :: %s\n\n", msg);
	print_usage_bernoulli();
	printf("EXITING\n\n");
	exit(0);
}


void 
killit(int n)
{
	char            tmp[512];

	sprintf(tmp, "\nERROR ::\nEXITING\n %d", n);	/* BT 4/25/97 */
	printf("%s", tmp);

	abort();
}


void 
p_warning(char *msg)
{
	char            tmp[512];

	sprintf(tmp, "\nWARNING :: %s\n", msg);
	printf("%s", tmp);
#ifdef _GUI_
	if (GL_GUI) {
		message(tmp);
		XFlush(idispl);
	}
#endif
}


void 
ErrorMsg(char *msg, char *filename, int line)
{
	char            tmp[FILENAME_MAX + 1024];

	sprintf(tmp, "\nERROR :: %s in %s at line %d \nEXITING\n", msg, filename, line);
	printf("%s", tmp);
#ifdef _GUI_
	if (GL_GUI) {
		message(tmp);
		XFlush(idispl);
	}
#endif

#ifdef _MPI_
	MPI_Abort(MPI_COMM_WORLD, -1);
#endif
	exit(1);
}
