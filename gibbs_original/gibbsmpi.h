/****************************************************************/
/* Gibbs - A program for detecting subtle sequence signals      */
/*                                                              */
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
/*                                                              */
/* Copyright (C) 2006   Health Research Inc.                    */
/* HEALTH RESEARCH INCORPORATED (HRI),                          */
/* ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.             */
/* Email:  gibbsamp@wadsworth.org                               */
/*                                                              */
/****************************************************************/
/*                                                              */
/* Changes Copyright (C) 2007   Brown University                */
/* Brown University                                             */
/* Providence, RI 02912                                         */
/* Email:  gibbs@brown.edu                                      */
/*                                                              */
/* For the Centroid sampler, please site,                       */
/* Thompson, W.A., Newberg, L., Conlan, S.P., McCue, L.A. and   */
/* Lawrence, C.E. (2007) The Gibbs Centroid Sampler             */
/* Nucl. Acids Res., doi:10.1093/nar/gkm265                     */
/*                                                              */
/* For the Phylogenetic Gibbs Sampler, please site,             */
/* Newberg, L., Thompson, W.A., Conlan, S.P., Smith, T.M.,      */
/* McCue, L.A. and Lawrence, C.E. (2007) A phylogenetic Gibbs   */
/* sampler that yields centroid solutions for cis regulatory    */
/* site prediction., Bioinformatics,                            */
/* doi:10.1093/bioinformatics/btm241.                           */
/*                                                              */
/****************************************************************/
/*                                                              */
/* This program is free software; you can redistribute it       */
/* and/or modify it under the terms of the GNU General Public   */
/* License as published by the Free Software Foundation;        */
/* either version 2 of the License, or (at your option)         */
/* any later version.                                           */
/*                                                              */
/* This program is distributed in the hope that it will be      */
/* useful, but WITHOUT ANY WARRANTY; without even the implied   */
/* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR      */
/* PURPOSE. See the GNU General Public License for more         */
/* details.                                                     */
/*                                                              */
/* You should have received a copy of the GNU General Public    */
/* License along with this program; if not, write to the        */
/* Free Software Foundation, Inc., 51 Franklin Street,          */
/* Fifth Floor, Boston, MA  02110-1301, USA.                    */
/****************************************************************/
#ifndef _GIBBS_MPI
#define _GIBBS_MPI


#include <stdio.h>
#include <unistd.h>
#include <stdarg.h>
#include "common.h"
#ifdef _GUI_
#include "GUI_common.h"
#endif
#include "mem_mgmt.h"
#include "process.h"
#include "recur.h"
#include "prior.h"
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <time.h>
#include <fcntl.h>
#include "sampling.h"
#include "parameters.h"
#include "motif.h"
#include "bayes.h"

#include <mpi.h>
#include <math.h>
#include <stdlib.h>

#ifndef _RPI_
#include <mpe.h>     /* BT 12/09/99 */
#endif

#define G_MPI_DONE             0
#define G_MPI_SEED             1
#define G_MPI_DATA             2
#define G_MPI_TEMP             4
#define G_MPI_COUNT_DATA       8
#define G_MPI_COUNT_DATA0      16
#define G_MPI_FINAL_TEMP       32
#define G_MPI_RAND_SEED        64
#define G_MPI_COMMAND_LINE     128
#define G_MPI_SUBOPT_DONE      256
#define G_MPI_FINISH           512
#define G_MPI_ALL_DONE         1024
#define G_MPI_QUIT             2048
#define G_MPI_CONTINUE         4096
#define G_MPI_SEED_DONE        8192
#define G_MPI_BAYES            16384
#define G_MPI_BAYES_DATA       32768
#define G_MPI_BAYES_SAMPLE     65536
#define G_MPI_BAYES_MH         131072 
#define G_MPI_MESSAGE_SIZE     100
#define G_MPI_TRANSFER_SIZE    1024*1024

#define EXCHNG_COUNT       1

typedef struct MPI_temp
{
  double          temp;
  double          currProb;
  double          seed;
  int             source;
} MPITemp;

#ifndef localtime_r
extern struct tm *localtime_r(const time_t *, struct tm *);
#endif

void Gibbs_MPI( int argc, char **argv );
int Gibbs_MPI_Send( Model B, void *buf, int count, MPI_Datatype datatype, int dest, 
		    int tag, MPI_Comm comm );
int Gibbs_MPI_Recv( Model B, void *buf, int count, MPI_Datatype datatype, int source, 
		    int tag, MPI_Comm comm, MPI_Status *status );

#endif
