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
#include "exchange.h"

void CopyMaxResults( MaxResults *destMax, MaxResults *srcMax, Model B )
{
  IPtype IP;
  int    maxlen;
  int    t;
  int    i;
  int    j;
  
  IP = B->IP;
 
  maxlen = srcMax->nMaxLen;
  NEW( destMax->nNumMotifs, IP->nNumMotifTypes, int);
  NEW(destMax->nMotifLen, IP->nNumMotifTypes, int);     /* BT 5/28/97 */
  NEW(destMax->dMap, IP->nNumMotifTypes, double);       /* BT 3/4/98 */
  NEWP(destMax->nMotifLoc, maxlen, int);
  NEWP(destMax->dvMotifProb, maxlen, double);
  NEWP(destMax->RevComp, maxlen, short);
  NEW(destMax->F, 1, FragStruct);

   for (i = 0; i < maxlen; i++) 
     {
       NEW(destMax->nMotifLoc[i], IP->nNumMotifTypes, int);
       NEW(destMax->dvMotifProb[i], IP->nNumMotifTypes, double);
       NEW(destMax->RevComp[i], IP->nNumMotifTypes, short);
     }

  destMax->F->nMaxLen = NULL;
  destMax->F->nvFragCnts = NULL;

  if(!IP->is_defined[cl_F]) 
    {
      NEW(destMax->F->FragWidth, IP->nNumMotifTypes, int);
      NEW(destMax->F->nOldFragWidth, IP->nNumMotifTypes, int);
      NEWP(destMax->F->nColMask, IP->nNumMotifTypes, int);
      NEWP(destMax->F->fragPos, IP->nNumMotifTypes, int);
      for(t = 0; t < IP->nNumMotifTypes; t++)  
	{
	  destMax->F->FragWidth[t] = srcMax->F->FragWidth[t];
	  destMax->F->nOldFragWidth[t] = srcMax->F->FragWidth[t];
	  NEW( destMax->F->nColMask[t], srcMax->F->nMaxLen[t] , int);
	  NEW( destMax->F->fragPos[t], IP->nMotifLen[t] , int);
	  for(i = 0; i <  srcMax->F->nMaxLen[t]; i++) 	    
	    destMax->F->nColMask[t][i] = srcMax->F->nColMask[t][i];
	  for(i = 0; i < IP->nMotifLen[t]; i++) 	    
	    destMax->F->fragPos[t][i] = srcMax->F->fragPos[t][i];
	}
    }

  destMax->dProbability = srcMax->dProbability;
  destMax->nIterationNum = srcMax->nIterationNum;
  destMax->nseed = srcMax->nseed;

  for(t = 0; t < IP->nNumMotifTypes; t++) 
    {
      destMax->nNumMotifs[t] = srcMax->nNumMotifs[t];
      destMax->nMotifLen[t] = srcMax->nMotifLen[t];     
    }

   for(t=0;t < IP->nNumMotifTypes; t++)
     {
       for( j = 0; j < destMax->nNumMotifs[t]; j++ )
	 {
	   destMax->nMotifLoc[j][t] = srcMax->nMotifLoc[j][t];
	   destMax->RevComp[j][t]= srcMax->RevComp[j][t];
	   destMax->dvMotifProb[j][t] = srcMax->dvMotifProb[j][t];
	 }
        destMax->dMap[t] = srcMax->dMap[t];
     }
}


#ifdef _MPI_

void ExchangeTemps( Model B, double currProb, int nSeed, int iter )
{
  IPtype      IP;
  AnlType     AN;
  MPI_Status  status;
  double      dTempInfo[4];
  struct tm   *currTime;
  time_t      tloc;

  IP = B->IP;
  AN = B->AN;

  dTempInfo[0] = AN->currTemp;
  dTempInfo[1] = currProb;
  dTempInfo[2] = (double) nSeed;
  dTempInfo[3] = (double) TotalNumMotifs( B );
  time( &tloc );
  currTime = localtime( &tloc );
  PrintTempOut( IP->Datafiles->mpiTemp_fpt, 
		"\n%d is sending temp = %f prob = %f seed = %d motifs = %d at %d:%d:%d iter = %d\n",
		IP->nRank, dTempInfo[0], dTempInfo[1], (int) dTempInfo[2], (int) dTempInfo[3], 
		currTime->tm_hour, currTime->tm_min, currTime->tm_sec, iter );
  Gibbs_MPI_Send( B, dTempInfo, 4, MPI_DOUBLE, 0, G_MPI_TEMP, MPI_COMM_WORLD );
  
  PrintTempOut( IP->Datafiles->mpiTemp_fpt, "Waiting on %d to send new temp data at %d:%d:%d\n", 0,
		currTime->tm_hour, currTime->tm_min, currTime->tm_sec );

  Gibbs_MPI_Recv( B, dTempInfo, 4, MPI_DOUBLE, 0, G_MPI_TEMP, 
	    MPI_COMM_WORLD, &status);
  time( &tloc );
  currTime = localtime( &tloc );
  PrintTempOut(  IP->Datafiles->mpiTemp_fpt, "%d has sent temp = %f prob = %f seed = %d at %d:%d:%d\n", 
		 0, dTempInfo[0], dTempInfo[1], (int) dTempInfo[2],
		 currTime->tm_hour, currTime->tm_min, currTime->tm_sec );
  
  AN->currTemp = dTempInfo[0];
}
#endif


void DumpMotifsForCluster( Model B, PoSition **Pos, int iter, int nSeed, double currProb )
{
  int      t;
  int      i;
  int      j;
  int      width;
  IPtype   IP;
  char     fileName[] = "/home/thompson/Muscle/Hum.Mouse/Pairs/3kb/mask/050103/bmc/cluster.test.120303.1.txt";  /* temporary for testing */
  FILE     *fpt;
  PoSition *pos_t;
  char     **desc;
  char     *seqStr;
  int      seqStart;
  int      siteStart;
  int      seqEnd;
  int      siteEnd;
  int      offset;
  char     dummyChars[4] = "ATCG";
  int      seq;

  IP = B->IP;

  fpt = fopen( fileName, "a" );

  NEWP( desc, B->IP->nNumSequences, char );     /* BT 5/21/97 */
  for( j = 0; j < B->IP->nNumSequences; j++ )
    {
      NEW( desc[j], 64, char );
      strncpy( desc[j], B->IP->fastaHeader[j], 50 );
    }
   
  for(t = 0; t < IP->nNumMotifTypes; t++) 
    {
      width = MotifWidth( B, t );
      pos_t = Pos[t];
      NEW( seqStr, width + 11, char );

      for( i = 0; i < IP->nSeqLen; i ++ )
	{
	  if( pos_t[i].nMotifStartPos )
	    {
	      seq = SequenceFromPosition( B, i );
	      fprintf( fpt, ">motif-%d-%d-%d-%d N %d %d test %c %d %g %s\n", 
		       t, nSeed, iter, IP->nRank, width, width + 10, 
		       (pos_t[i].RevComp ? 'R' : 'F' ),
		       seq, currProb, 
		       desc[seq] );

	      offset = 0;
	      seqStart = SequenceStartPos( B, seq );
	      seqEnd = SequenceEndPos( B, seq );

	      if( ! pos_t[i].RevComp )
		{
		  siteStart = i - 5;
		  siteEnd = i + width + 5 - 1;
		  
		  for( j = siteStart; j <= siteEnd; j++ )
		    {		  
		      if( j < seqStart )		    
			seqStr[offset] = dummyChars[j % 4];
		      else if( j > seqEnd )
			seqStr[offset] = dummyChars[j % 4];
		      else if( j < i || j >= i + width )
			seqStr[offset] = curr_ch((*B->Seq->R)[j], B->IP->nAlphaLen, FALSE);
		      else if( j >= i && j <= i + width - 1 )
			seqStr[offset] = curr_ch((*B->Seq->R)[j], B->IP->nAlphaLen, TRUE);
		      offset++;
		    }
		}
	      else
		{
		  siteStart = i + width + 5 - 1;
		  siteEnd = i - 5;

		  for( j = siteStart; j >= siteEnd; j-- )
		    {		  
		      if( j < seqStart )		    
			seqStr[offset] = dummyChars[j % 4];
		      else if( j > seqEnd )
			seqStr[offset] = dummyChars[j % 4];
		      else if( j < i || j >= i + width )
			seqStr[offset] = curr_ch( complement( (*B->Seq->R)[j] ), B->IP->nAlphaLen, FALSE);
		      else if( j >= i && j <= i + width - 1 )
			seqStr[offset] = curr_ch( complement( (*B->Seq->R)[j] ), B->IP->nAlphaLen, TRUE);
		      offset++;
		    }
		}
	      seqStr[offset+1] = '\0';

	      fprintf( fpt, "%s\n", seqStr );		  
	    }
	}
      fprintf( fpt, "\n" );		  
      free( seqStr );
    }

  fclose( fpt );

  FREEP( desc, IP->nNumSequences );   
}

