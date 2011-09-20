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
/***************************************************************************/
/* $Id: mem_mgmt.c,v 1.17 2009/04/23 18:43:54 Bill Exp $           */
/*                                                                         */
/* Author:         Eric C. Rouchka, 1996.6                                 */
/*                 Jun Zhu, 1996.8                                         */
/*		   Bill Thompson 1/27/97                                   */ 
/*                                                                         */
/***************************************************************************/

#include "mem_mgmt.h"

/******************* alloc_Model ***************************************/
/*                                                                     */
/* DESCRIPTION : Allocates the memory for the basic Model and returns  */
/*               a structure                                           */
/***********************************************************************/

Model alloc_Model(void)

{
   Model B;

   NEW(B, 1, Modelstruct);
   NEW(B->IP, 1, InputParams);
   NEW(B->C, 1, Counts);
   NEW(B->First, 1, Counts);
   NEW(B->Seq, 1, Stringstruct);
   NEWP(B->Seq->R, 1, char);
   NEWP(B->Seq->nvEndLocs, 1, int);
   NEWP(B->Seq->ProcessedSTR, 1, char); 
   NEW( B->Phylo, 1, PhyloStruct );
   return B;
}
 
/********************  alloc_Counts   *******************************/

void alloc_Counts(Model B)

{
   int    i,j,t,maxlen;
   IPtype IP;           /* for quick access data */
   Ctype  First, C;     /* for quick access data */
   int    k;
 
   IP=B->IP;
   First=B->First;
   C=B->C;

   /* maximum motif length */
   maxlen = findMaxMotifLen(IP);

   /* C->dSumPseudo[MotifType][BG/MOTIF][nAlphaLen] */
   NEWPP(C->dSumPseudo,  IP->nNumMotifTypes, double);               
   NEWPP(First->dSumPseudo, IP->nNumMotifTypes, double);

   /* C->dvPseudoCounts[MotifType][PositionInMotif][nAlphaLen]  */
   NEWPP(C->dPseudoCounts,IP->nNumMotifTypes , double);
   NEWPP(First->dPseudoCounts,IP->nNumMotifTypes, double);

   for(t=0;t< IP->nNumMotifTypes;t++){
     NEWP(First->dSumPseudo[t],2,double);
     NEWP(C->dSumPseudo[t],2,double);

     NEW(C->dSumPseudo[t][BG], IP->nAlphaLen, double);  
     NEW(C->dSumPseudo[t][MOTIF], IP->nAlphaLen, double);
     NEW(First->dSumPseudo[t][BG], IP->nAlphaLen, double); 
     NEW(First->dSumPseudo[t][MOTIF], IP->nAlphaLen, double);
     
     NEWP(C->dPseudoCounts[t],maxlen+1, double);
     NEWP(First->dPseudoCounts[t],maxlen+1, double);    
     for(j = 0; j <= maxlen; j++) {
         NEW(C->dPseudoCounts[t][j], IP->nAlphaLen, double);
         NEW(First->dPseudoCounts[t][j],IP->nAlphaLen, double);
      }
   }
   First->nTotBack = 0;
   NEWP(First->dmodel_sites, IP->nNumMotifTypes, double);
   NEWP(C->dmodel_sites, IP->nNumMotifTypes, double);
   NEWP(First->dmodel_pseudo, IP->nNumMotifTypes, double);
   NEWP(C->dmodel_pseudo, IP->nNumMotifTypes, double);
   NEW(First->dtot_sites, IP->nNumMotifTypes, double);
   NEW(C->dtot_sites, IP->nNumMotifTypes, double);
   NEW(First->dtot_pseudo, IP->nNumMotifTypes, double);
   NEW(C->dtot_pseudo, IP->nNumMotifTypes, double);
   NEW(First->dTot, IP->nNumMotifTypes, double);
   NEW(C->dTot, IP->nNumMotifTypes, double);
   NEW(First->dbg_pseudo, IP->nNumMotifTypes, double);
   NEW(C->dbg_pseudo, IP->nNumMotifTypes, double);

   /* First->fCounts[MotifType][PositionInMotif][nAlphabet] */
   NEWPP(First->fCounts,IP->nNumMotifTypes, double);
   NEWPP(First->wCounts,IP->nNumMotifTypes, double);
   NEWPP(C->fCounts, IP->nNumMotifTypes, double);
   NEWPP(C->wCounts, IP->nNumMotifTypes, double);
   for(t=0;t<IP->nNumMotifTypes;t++)
     {
       NEWP(First->fCounts[t],maxlen+1, double);
       NEWP(First->wCounts[t],maxlen+1, double);
       NEWP(C->fCounts[t], maxlen+1, double);
       NEWP(C->wCounts[t], maxlen+1, double);
       for(i = 0; i <= maxlen; i++) 
	 { 
      /* we also allocate space for symbol represent any A.A's 
         or nucleotides                                          */
	   NEW(First->fCounts[t][i], IP->nAlphaLen+1, double);
	   NEW(First->wCounts[t][i], IP->nAlphaLen+1, double);
	   NEW(C->fCounts[t][i], IP->nAlphaLen+1, double);
	   NEW(C->wCounts[t][i], IP->nAlphaLen+1, double);
	 }
     }
   for(t = 0; t < IP->nNumMotifTypes; t++) {
      NEW(First->dmodel_sites[t], IP->RevComplement + 1, double);
      NEW(C->dmodel_sites[t], IP->RevComplement + 1, double);
      NEW(First->dmodel_pseudo[t], IP->RevComplement + 1, double);
      NEW(C->dmodel_pseudo[t], IP->RevComplement + 1, double);
   }

   IP->nMaxSeqLen = -INT_MAX;
   for( i = 0; i < IP->nNumSequences; i++ )
     {
       if( SequenceLength( B, i ) > IP->nMaxSeqLen )
	 IP->nMaxSeqLen = SequenceLength( B, i );
     }

   if( IP->is_defined[cl_Q] )    
     {
       NEWP3( IP->nAlignCnts, IP->nSeeds + 1, int );   
       for( k = 0; k <= IP->nSeeds; k++ )
	 {
	   NEWPP( IP->nAlignCnts[k], IP->nNumMotifTypes, int );   
	   for( j = 0; j < IP->nNumMotifTypes; j++ )
	     {
	       NEWP( IP->nAlignCnts[k][j], IP->nNumSequences, int );   
	       for( i = 0; i < IP->nNumSequences; i++ )
		 NEW( IP->nAlignCnts[k][j][i], SequenceLength( B, i ), int );
	     }
	 }
     }

   if( IP->is_defined[cl_B] )
     {
       NEW( B->BP, 1, BkgndStruct );
       NEWPP( B->BP->dBkgndProb, IP->nNumSequences, double );   
       for( i = 0; i < IP->nNumSequences; i++ )
	 {
	   NEWP( B->BP->dBkgndProb[i], SequenceLength( B, i ), double );
	   for( j = 0; j < SequenceLength( B, i ); j++ )
	     NEW( B->BP->dBkgndProb[i][j], IP->nAlphaLen, double );
	 }
     }
}


/*****************  init_maxdata ********************/

void init_maxdata(MaxResults *tmp)
{
   tmp->F = NULL;
   tmp->nNumMotifs = NULL;
   tmp->nMotifLoc = NULL;
   tmp->RevComp = NULL;
   tmp->dvMotifProb = NULL;
   tmp->frequency = NULL;
   tmp->nMotifLen = NULL;    /* BT 5/23/97 */
   tmp->dMap = NULL;         /* 3/4/98 */
   tmp->nMaxLen = 0;
   tmp->dFragMap = NULL;
}

/***************  FreeData(Model B) **************************/
/*                                                           */
/*       free space in a model                               */
/*************************************************************/

void FreeData(Model B)
{
   int i, t, maxlen;

   maxlen = findMaxMotifLen(B->IP)+1;

   FreeRProb( B );    /* BT 12/17/97 */

   if( B->InitPos != NULL )
     FREEP( B->InitPos, B->IP->nNumMotifTypes );   /* BT 8/5/98 */

   FREEPP(B->C->dSumPseudo,  B->IP->nNumMotifTypes,2);
   FREEPP(B->First->dSumPseudo,  B->IP->nNumMotifTypes,2);
   FREEPP(B->C->dPseudoCounts, B->IP->nNumMotifTypes,maxlen);
   FREEPP(B->First->dPseudoCounts, B->IP->nNumMotifTypes,maxlen);
   FREEP(B->First->dmodel_sites, B->IP->nNumMotifTypes);
   FREEP(B->C->dmodel_sites, B->IP->nNumMotifTypes);
   FREEP(B->First->dmodel_pseudo, B->IP->nNumMotifTypes);
   FREEP(B->C->dmodel_pseudo, B->IP->nNumMotifTypes);
   free(B->First->dtot_sites);
   free(B->C->dtot_sites);
   free(B->First->dtot_pseudo);
   free(B->C->dtot_pseudo);
   free(B->First->dTot);
   free(B->C->dTot);
   free(B->First->dbg_pseudo);
   free(B->C->dbg_pseudo);
   FREEPP(B->First->fCounts, B->IP->nNumMotifTypes,maxlen);
   FREEPP(B->First->wCounts, B->IP->nNumMotifTypes,maxlen);
   FREEPP(B->C->fCounts, B->IP->nNumMotifTypes,maxlen);
   FREEPP(B->C->wCounts, B->IP->nNumMotifTypes,maxlen);

   if(B->IP->Datafiles->output_filename != NULL)
   {
      fclose( B->IP->Datafiles->out_fpt );		/* BT 1/27/97 */
      free(B->IP->Datafiles->output_filename);
   } 
   if(B->IP->Datafiles->prior_filename != NULL)
      free(B->IP->Datafiles->prior_filename);
   if(B->IP->Datafiles->WilcoxonFilename != NULL)     /* BT 7/23/97 */
     {
       remove( B->IP->Datafiles->WilcoxonFilename );
       free(B->IP->Datafiles->WilcoxonFilename);
     }
   if(B->IP->Datafiles->InputFilename != NULL)        /* BT 7/2/97 */
      free(B->IP->Datafiles->InputFilename);

   if(B->IP->Datafiles->ScanFilename != NULL)        /* BT 9/17/97 */
     {
       fclose( B->IP->Datafiles->sn_fpt );
       free(B->IP->Datafiles->ScanFilename);
     }
   if(B->IP->Datafiles->PriorOutFilename != NULL)        /* BT 9/26/97 */
     {
       fclose( B->IP->Datafiles->prout_fpt );
       free(B->IP->Datafiles->PriorOutFilename);
     }
   if(B->IP->Datafiles->ControlFileName != NULL)        /* BT 7/2/97 */
     free(B->IP->Datafiles->ControlFileName);

   if(B->IP->Datafiles->SampleFileName != NULL)        
     {
       if( B->IP->nRank == 0 )
	 fclose( B->IP->Datafiles->sample_fpt );
       else
	 remove( B->IP->Datafiles->SampleFileName );
	 
       free(B->IP->Datafiles->SampleFileName);
     }

   if(B->IP->Datafiles->OccurFileName != NULL)        /* BT 4/3/98 */
     {
       if( B->IP->nRank == 0 )
	 fclose( B->IP->Datafiles->occur_fpt );
       else
	 remove( B->IP->Datafiles->OccurFileName );
       free(B->IP->Datafiles->OccurFileName);
     }
   if(B->IP->Datafiles->NearOccurFileName != NULL)        /* BT 4/3/98 */
     {
       if( B->IP->nRank == 0 )
	 fclose( B->IP->Datafiles->near_occur_fpt );
       else
	 remove( B->IP->Datafiles->NearOccurFileName );
       free(B->IP->Datafiles->NearOccurFileName);
     }
   if(B->IP->Datafiles->SitesFileName != NULL)        /* BT 4/3/98 */
     {
       fclose( B->IP->Datafiles->sites_fpt );
       free(B->IP->Datafiles->SitesFileName);
     }
   if(B->IP->Datafiles->SitesPrFileName != NULL)        /* BT 4/3/98 */
     {
       if( B->IP->nRank == 0 )
	 fclose( B->IP->Datafiles->sitespr_fpt );
       else
	 remove( B->IP->Datafiles->SitesPrFileName );
       free(B->IP->Datafiles->SitesPrFileName);
     }
   if(B->IP->Datafiles->DistFileName != NULL)        /* BT 4/3/98 */
     {
       fclose( B->IP->Datafiles->dist_fpt );
       free(B->IP->Datafiles->DistFileName);
     }
   if(B->IP->Datafiles->ProbFileName != NULL)        
     {
       fclose( B->IP->Datafiles->prob_fpt );
       free(B->IP->Datafiles->ProbFileName);
     }
   if(B->IP->Datafiles->BkgndFileName != NULL)  
     {
       fclose( B->IP->Datafiles->bkgnd_fpt );
       free(B->IP->Datafiles->BkgndFileName);
     }
   if(B->IP->Datafiles->SpacingFileName != NULL)  
     {
       fclose( B->IP->Datafiles->spacing_fpt );
       free(B->IP->Datafiles->SpacingFileName);
     }
   if(B->IP->Datafiles->WeightFileName != NULL)  
     {
       fclose( B->IP->Datafiles->weight_fpt );
       free(B->IP->Datafiles->WeightFileName);
     }
   if(B->IP->Datafiles->mpiTempFileName != NULL)  
     {
       fclose( B->IP->Datafiles->mpiTemp_fpt );
       free(B->IP->Datafiles->mpiTempFileName);
     }
   if(B->IP->Datafiles->mpiSampleFileName != NULL)  
     {
       /* remove( B->IP->Datafiles->mpiSampleFileName );  */
       free(B->IP->Datafiles->mpiSampleFileName);
     }
   if( B->IP->Datafiles->hierTempFileName != NULL )
     {
       fclose( B->IP->Datafiles->hier_fpt );
       remove( B->IP->Datafiles->hierTempFileName );
       free(B->IP->Datafiles->hierTempFileName);
     }
   free(B->IP->Datafiles->filename);

   free(B->IP->Datafiles);
   free(B->IP->dposterior_prob);
   free(B->IP->nMotifLen);
   free(B->IP->nMinMotifLen);    /* BT 5/9/97 */
   free(B->IP->nMaxMotifLen);    /* BT 5/9/97 */
   free(B->IP->nInputMotifLen);    /* BT 5/23/97 */
   free(B->IP->DOF);
   free(B->IP->nPossSites);
   free( B->IP->nMaxFragWidth );
   free( B->IP->col_shift );

   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       free(B->IP->nNumMotifs[t]);
       if(B->IP->AltModel->Collapsed != NULL)
         free(B->IP->AltModel->Collapsed[t]);
       if(B->IP->AltModel->Palandromic != NULL)
         free(B->IP->AltModel->Palandromic[t]);
       if(B->IP->AltModel->Repeat != NULL)
         free(B->IP->AltModel->Repeat[t]);
       if(B->IP->AltModel->Concen != NULL)
         free(B->IP->AltModel->Concen->Concentrated[t]);
       if(B->F != NULL)  
	 {
	   free(B->F->nColMask[t]);
	   free(B->F->nOldColMask[t]); 
	   free(B->F->fragPos[t]); 
	   free(B->F->fragInitMask[t]); 
	 } 
     }
   
   if( B->IP->nWidthCnts != NULL )
     FREEPP( B->IP->nWidthCnts, B->IP->nNumMotifTypes, 2 );   /* BT 2/18/98 */

   if( B->IP->nAlignCnts != NULL )         /* BT 3/27/98 */
     {
       for( i = 0; i <= B->IP->nSeeds; i++ )
	 FREEPP( B->IP->nAlignCnts[i],  B->IP->nNumMotifTypes,  B->IP->nNumSequences );
       free( B->IP->nAlignCnts );
     }

   if( B->IP->is_defined[cl_B] )
     {
       for( i = 0; i < B->IP->nNumSequences; i++ )
	 {
	   FREEP( B->BP->dBkgndProb[i], SequenceLength( B, i ) );
	 }
       free( B->BP->dBkgndProb );
       free( B->BP );
     }
      
   if( B->IP->is_defined[cl_X] )
     {
       free( B->AN->dTemp );

       for( i = 0; i < B->AN->nTemps; i++ )
	 free_maxdata( &(B->AN->results[i]), B->IP);

       free( B->AN->results );
       free( B->AN );
     }

   if( B->WT != NULL )
     {
       if(  B->WT->weightTable != NULL )
	 {
	   for( i = 0; i < B->WT->nTableSize; i++ )
	     {
	       free( B->WT->weightTable[i].weight );
	     }
	   free( B->WT->weightTable );
	 }
       free( B->WT->seqCodes );
       free( B->WT->weight );
       free( B->WT );
     }

   if(B->IP->AltModel->Collapsed != NULL)
      free(B->IP->AltModel->Collapsed);
   if(B->IP->AltModel->Palandromic != NULL)
      free(B->IP->AltModel->Palandromic);
   if(B->IP->AltModel->Repeat != NULL)
      free(B->IP->AltModel->Repeat);
   if(B->IP->AltModel->Concen != NULL) {
      free(B->IP->AltModel->Concen->Concentrated);
      free(B->IP->AltModel->Concen->NumConcen);
      free(B->IP->AltModel->Concen);
   }
   free(B->IP->AltModel);    /* BT 8/27/97 */

   if(B->F != NULL) {
      free(B->F->nColMask);
      free(B->F->nOldColMask);
      free(B->F->fragPos);
      free(B->F->fragInitMask);
      free(B->F->nMaxLen);
      free(B->F->FragWidth);
      free(B->F->nOldFragWidth);
      free(B->F->shift);
      free(B->F); 
   } 

   if( *B->Seq->ProcessedSTR != *B->Seq->R )       /* BT 10/20/97 */
     free(*B->Seq->ProcessedSTR);
   free( B->Seq->ProcessedSTR );
   
   if( B->IP->comment != NULL )
     {
       free( B->IP->comment );
     }

   if( B->Phylo->treeCount )
     {
       for( i = 0; i < B->Phylo->treeCount; i++ )
	 FreePhylTree( B, B->Phylo->phyloTree[i] );
       free( B->Phylo->phyloTree );
     }
   if( B->Phylo->phyloSeq )
     free( B->Phylo->phyloSeq );
   if( B->Phylo->phyloSpecies )
     free( B->Phylo->phyloSpecies );
   if( B->Phylo->phyloSpeciesSample )
     free( B->Phylo->phyloSpeciesSample );
   if( B->Phylo->phyloTreeStartSeq )
     free( B->Phylo->phyloTreeStartSeq );
   if( B->Phylo->phyloIndex )
     free( B->Phylo->phyloIndex );
   if( B->Phylo->phyloBkgnd )
     FREEP( B->Phylo->phyloBkgnd, B->Phylo->phyloSpeciesCount );
   free( B->Phylo );

   if( B->IP->programName )
     free( B->IP->programName );

   for( i = 0; i < B->IP->nNumSequences; i++ )
     free( B->IP->fastaHeader[i] );
   free( B->IP->fastaHeader );

   if( B->IP->freqSeqs )
     {
       FreeNMerArrays( B );
       free( B->IP->freqSeqs );
       free( B->IP->nMerSize );
     }

   free(*B->Seq->R);
   free(B->IP->nNumMotifs);
   free(*B->Seq->nvEndLocs);    /* BT 8/27/97 */
   free(B->Seq->nvEndLocs);    /* BT 8/27/97 */
   free(B->Seq->R);
   /* free individual sequences */
   /*   for(i=0;i<=B->IP->nNumSequences-1;i++) */
   for(i=0; i < B->IP->nSeqAllocCnt; i++ )    /* BT 8/27/97 */
     free(B->Seq->Orig[i]);
   free(B->Seq->Orig);
   free(B->Seq->SeqLen);
   free(B->First);
   free(B->C); 
   free(B->Seq);
   free(B->IP);
   FreeAlignmentCounts( B );      /* BT 1/16/98 */

   free(B);
   B = NULL;
}



/************************** free_maxdata  ************************/

void free_maxdata(MaxResults *tmp, IPtype IP)
 
{
   int maxlen=tmp->nMaxLen;

   if(tmp->F != NULL) 
     {
       if( tmp->F->nMaxLen != NULL )
	 free(tmp->F->nMaxLen);
       if( tmp->F->FragWidth != NULL )
	 free(tmp->F->FragWidth);
       if( tmp->F->nOldFragWidth != NULL )
	 free(tmp->F->nOldFragWidth);
       if( tmp->F->shift != NULL )
	 free(tmp->F->shift);
       if(tmp->F->nColMask != NULL)
         FREEP(tmp->F->nColMask, IP->nNumMotifTypes);
       if(tmp->F->nOldColMask != NULL)
         FREEP(tmp->F->nOldColMask, IP->nNumMotifTypes);
       if(tmp->F->fragPos != NULL)
         FREEP(tmp->F->fragPos, IP->nNumMotifTypes);
       free(tmp->F);
     }
   
   if(tmp->nNumMotifs != NULL) 
     free(tmp->nNumMotifs);
   
   if(tmp->nMotifLoc != NULL)
      FREEP(tmp->nMotifLoc, maxlen);
   if(tmp->RevComp != NULL)
      FREEP(tmp->RevComp, maxlen);
   if(tmp->dvMotifProb != NULL)
      FREEP(tmp->dvMotifProb, maxlen);
   if(tmp->frequency != NULL)
      FREEP(tmp->frequency, IP->nSeqLen);
   if( tmp->nMotifLen != NULL )           /* BT 5/23/97 */
     free( tmp->nMotifLen );
   if( tmp->dMap != NULL )           /* BT 3/4/98 */
     free( tmp->dMap );
   if( tmp->dFragMap != NULL )           /* BT 3/4/98 */
     free( tmp->dFragMap );
}











