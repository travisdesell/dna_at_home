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
#include "phyllip.h"

PhyloTree AllocPhyloNode( Model B );
PhyloTree PhylTree( Model B );
PhyloTree PhylExpr( Model B );
PhyloTree PhylSubtree( Model B );
void      PhyloScan();
void      GetToken();
char      GetChar();
void      PutChar();
void      CalcSubsMatrix2( Model B, PhyloTree tree, double k, double *prob );
void      DumpPhylTree2( PhyloTree tree, int indent );

#define GAMMA  1.0
#define ALPHA  (1.0/3.0) 

/* some stuff local to this file to make life easier for the parser
   we don't have to pass so many pointers around.
   
   store the current token and the string here 
   strPos keeps track of where the parser is in the string 

   speciesCnt keeps track of the species offset
   sequences in the fasta file must be in the same order as the tree
*/

static struct PHYLOTOKEN
{
  char str[MAX_PHY_STRING];
  int  type;
} token;

static char *treeStr;
static int  strPos;
static int  speciesCnt;


PhyloTree AllocPhyloNode( Model B )
{
  PhyloTree  node;
  int        i;

  NEW( node, 1, PhyloNode );
  NEW( node->likelihood, B->IP->nAlphaLen, double );
  NEWP( node->subs, B->IP->nAlphaLen, double );
  for( i = 0; i < B->IP->nAlphaLen; i++ )
    NEW( node->subs[i], B->IP->nAlphaLen, double );

  return( node );
}


int ReadPhylTree( Model B, char *phyloTreeStr, PhyloTree *tree )
{
  if( ! B->IP->is_defined[cl_E] )
    p_error( "ReadPhylTree: Phylogenetic trees can only be used with the -E option." );

  treeStr = phyloTreeStr;
  strPos = 0;
  speciesCnt = 0;
  
  PhyloScan();
  *tree = PhylTree( B );

#ifdef _DEBUG_
  DumpPhylTree( *tree );
#endif

  return( speciesCnt );
}


/* Parse the tree
   Assumes the grammar of the Phyllip tree is
   phyl_tree ==> descendant_list [ : branch_length ] ;

   descendant_list ==> ( subtree { , subtree } )

   subtree ==> descendant_list [: branch_length]
           ==> leaf_label [: branch_length]
  
   phyl_tree := phyl_expr;
   phyl_expr := (subtree {, subtree})
   subtree   := species:len |
                phyl_expr {:len} 
     
   I'm not sure this is correct but it conforms to the example
*/

/* Simple recursive descent parser without much error checking
   returns a pointer to the root node of the tree
*/
PhyloTree PhylTree( Model B )
{
  PhyloTree tree;
  
  tree = PhylExpr( B );
  
  if( strcmp( token.str,  ";" ) != 0  )
    {
      p_error( "PhylTree: Missing ','" );
    }
 
  return( tree );
}


PhyloTree PhylExpr( Model B )
{
  PhyloTree tree = NULL;
  PhyloTree node1 = NULL;
  PhyloTree node2 = NULL;

  if( strcmp( token.str, "(" ) == 0 )
    {
      PhyloScan();
      tree = PhylSubtree( B );
      while( strcmp( token.str,  "," ) == 0  )
	{
	  node1 = tree;
	  PhyloScan();
	  node2 = PhylSubtree( B );	  
	  tree = AllocPhyloNode( B );
	  tree->leftChild = node1;
	  tree->rightChild = node2;
	}

      if( strcmp( token.str,  ")" ) != 0  )
	{
	  p_error( "PhylTree: Missing ')'" );
	}
      PhyloScan();
    }
  else
    {
      p_error( "PhylTree: Missing '('" );
    }

  return( tree );
  }


PhyloTree PhylSubtree( Model B )
{
  PhyloTree tree = NULL;
  char      *str;

  if( strcmp( token.str, "(" ) == 0 )
    {
      tree = PhylExpr( B );
      if( strcmp( token.str, ":" ) == 0 )
	{
	  PhyloScan();
	  tree->length = strtod( token.str, NULL );
	  PhyloScan();	  
	}
    }
  else
    {
      tree = AllocPhyloNode( B );
      
      str = strtok( token.str, ":" );
      if( ! str )
	p_error( "PhylExpr: invalid species name in phylogenetic tree." );

      NEW( tree->species, strlen( str ) + 1, char );
      strcpy( tree->species, str );
      str = strtok( NULL, ":" );
      if( ! str )
	p_error( "PhylExpr: invalid species name in phylogenetic tree." );

      tree->length = strtod( str, NULL );
      tree->offset = speciesCnt;
      speciesCnt++;
      PhyloScan();
    }

  return( tree );
}


/* Get a token from treeStr */
/* and an infinite supply of )'s at end of string  */
void PhyloScan()
{
  GetToken();
  if( token.type == PHY_EOF )
    strcpy( token.str, ")" );
}


/* Set the package variable token with the current token from treestr
   token contains the token string and the token type
   tokens are: ( ) string:number : or number 
*/
void GetToken()
{
  char  ch;
  int   chCount;
  
  token.str[0] = '\0';

  ch = GetChar();
      
  while( ch != PHY_EOF && isspace( (int) ch ) )
      {
	ch = GetChar();
      }
  
  if( ch == PHY_EOF )
    {
      token.type = PHY_EOF;
    }
  else if( ch == ':' )
    {
      token.str[0] = ch;
      token.str[1] = '\0';
      token.type = PHY_LENGTH;
    }
  else if( ch == ';' )
    {
      token.str[0] = ch;
      token.str[1] = '\0';
      token.type = PHY_TERMINATOR;
    }
  else if( isdigit( (int) ch ) )
    {
      chCount = 0;
      token.type = PHY_NUMBER;
      token.str[chCount] = ch;
      chCount++;
      ch = GetChar();
      while( isdigit( (int) ch ) || ch == '.'  )
	{
	  token.str[chCount] = ch;
	  chCount++;
	  ch = GetChar();
	}
      token.str[chCount] = '\0';      
      PutChar( ch );
    }
  else if( isalpha( (int) ch ) )
    {
      chCount = 0;
      token.type = PHY_STRING;
      token.str[chCount] = ch;
      chCount++;
      ch = GetChar();
      while( isalnum( (int) ch ) || ch == ':' || isdigit( (int) ch ) || 
	     ch == '.' || ch == '_' || ch == '-' || ch == '|' )
	{
	  token.str[chCount] = ch;
	  chCount++;
	  if( chCount >= MAX_PHY_STRING - 1 )
	    p_error( "GetToken: species name exceeds the maximum length." );
	  ch = GetChar();
	}
      token.str[chCount] = '\0'; 
      PutChar( ch );
    }
  else
    {
      token.str[0] = ch;
      token.str[1] = '\0';
      token.type = PHY_DELIMITER;	
    }

#ifdef _DEBUG_
  printf( "token = %s\n", token.str );
  fflush(stdout);
#endif
}


/* pull a character from treeStr */
char GetChar()
  {
    char ch;

    if( strPos >= strlen( treeStr ) )
      {
	ch = PHY_EOF;
      }
    else
      {
	ch = treeStr[strPos];
	strPos++;
      }

    return( ch );
  }


void PutChar()
{
  if( strPos > 0 )
    strPos--;
}


void CalcSubsMatrix( Model B, PhyloTree tree, double *prob )
{
  int    i;
  double k;

  for( k = 1.0, i = 0; i < B->IP->nAlphaLen; i++ )
    k -= prob[i] * prob[i];
  k = 1.0 / k;
  
  if( k < 0.0 )
    k = fabs( k );

  CalcSubsMatrix2( B, tree, k, prob );
}


void CalcSubsMatrix2( Model B, PhyloTree tree, double k, double *prob )
{
  int    i;
  int    j;
  double d;
  double d2;
  double *trSubs; 

  if( tree )
    {
      CalcSubsMatrix2( B, tree->leftChild, k, prob );
      CalcSubsMatrix2( B, tree->rightChild, k, prob );

      if( k > MAX_K )
	d = 0;
      else
	d = exp( -k * tree->length );

      d2 = 1 - d;

      for( i = 0; i < B->IP->nAlphaLen; i++ )
	{
	  trSubs = tree->subs[i];
	  for( j = 0; j < B->IP->nAlphaLen; j++ )
	    {
	      if( i == j )
		trSubs[j] = d + d2 * prob[j];
	      else
		trSubs[j] = d2 * prob[j]; 
	    }
	} 
    }
}


void FreePhylTree( Model B, PhyloTree tree )
{
  if( tree )
    {
      FreePhylTree( B, tree->leftChild );
      FreePhylTree( B, tree->rightChild );
      if( tree->species)
	free( tree->species );
      if( tree->likelihood )
	free( tree->likelihood );
      if( tree->subs )
	FREEP( tree->subs, B->IP->nAlphaLen );
      free( tree );
    }
}


void DumpPhylTree( PhyloTree tree )
{
  DumpPhylTree2( tree, 0 );
  fflush( stdout );
}


void DumpPhylTree2( PhyloTree tree, int indent )
{
  int  i;

  if( tree )
    {
      for( i = 0; i < indent; i++ )
	printf( " " );
      printf( "%lx ", (unsigned long int) tree );
      if( tree->species )
	printf( "%s ", tree->species );
      else
	printf( "- " );
      printf( "%6.4f\n", tree->length );
      DumpPhylTree2( tree->leftChild, indent + 2 );
      DumpPhylTree2( tree->rightChild, indent + 2 );
    }
}


