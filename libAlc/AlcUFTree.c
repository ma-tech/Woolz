/* TODO Written but needs to be debuged tested */

#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlcUFTree_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         AlcUFTree.c
* \author       Bill Hill
* \date         September 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief	A general purpose union tree based on Sedgewick's
* 		Weighted Quick Union Find, see
* 		Robert Sedgewick, Kevin Wayne "Algorithms (4th Edition)".
* 		There is very little parameter checking in these functions
* 		and all given parameters must be valid.
* \ingroup	AlcUFTree
*/

#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/*!
* \ingroup	AlcUFTree
* \brief	Free a union find tree data structure.
* \param	uft			Given union find tree.
*/
void				AlcUFTreeFree(
				  AlcUFTree *uft)
{
  if(uft)
  {
    AlcFree(uft->pr);	      /* uft->sz allocated with uft->pr then offset. */
    AlcFree(uft);
  }
}

/*!
* \ingroup	AlcUFTree
* \brief	Initialises a union find tree data structure to allow reuse.
* \param	uft			The union find tree.
* \param	nNod			Number of nodes, which must be less
* 					than the value of maxNod in the data
* 					structure..
*/
void				AlcUFTreeInit(
				  AlcUFTree *uft,
				  int n)
{
  int		i;

  uft->nNod = n;
  uft->nCmp = n;
  for(i = 0; i < n; ++i)
  {
    uft->pr[i] = i;
    uft->sz[i] = 1;
  }
}

/*!
* \return	Union find tree, or NULL on error.
* \ingroup	AlcUFTree
* \brief	Creates a union find tree data structure which is required
*		by all the other AlcUFTree functions.
* \param	maxNod 			Maximum number of nodes space
* 					allocated for.
* \param	nNod			Number of nodes.
*/
AlcUFTree			*AlcUFTreeNew(
  				  int maxNod,
				  int nNod)
{
  AlcUFTree	*uft = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if(((uft = (AlcUFTree *)AlcCalloc(1, sizeof(AlcUFTree))) == NULL) ||
     ((uft->pr = (int *)AlcMalloc(maxNod * 2 * sizeof(int))) == NULL))
  {
    AlcFree(uft);
    uft = NULL;
  }
  else
  {
    uft->sz = uft->pr + maxNod;         /* Add offset to uft->pr for uft->pr */
    uft->maxNod = maxNod;
    AlcUFTreeInit(uft, nNod);
  }
  return(uft);
}


/*!
* \return	The component containing the given node.
* \ingroup	AlcUFTree
* \brief	Finds the component containing the given node.
* \param	uft			The union find tree.
* \param	maxNod 			Maximum number of nodes space
* 					allocated for.
* \param	nNod			Number of nodes.
*/
int				AlcUFTreeFind(
				  AlcUFTree *uft,
				  int p)
{
  int		r,
  		s;

  r = p;
  while(r != (s = uft->pr[uft->pr[r]]))
  {
    r = s;
  }
  return(r);
}

/*!
* \return	Non-zero if the nodes are connected.
* \ingroup	AlcUFTree
* \brief	Returns non-zero if the two given nodes are connected,
* 		ie if they are in the same component.
* \param	uft			The union find tree.
* \param	p			Node in first component.
* \param	q			Node in second component.
*/
int				AlcUFTreeConnected(
				  AlcUFTree *uft,
				  int p,
				  int q)
{
  int		con;

  con = (AlcUFTreeFind(uft, p) == AlcUFTreeFind(uft, q));
  return(con);
}
  
/*!
* \ingroup	AlcUFTree
* \brief	If the two given nodes have different components then
* 		their components are merged to form one.
* \param	uft			The union find tree.
* \param	p			Node in first component.
* \param	q			Node in second component.
*/
void				AlcUFTreeUnion(
				  AlcUFTree *uft,
				  int p,
				  int q)
{
  int		rP,
  		rQ;

  rP = AlcUFTreeFind(uft, p);
  rQ = AlcUFTreeFind(uft, q);
  if(rP != rQ)
  {
    if(uft->sz[rP] < uft->sz[rQ])
    {
      uft->pr[rP] = rQ;
      uft->sz[rQ] += uft->sz[rP];
    }
    else
    {
      uft->pr[rQ] = rP;
      uft->sz[rP] += uft->sz[rQ];
    }
    --(uft->nCmp);
  }
}

#ifdef ALC_UFTREE_MAIN
int		main(int argc, char *argv[])
{
  int		i,
  		n = 10;
  AlcUFTree 	*uft;

  uft = AlcUFTreeNew(10,10);
  for(i = 0; i + 2 < n; i += 4)
  {
    AlcUFTreeUnion(uft, i, i + 2);
  }
  for(i = 1; i + 2 < n; i += 4)
  {
    AlcUFTreeUnion(uft, i, i + 2);
  }
  AlcUFTreeUnion(uft, 1, 5);
  AlcUFTreeUnion(uft, 2, 4);
  AlcUFTreeUnion(uft, 8, 9);
  
  (void )printf("nNod   %d\n", uft->nNod);
  (void )printf("nCmp   %d\n", uft->nCmp);
  (void )printf("maxNod %d\n", uft->maxNod);
  (void )printf("nodes:       ");
  for(i = 0; i < n; ++i)
  {
    (void )printf(" % 2d", i);
  }
  (void )printf("\n");
  (void )printf("node parents:");
  for(i = 0; i < n; ++i)
  {
    (void )printf(" % 2d", uft->pr[i]);
  }
  (void )printf("\n");
  (void )printf("Comp sizes:  ");
  for(i = 0; i < n; ++i)
  {
    (void )printf(" % 2d", uft->sz[i]);
  }
  (void )printf("\n");
  AlcUFTreeFree(uft);
}
#endif /* ALC_UFTREE_MAIN */
