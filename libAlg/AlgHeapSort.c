#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgHeapSort_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/AlgHeapSort.c
* \author       Bill Hill
* \date         July 2000
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief        General purpose heap sort algorithms.
* \ingroup      AlgSort
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Alg.h>

static void 			AlgHeapSiftDown(
			  	  char *data,
				  int width,
				  int (*cmpFn)(void *, void *),
				  int idL,
				  int idU);
static void 			AlgHeapSiftDownIdx(
				  void *data,
				  int *idx,
				  int (*cmpFn)(void *, int *, int, int),
				  int idL,
				  int idU);

/*!
* \return	Error code.
* \ingroup	AlgSort
* \brief	Sorts given data using a heapsort algorithm.
* \param	data			The given data to sort.
* \param	nElm			Number of data elements.
* \param	elmSz			Size of the data elements.
* \param	cmpFn			Given element comparison
*					function which must return
*					an integer < 0, == 0 or > 0 to
*					indicate that the first entry
*					is <, == or > the second.
*/
AlgError	AlgHeapSort(void *data, unsigned nElm, unsigned elmSz,
			     int (*cmpFn)(void *, void *))
{
  int 		idI,
  		idL;
  char		*dataC;
  AlgError	algErr = ALG_ERR_NONE;

  if((data == NULL) || (cmpFn == NULL) || (elmSz < 1))
  {
    algErr = ALG_ERR_FUNC;
  }
  else if(nElm > 1)
  {
    idL = (nElm - 1) * elmSz;
    idI = ((idL / 2) + 1) * elmSz;
    dataC = (char *)data;
    /* Build the heap. */
    while(idI > 0)
    {
      idI -= elmSz;
      AlgHeapSiftDown(dataC, elmSz, cmpFn, idI, idL);
    }
    idI = idL;
    /* Repeatedly swap the top value and then sift down the lower
     * value. */
    while(idI > 0)
    {
      AlgHeapElmSwap(dataC + idI, dataC, elmSz);
      idI -= elmSz;
      AlgHeapSiftDown(dataC, elmSz, cmpFn, 0, idI);
    }
  }
  return(algErr);
}

/*!
* \return	Error code.
* \ingroup	AlgSort
* \brief	Sorts given data using a heapsort algorithm only the
*		indicies are modified.
* \param	data			The given data.
* \param	idx			Data indicies to sort.
* \param	nElm			Number of data elements.
* \param	cmpFn			Given indexed
*					element comparison function
*					which must return an integer
*					that is < 0, == 0 or > 0 to
*					indicate that the first entry
*					is <, == or > the second indexed
*					entry.
*/
AlgError	AlgHeapSortIdx(void *data, int *idx,
			       unsigned nElm,
			       int (*cmpFn)(void *, int *, int, int))
{
  int		idI,
  		idL,
		swapI;
  AlgError	algErr = ALG_ERR_NONE;

  if((data == NULL) || (idx == NULL) || (cmpFn == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else if(nElm > 1)
  {
    idL = nElm - 1;
    idI = (idL / 2) + 1;
    /* Build the heap. */
    while(idI > 0)
    {
      --idI;
      AlgHeapSiftDownIdx(data, idx, cmpFn, idI, idL);
    }
    idI = idL;
    /* Repeatedly swap the top value's index and then sift down the lower
     * value's index. */
    while(idI > 0)
    {
      swapI = *(idx + idI); *(idx + idI) = *(idx + 0); *(idx + 0) = swapI;
      --idI;
      AlgHeapSiftDownIdx(data, idx, cmpFn, 0, idI);
    }
  }
  return(algErr);
}

/*!
* \return	void
* \ingroup	AlgSort
* \brief	Swaps two data elements.
* \param	elm0			Ptr to first element.
* \param	elm1			Ptr to second element.
* \param	cnt			Element size.
*/
void 		AlgHeapElmSwap(void *elm0, void *elm1, int cnt)
{
  char		tC0;
  char		*ePC0,
  		*ePC1;

  ePC0 = (char *)elm0;
  ePC1 = (char *)elm1;
  while(cnt-- > 0)
  {
    tC0 = *ePC1; *ePC1++ = *ePC0; *ePC0++ = tC0;
  }
}

/*!
* \return	void
* \ingroup	AlgSort
* \brief	Sifts a data element down through the heap.
* \param	data			Data elements.
* \param	elmSz			Size of a data element.
* \param	cmpFn			Given element comparison function.
* \param	idL			Lower data index.
* \param	idu			Upper data index.
*/
static void 	AlgHeapSiftDown(char *data, int elmSz,
			        int (*cmpFn)(void *, void *),
				int idL, int idU)
{
  int 		idI;

  idI = idL * 2;
  while(idI <= idU)
  {
    if((idI < idU) && ((*cmpFn)(data + idI + elmSz, data + idI) > 0))
    {
      idI += elmSz;
    }
    if((*cmpFn)(data + idL, data + idI) >= 0)
    {
      idI = idU + 1;
    }
    else
    {
      AlgHeapElmSwap(data + idI, data + idL, elmSz);
      idL = idI;
      idI = idL * 2;
    }
  }
}

/*!
* \return	void
* \ingroup	AlgSort
* \brief	Sifts an indexed data element down through the heap.
* \param	data			Data elements.
* \param	idx			Indicies to data elements.
* \param	cmpFn			Given element comparison function.
* \param	idL			Lower data index.
* \param	idu			Upper data index.
*/
static void 	AlgHeapSiftDownIdx(void *data, int *idx,
				   int (*cmpFn)(void *, int *, int, int),
				   int idL, int idU)
{
  int		idI,
  		swap;

  idI = idL * 2;
  while(idI <= idU)
  {
    if((idI < idU) && ((*cmpFn)(data, idx, idI + 1, idI) > 0))
    {
      ++idI;
    }
    if((*cmpFn)(data, idx, idL, idI) >= 0)
    {
      idI = idU + 1;
    }
    else
    {
      swap = *(idx + idI); *(idx + idI) = *(idx + idL); *(idx + idL) = swap;
      idL = idI;
      idI = idL * 2;
    }
  }
}


/*!
* \return	*((char *)datum0) - *((char *)datum1)
* \ingroup	AlgSort
* \brief	Char comparison function for AlgHeapSort().
*		Sorted data will have smallest entry first and greatest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortCmpCFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((char *)datum0) - *((char *)datum1);
  return(cmp);
}

/*!
* \return	*((unsigned char *)datum0) - *((unsigned char *)datum1)
* \ingroup	AlgSort
* \brief	Unsigned char comparison function for AlgHeapSort().
*		Sorted data will have smallest entry first and greatest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortCmpUFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((unsigned char *)datum0) - *((unsigned char *)datum1);
  return(cmp);
}

/*!
* \return	*((short *)datum0) - *((short *)datum1)
* \ingroup	AlgSort
* \brief	Short comparison function for AlgHeapSort().
*		Sorted data will have smallest entry first and greatest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortCmpSFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((short *)datum0) - *((short *)datum1);
  return(cmp);
}

/*!
* \return	*((int *)datum0) - *((int *)datum1)
* \ingroup	AlgSort
* \brief	Int comparison function for AlgHeapSort().
*		Sorted data will have smallest entry first and greatest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortCmpIFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((int *)datum0) - *((int *)datum1);
  return(cmp);
}

/*!
* \return	*((long *)datum0) - *((long *)datum1)
* \ingroup	AlgSort
* \brief	Long comparison function for AlgHeapSort().
*		Sorted data will have smallest entry first and greatest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortCmpLFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((long *)datum0) - *((long *)datum1);
  return(cmp);
}

/*!
* \return	*((float *)datum0) - *((float *)datum1)
* \ingroup	AlgSort
* \brief	Float comparison function for AlgHeapSort().
*		Sorted data will have smallest entry first and greatest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortCmpFFn(void *datum0, void *datum1)
{
  int		cmp = 0;
  float		cmpF;

  cmpF = *((float *)datum0) - *((float *)datum1);
  if(cmpF < 0.0)
  {
    cmp = -1;
  }
  else if(cmpF > 0.0)
  {
    cmp = 1;
  }
  return(cmp);
}

/*!
* \return	*((double *)datum0) - *((double *)datum1)
* \ingroup	AlgSort
* \brief	Double comparison function for AlgHeapSort().
*		Sorted data will have smallest entry first and greatest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortCmpDFn(void *datum0, void *datum1)
{
  int		cmp = 0;
  double	cmpD;

  cmpD = *((double *)datum0) - *((double *)datum1);
  if(cmpD < 0.0)
  {
    cmp = -1;
  }
  else if(cmpD > 0.0)
  {
    cmp = 1;
  }
  return(cmp);
}

/*!
* \return	*((char *)datum1) - *((char *)datum0)
* \ingroup	AlgSort
* \brief	Inverse char comparison function for AlgHeapSort().
*		Sorted data will have greatest entry first and smallest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortInvCmpCFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((char *)datum1) - *((char *)datum0);
  return(cmp);
}

/*!
* \return	*((unsigned char *)datum1) - *((unsigned char *)datum0)
* \ingroup	AlgSort
* \brief	Inverse unsigned char comparison function for AlgHeapSort().
*		Sorted data will have greatest entry first and smallest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortInvCmpUFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((unsigned char *)datum1) - *((unsigned char *)datum0);
  return(cmp);
}

/*!
* \return	*((short *)datum1) - *((short *)datum0)
* \ingroup	AlgSort
* \brief	Inverse short comparison function for AlgHeapSort().
*		Sorted data will have greatest entry first and smallest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortInvCmpSFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((short *)datum1) - *((short *)datum0);
  return(cmp);
}

/*!
* \return	*((int *)datum1) - *((int *)datum0)
* \ingroup	AlgSort
* \brief	Inverse int comparison function for AlgHeapSort().
*		Sorted data will have greatest entry first and smallest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortInvCmpIFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((int *)datum1) - *((int *)datum0);
  return(cmp);
}

/*!
* \return	*((long *)datum1) - *((long *)datum0)
* \ingroup	AlgSort
* \brief	Inverse long comparison function for AlgHeapSort().
*		Sorted data will have greatest entry first and smallest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortInvCmpLFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = *((long *)datum1) - *((long *)datum0);
  return(cmp);
}

/*!
* \return	*((float *)datum1) - *((float *)datum0)
* \ingroup	AlgSort
* \brief	Inverse float comparison function for AlgHeapSort().
*		Sorted data will have greatest entry first and smallest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortInvCmpFFn(void *datum0, void *datum1)
{
  int		cmp = 0;
  float		cmpF;

  cmpF = *((float *)datum1) - *((float *)datum0);
  if(cmpF < 0.0)
  {
    cmp = -1;
  }
  else if(cmpF > 0.0)
  {
    cmp = 1;
  }
  return(cmp);
}

/*!
* \return	*((double *)datum1) - *((double *)datum0)
* \ingroup	AlgSort
* \brief	Inverse double comparison function for AlgHeapSort().
*		Sorted data will have greatest entry first and smallest last.
* \param	datum0			First data pointer.
* \param	datum1			Second data pointer.
*/
int		AlgHeapSortInvCmpDFn(void *datum0, void *datum1)
{
  int		cmp = 0;
  double	cmpD;

  cmpD = *((double *)datum1) - *((double *)datum0);
  if(cmpD < 0.0)
  {
    cmp = -1;
  }
  else if(cmpD > 0.0)
  {
    cmp = 1;
  }
  return(cmp);
}

/*!
* \return	*((char *)data + *(idx + id0)) -
*		*((char *)data + *(idx + id1))
* \ingroup	AlgSort
* \brief	Char comparison function for AlgHeapSortIdx().
*		Sorted data will have smallest entry first and greatest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortCmpIdxCFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((char *)data + *(idx + id0)) - *((char *)data + *(idx + id1));
  return(cmp);
}

/*!
* \return	*((unsigned char *)data + *(idx + id0)) -
*		*((unsigned char *)data + *(idx + id1))
* \ingroup	AlgSort
* \brief	Unsigned char comparison function for AlgHeapSortIdx().
*		Sorted data will have smallest entry first and greatest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortCmpIdxUFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((unsigned char *)data + *(idx + id0)) -
        *((unsigned char *)data + *(idx + id1));
  return(cmp);
}

/*!
* \return	*((short *)data + *(idx + id0)) -
*		*((short *)data + *(idx + id1))
* \ingroup	AlgSort
* \brief	Short comparison function for AlgHeapSortIdx().
*		Sorted data will have smallest entry first and greatest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortCmpIdxSFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((short *)data + *(idx + id0)) - *((short *)data + *(idx + id1));
  return(cmp);
}

/*!
* \return	*((int *)data + *(idx + id0)) -
*		*((int *)data + *(idx + id1))
* \ingroup	AlgSort
* \brief	Int comparison function for AlgHeapSortIdx().
*		Sorted data will have smallest entry first and greatest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortCmpIdxIFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((int *)data + *(idx + id0)) - *((int *)data + *(idx + id1));
  return(cmp);
}

/*!
* \return	*((long *)data + *(idx + id0)) -
*		*((long *)data + *(idx + id1))
* \ingroup	AlgSort
* \brief	Long comparison function for AlgHeapSortIdx().
*		Sorted data will have smallest entry first and greatest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortCmpIdxLFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((long *)data + *(idx + id0)) - *((long *)data + *(idx + id1));
  return(cmp);
}

/*!
* \return	*((float *)data + *(idx + id0)) -
*		*((float *)data + *(idx + id1))
* \ingroup	AlgSort
* \brief	Float comparison function for AlgHeapSortIdx().
*		Sorted data will have smallest entry first and greatest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortCmpIdxFFn(void *data, int *idx, int id0, int id1)
{
  int		cmp = 0;
  float		cmpF;

  cmpF = *((float *)data + *(idx + id0)) - *((float *)data + *(idx + id1));
  if(cmpF < 0.0)
  {
    cmp = -1;
  }
  else if(cmpF > 0.0)
  {
    cmp = 1;
  }
  return(cmp);
}

/*!
* \return	*((double *)data + *(idx + id0)) -
*		*((double *)data + *(idx + id1))
* \ingroup	AlgSort
* \brief	Double comparison function for AlgHeapSortIdx().
*		Sorted data will have smallest entry first and greatest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortCmpIdxDFn(void *data, int *idx, int id0, int id1)
{
  int		cmp = 0;
  double	cmpD;

  cmpD = *((double *)data + *(idx + id0)) - *((double *)data + *(idx + id1));
  if(cmpD < 0.0)
  {
    cmp = -1;
  }
  else if(cmpD > 0.0)
  {
    cmp = 1;
  }
  return(cmp);
}

/*!
* \return	*((char *)data + *(idx + id1)) -
*		*((char *)data + *(idx + id0))
* \ingroup	AlgSort
* \brief	Inverse char comparison function for AlgHeapSortIdx().
*		Sorted data will have greatest entry first and smallest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortInvCmpIdxCFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((char *)data + *(idx + id1)) - *((char *)data + *(idx + id0));
  return(cmp);
}

/*!
* \return	*((unsigned char *)data + *(idx + id1)) -
*		*((unsigned char *)data + *(idx + id0))
* \ingroup	AlgSort
* \brief	Inverse unsigned char comparison function for AlgHeapSortIdx().
*		Sorted data will have greatest entry first and smallest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortInvCmpIdxUFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((unsigned char *)data + *(idx + id1)) -
        *((unsigned char *)data + *(idx + id0));
  return(cmp);
}

/*!
* \return	*((short *)data + *(idx + id1)) -
*		*((short *)data + *(idx + id0))
* \ingroup	AlgSort
* \brief	Inverse short comparison function for AlgHeapSortIdx().
*		Sorted data will have greatest entry first and smallest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortInvCmpIdxSFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((short *)data + *(idx + id1)) - *((short *)data + *(idx + id0));
  return(cmp);
}

/*!
* \return	*((int *)data + *(idx + id1)) -
*		*((int *)data + *(idx + id0))
* \ingroup	AlgSort
* \brief	Inverse int comparison function for AlgHeapSortIdx().
*		Sorted data will have greatest entry first and smallest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortInvCmpIdxIFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((int *)data + *(idx + id1)) - *((int *)data + *(idx + id0));
  return(cmp);
}

/*!
* \return	*((long *)data + *(idx + id1)) -
*		*((long *)data + *(idx + id0))
* \ingroup	AlgSort
* \brief	Inverse long comparison function for AlgHeapSortIdx().
*		Sorted data will have greatest entry first and smallest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortInvCmpIdxLFn(void *data, int *idx, int id0, int id1)
{
  int		cmp;

  cmp = *((long *)data + *(idx + id1)) - *((long *)data + *(idx + id0));
  return(cmp);
}

/*!
* \return	*((float *)data + *(idx + id1)) -
*		*((float *)data + *(idx + id0))
* \ingroup	AlgSort
* \brief	Inverse float comparison function for AlgHeapSortIdx().
*		Sorted data will have greatest entry first and smallest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortInvCmpIdxFFn(void *data, int *idx, int id0, int id1)
{
  int		cmp = 0;
  float		cmpF;

  cmpF = *((float *)data + *(idx + id1)) - *((float *)data + *(idx + id0));
  if(cmpF < 0.0)
  {
    cmp = -1;
  }
  else if(cmpF > 0.0)
  {
    cmp = 1;
  }
  return(cmp);
}

/*!
* \return	*((double *)data + *(idx + id1)) -
*		*((double *)data + *(idx + id0))
* \ingroup	AlgSort
* \brief	Inverse double comparison function for AlgHeapSortIdx().
*		Sorted data will have greatest entry first and smallest last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
int		AlgHeapSortInvCmpIdxDFn(void *data, int *idx, int id0, int id1)
{
  int		cmp = 0;
  double	cmpD;

  cmpD = *((double *)data + *(idx + id1)) - *((double *)data + *(idx + id0));
  if(cmpD < 0.0)
  {
    cmp = -1;
  }
  else if(cmpD > 0.0)
  {
    cmp = 1;
  }
  return(cmp);
}

