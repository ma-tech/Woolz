#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzBoundaryUtils.c
* \author       Bill Hill
* \date         September 2002
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Simple functions that operate on boundary lists.
* \ingroup	WlzBoundary
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>

static int			WlzBoundPolyCountFn(
				  WlzBoundList *bnd);
static void			WlzBoundObjToPolyFillArray(
				  WlzBoundList *bnd,
				  WlzPolygonDomain **polyAry,
				  int *idx);

/*!
* \return	Array of polygon domains.
* \ingroup	WlzBoundary
* \brief	Given a boundary list object returns a simple array of
*		polygon domains.
* \param	bndObj			Given boundary list object.
* \param	dstArySz		Destination ptr for array size.
* \param	dstPolyAry		Destination ptr for the array.
*/
WlzErrorNum	WlzBoundObjToPolyDomArray(WlzObject *bndObj, int *dstArySz,
					  WlzPolygonDomain ***dstPolyAry)
{
  int		idx,
  		polyCnt;
  WlzPolygonDomain **polyAry = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bndObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(bndObj->type != WLZ_BOUNDLIST)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(bndObj->domain.b == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dstArySz == NULL) || (dstPolyAry == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  /* Count number of polydomains in the boundary list. */
  polyCnt = WlzBoundPolyCount(bndObj->domain.b, &errNum);
  /* Allocate array for polygon domain pointers. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((polyAry = (WlzPolygonDomain **)AlcMalloc(sizeof(WlzPolygonDomain *) *
    						 polyCnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Fill in the Array. */
  if(errNum == WLZ_ERR_NONE)
  {
    idx = 0;
    WlzBoundObjToPolyFillArray(bndObj->domain.b, polyAry, &idx);
    *dstArySz = polyCnt;
    *dstPolyAry = polyAry;
  }
  return(errNum);
}

/*!
* \return	The number of polygon domains in the given boundary list.
* \ingroup	WlzBoundary
* \brief	Count the number of polygon domains in a boundary list.
* \param	bnd			Given boundary list.
* \param	dstErr			Destination ptr for error, may be NULL.
*/
int		WlzBoundPolyCount(WlzBoundList *bnd, WlzErrorNum *dstErr)
{
  int		cnt = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bnd == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    cnt = WlzBoundPolyCountFn(bnd);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cnt);
}

/*!
* \return	The number of polygon domains in the given boundary list
*		from the given entry point downwards..
* \ingroup	WlzBoundary
* \brief	Recursive function which actualy does the counting of the
* 		number of polygon domains in the given boundary list and
*		its children.
* \param	bnd
*/
static int	WlzBoundPolyCountFn(WlzBoundList *bnd)
{
  int		cnt = 0;

  if(bnd)
  {
    do
    {
      ++cnt;
      cnt += WlzBoundPolyCountFn(bnd->down);
      bnd = bnd->next;
    } while(bnd != NULL);
  }
  return(cnt);
}

/*!
* \return
* \ingroup	WlzBoundary
* \brief	Sets polygon domains from a boundary list to the given
*		array, which must be large enough. This is a recursive
*		function.
* \param	bnd			Given boundary list.
* \param	polyAry			The polygon domains array.
* \param	idx			Index incremented for each assignment.
*/
static void	WlzBoundObjToPolyFillArray(WlzBoundList *bnd,
					   WlzPolygonDomain **polyAry,
					   int *idx)
{
  if(bnd)
  {
    do
    {
      *(polyAry + *idx) = bnd->poly;
      ++*idx;
      WlzBoundObjToPolyFillArray(bnd->down, polyAry, idx);
      bnd = bnd->next;
    } while(bnd != NULL);
  }
}
