#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreySetHilbert_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreySetHilbert.c
* \author       Bill Hill
* \date         September 2011
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
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
* \brief	Functions for creating objects with integral grey values
* 		that are related to their Hilbert indices.
* \ingroup	WlzValuesUtils
*/
#include <Wlz.h>

static WlzObject 		*WlzGreyNewHilbertRankValues2D(
				  WlzObject *in,
				  unsigned int *rVal,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzGreyNewHilbertRankValues3D(
				  WlzObject *in,
				  unsigned int *rVal,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzGreySetHilbertRankValues2D(
				  WlzObject *obj,
				  unsigned int *rVal);
static WlzErrorNum		WlzGreySetHilbertRankValues3D(
				  WlzObject *obj,
				  unsigned int *rVal);
static void			WlzGreyNewHilbertSortUI(
				  unsigned int *iVal,
				  int nVal);
static int			WlzGreyHilbertRankFn0(
				  const void *v0,
				  const void *v1);
static int			WlzGreyHilbertRankFn1(
				  const void *v0,
				  const void *v1);

/*!
* \return       New grey value domain object with grey values set to the
* 		Hilbert rank of the voxel coordinates.
* \ingroup      WlzValuesUtils
* \brief	Creates a new domain object with integer values that increment
* 		throughout the object in Hilbert rank order. The input object's
* 		domain is used in the output object. Output object values
* 		start from 1.
* \param	in			Input domain object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzGreyNewHilbertRankValues(WlzObject *in,
					     WlzErrorNum *dstErr)
{
  unsigned int  val = 1;
  WlzObject     *out = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(in == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(in->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(in->type)
    {
      case WLZ_2D_DOMAINOBJ:
        out = WlzGreyNewHilbertRankValues2D(in, &val, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
        out = WlzGreyNewHilbertRankValues3D(in, &val, &errNum);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(out);
}

/*!
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief	Sets the values of domain object with int values so that they
* 		increment throughout the object in Hilbert rank order.
* \param	obj			Input domain object.
* \param	rVal			Minimum rank value on input and maximum
* 					rank value on output, NULL is
* 					equivalent to zero minimum rank value.
*/
WlzErrorNum	WlzGreySetHilbertRankValues(WlzObject *obj, unsigned int *rVal)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        errNum = WlzGreySetHilbertRankValues2D(obj, rVal);
        break;
      case WLZ_3D_DOMAINOBJ:
        errNum = WlzGreySetHilbertRankValues3D(obj, rVal);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  return(errNum);
}

/*!
* \return       New grey value domain object with grey values set to the
* 		Hilbert rank of the pixel coordinates.
* \ingroup      WlzValuesUtils
* \brief	Creates a new domain object with integer values that increment
* 		throughout the object in Hilbert rank order. The input object's
* 		domain is used in the output object.
* 		This function assumes that the pixels of an integer power
* 		of two sided square large enough to enclose the input object
* 		can be enumerated using a native integer.
* \param	in			Input domain object.
* \param	rVal			Minimum rank value on input and maximum
* 					rank value on output, NULL is
* 					equivalent to zero minimum rank value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzGreyNewHilbertRankValues2D(WlzObject *in,
 				                unsigned int *rVal,
				                WlzErrorNum *dstErr)
{
  WlzObject     *out = NULL;
  WlzObjectType gTT;
  WlzPixelV     bgd;
  WlzValues     values;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  bgd.type = WLZ_GREY_INT;
  bgd.v.inv = 0;
  values.core = NULL;
  gTT = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    values.v = WlzNewValueTb(in, gTT, bgd, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    out = WlzMakeMain(in->type, in->domain, values, in->plist, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGreySetHilbertRankValues2D(out, rVal);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(out != NULL)
    {
      (void )WlzFreeObj(out);
      out = NULL;
    }
    else if(values.core != NULL)
    {
      (void )WlzFreeValues(values);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(out);
}

/*!
* \return       New grey value domain object with grey values set to the
* 		Hilbert rank of the voxel coordinates.
* \ingroup      WlzValuesUtils
* \brief	Creates a new domain object with integer values that increment
* 		throughout the object in Hilbert rank order. The input object's
* 		domain is used in the output object.
* 		This function assumes that the voxels of an integer power
* 		of two sided cube large enough to enclose the input object
* 		can be enumerated using a native integer.
* \param	in			Input domain object.
* \param	rVal			Minimum rank value on input and maximum
* 					rank value on output, NULL is
* 					equivalent to zero minimum rank value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzGreyNewHilbertRankValues3D(WlzObject *in,
				                unsigned int *rVal,
				                WlzErrorNum *dstErr)
{
  WlzObjectType gTT;
  WlzObject     *out = NULL;
  WlzPixelV     bgd;
  WlzValues     values;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  bgd.type = WLZ_GREY_INT;
  bgd.v.inv = 0;
  values.core = NULL;
  gTT = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, &errNum);
  values.vox = WlzNewValuesVox(in, gTT, bgd, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    out = WlzMakeMain(in->type, in->domain, values, in->plist, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGreySetHilbertRankValues3D(out, rVal);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(out)
    {
      (void )WlzFreeObj(out);
      out = NULL;
    }
    else if(values.core != NULL)
    {
      (void )WlzFreeVoxelValueTb(values.vox);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(out);
}

/*!
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief	Sets the values of the given object to integer values that
*		increment throughout the object in Hilbert rank order.
* 		This function assumes that the pixels of an integer power
* 		of two sided square large enough to enclose the given object
* 		can be enumerated using a native integer.
* \param	obj			Given object.
* \param	rVal			Minimum rank value on input and maximum
* 					rank value on output, NULL is
* 					equivalent to zero minimum rank value.
*/
static WlzErrorNum WlzGreySetHilbertRankValues2D(WlzObject *obj,
				                 unsigned int *rVal)
{
  int		nB,
		minVal,
  		nVal;
  unsigned int  i;
  unsigned int  *iV,
  		*iVal = NULL;
  WlzIntervalWSpace iWsp;
  WlzGreyWSpace gWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find number of foreground values required and then allocate a pair
   * of unsigned integers for each foreground value. */
  nVal = WlzArea(obj, &errNum);
  minVal = (rVal == NULL)? 0: *rVal;
  if(errNum == WLZ_ERR_NONE)
  {
    nB = ALG_MAX(obj->domain.i->lastln - obj->domain.i->line1,
                 obj->domain.i->lastkl - obj->domain.i->kol1) + 1;
    nB = AlgBitMostSigSet((unsigned long )nB);
    if(nB > (8 * sizeof(unsigned int) / 2))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((iVal = (unsigned int *)
                AlcMalloc(sizeof(unsigned int) * 2 * nVal)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* For each forground value compute it's Hilbert index within a
   * square that has a side an integer power of two and is larger
   * than the given object's domain. This gives an array of integer
   * pairs {Hilbert index, image index}. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitRasterScan(obj, &iWsp, WLZ_RASTERDIR_ILIC);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    i = 0;
    iV = iVal;
    while((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE)
    {
      int	k,
      		nK;
      unsigned int h[2],
      		   c[2];

      c[1] = iWsp.linpos;
      nK = iWsp.rgtpos - iWsp.lftpos + 1;
      for(k = 0; k < nK; ++k)
      {
        c[0] = iWsp.lftpos + k;
	AlgHilbertIndex(h, c, 2, nB);
	*iV++ = (h[1] << nB) | h[0];
	*iV++ = i++;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  /* Sort the integer pairs by Hilbert index, change the Hilbert index to
   * it's rank and then re-sort by image index. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzGreyNewHilbertSortUI(iVal, nVal);
  }
  /* Set object's values to be the Hilbert rank of their coordinate. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(obj, &iWsp, &gWsp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    i = 0;
    iV = iVal;
    while((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE)
    {
      int	k,
      		nK;
      int	*gP;

      gP = gWsp.u_grintptr.inp;
      nK = iWsp.rgtpos - iWsp.lftpos + 1;
      for(k = 0; k < nK; ++k)
      {
	*gP++ = minVal + *iV;
	iV += 2;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  AlcFree(iVal);
  if(errNum == WLZ_ERR_NONE)
  {
    if(rVal)
    {
      *rVal = minVal + nVal - 1;
    }
  }
  return(errNum);
}

/*!
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief	Sets the values of the given object to integer values that
*		increment throughout the object in Hilbert rank order.
* 		This function assumes that the voxels of an integer power
* 		of two sided cube large enough to enclose the given object
* 		can be enumerated using a native integer.
* \param	obj			Given object.
* \param	rVal			Minimum rank value on input and maximum
* 					rank value on output, NULL is
* 					equivalent to zero minimum rank value.
*/
static WlzErrorNum WlzGreySetHilbertRankValues3D(WlzObject *obj,
					         unsigned int *rVal)
{
  int		nB,
		nP,
  		nVal,
		minVal;
  unsigned int  i;
  unsigned int  *iV,
  		*iVal = NULL;
  WlzDomain	*dom;
  WlzValues	*val;
  WlzIntervalWSpace iWsp;
  WlzGreyWSpace gWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find number of foreground values required and then allocate a pair
   * of unsigned integers for each foreground value. */
  minVal = (rVal == NULL)? 0: *rVal;
  nVal = WlzVolume(obj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    nB = ALG_MAX3(obj->domain.p->lastpl - obj->domain.p->plane1,
                  obj->domain.p->lastln - obj->domain.p->line1,
                  obj->domain.p->lastkl - obj->domain.p->kol1) + 1;
    nB = AlgBitMostSigSet((unsigned long )nB);
    if(nB > (8 * sizeof(unsigned int) / 3))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((iVal = (unsigned int *)
                AlcMalloc(sizeof(unsigned int) * 2 * nVal)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* For each forground value compute it's Hilbert index within a
   * cube that has a side an integer power of two and is larger
   * than the given object's domain. This gives an array of integer
   * pairs {Hilbert index, image index}. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		p;
    unsigned int c[3];

    i = 0;
    iV = iVal;
    dom = obj->domain.p->domains;
    nP = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
    for(p =  0; p < nP; ++p)
    {
      if(dom[p].core != NULL)
      {
	WlzValues dumVal;
        WlzObject *obj2;

	dumVal.core = NULL;
	c[2] = iWsp.linpos;
	obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom[p], dumVal,
	                   NULL, NULL, &errNum);
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzInitRasterScan(obj2, &iWsp, WLZ_RASTERDIR_ILIC);
	}
	while((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE)
	{
	  int	    k,
		    nK;
	  unsigned int h[3];

	  c[1] = iWsp.linpos;
	  nK = iWsp.rgtpos - iWsp.lftpos + 1;
	  for(k = 0; k < nK; ++k)
	  {
	    c[0] = iWsp.lftpos + k;
	    AlgHilbertIndex(h, c, 3, nB);
	    *iV++ = ((h[2] << nB) | (h[1]) << nB) | h[0];
	    *iV++ = i++;
	  }
	}
        (void )WlzFreeObj(obj2);
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
	else if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
  }
  /* Sort the integer pairs by Hilbert index, change the Hilbert index to
   * it's rank and then re-sort by image index. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzGreyNewHilbertSortUI(iVal, nVal);
  }
  /* Set the object's values to be the Hilbert rank of their coordinate. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		p;

    iV = iVal;
    dom = obj->domain.p->domains;
    val = obj->values.vox->values;
    for(p =  0; p < nP; ++p)
    {
      if(dom[p].core != NULL)
      {
        WlzObject *obj2;

	obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom[p], val[p],
	                   NULL, NULL, &errNum);
        if(errNum == WLZ_ERR_NONE)
	{
          errNum = WlzInitGreyScan(obj2, &iWsp, &gWsp);
	}
	while((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE)
	{
	  int	k,
		nK;
	  int	*gP;

	  gP = gWsp.u_grintptr.inp;
	  nK = iWsp.rgtpos - iWsp.lftpos + 1;
	  for(k = 0; k < nK; ++k)
	  {
	    *gP++ = minVal + *iV;
	    iV += 2;
	  }
	}
        (void )WlzFreeObj(obj2);
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
	else if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }

  }
  AlcFree(iVal);
  if(errNum == WLZ_ERR_NONE)
  {
    if(rVal)
    {
      *rVal = minVal + nVal - 1;
    }
  }
  return(errNum);
}
/*!
* \ingroup	WlzValuesUtils
* \brief	Sorts  integer pairs by Hilbert index, changes the Hilbert
* 		index to it's rank and then re-sorts by image index.
* \param	iVal			The integer value pairs.
* \param	nVal			The number of value pairs.
*/
static void	WlzGreyNewHilbertSortUI(unsigned int *iVal, int nVal)
{
  int		i;
  unsigned int	*iV;

  iV = iVal;
  qsort(iVal, nVal, 2 * sizeof(unsigned int), WlzGreyHilbertRankFn0);
  for(i = 0; i < nVal; ++i)
  {
    *iV = i;
    iV += 2;
  }
  qsort(iVal, nVal, 2 * sizeof(unsigned int), WlzGreyHilbertRankFn1);
}

/*!
* \return	Signed int for qsort().
* \ingroup	WlzValuesUtils
* \brief	Does a comparison for qsort(). Compares unsigned integers
* 		with zero offset.
* \param	v0			To cast to unsigned int.
* \param	v1			To cast to unsigned int.
*/
static int	WlzGreyHilbertRankFn0(const void *v0, const void *v1)
{
  int		cmp = 0;
  unsigned int	*u0,
                *u1;

  u0 = (unsigned int *)v0;
  u1 = (unsigned int *)v1;
  cmp = u0[0] - u1[0];
  return(cmp);
}

/*!
* \return	Signed int for qsort().
* \ingroup	WlzValuesUtils
* \brief	Does a comparison for qsort(). Compares unsigned integers
* 		with an offset of one.
* \param	v0			To cast to unsigned int.
* \param	v1			To cast to unsigned int.
*/
static int	WlzGreyHilbertRankFn1(const void *v0, const void *v1)
{
  int		cmp = 0;
  unsigned int	*u0,
                *u1;

  u0 = (unsigned int *)v0;
  u1 = (unsigned int *)v1;
  cmp = u0[1] - u1[1];
  return(cmp);
}
