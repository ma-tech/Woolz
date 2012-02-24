#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzWindow_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzWindow.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Weights the grey values of domain objects according
* 		to their spatial distribution using a weighting
* 		function.
* \ingroup	WlzValuesFilters
*/

#include <stdlib.h>
#include <Wlz.h>

static const char wlzWinFnNameBlackman[]  = "blackman",
		wlzWinFnNameHamming[]     = "hamming",
		wlzWinFnNameHanning[]     = "hanning",
		wlzWinFnNameParzen[]      = "parzen",
		wlzWinFnNameRectangle[]   = "rectangle",
		wlzWinFnNameWelch[]       = "welch",
		wlzWinFnNameUnspecified[] = "unspecified";

/*!
* \return	void
* \ingroup	WlzValuesFilters
* \brief	Applies the specified window function to the given Woolz
*		object.
*
*		The 1D window functions are applied along lines from
*		center of the elipse to its edge, all data outside of
*		the elipse are set to zero. The distance from the
*		center of an elipse (at the origin) to its edge passing
*		through soome point \f$P\f$ is given by
*		\f[
		  R_x R_y  \sqrt{\frac{1 + {(\frac{p_x}{p_y})}^2}
				      {R_y^2 +
				       {(R_x \frac{p_x}{p_y}}^2)}}
		\f]
		and the distance from its center to the point \f$p\f$ is
*	  	\f$\sqrt{p_x^2 + p_y^2}\f$
*		where:
*		\f$R_x\f$ and \f$R_y\f$	are the horizontal and vertical radii
*		of the axis aligned ellipse; \f$p_x\f$ and \f$p_y\f$ are the
*		coordinates of the point \f$P\f$.
*		The 1D window functions are applied along these lines
*		using the above distances.
*		This function, which is only called by WlzWindow()
*		assumes that all its parameters are valid.
* \param	obj			Given woolz object which MUST have a
*					non-null rectangular value table as
*					produced by WlzCutObjToBox2D().
* \param	center			Centre of the window with respect to
*					the given object.
* \param	radius			Radius of the window function.
* \param	winFn			Required windowing function.
*/
static void	WlzWindowApplyFn(WlzObject *obj, WlzIVertex2 center,
				 WlzIVertex2 radius, WlzWindowFnType winFn)
{

  WlzRectValues	*rVTab;
  WlzGreyType	vType;
  WlzGreyP	rVData,
  		tData;
  WlzIBox2	vLimits;
  WlzIVertex2	pos;
  WlzIVertex2	quadOff[2];
  double	tD0,
		tD1,
		distE,
		distP;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzWindowApplyFn FE %p {%d %d} {%d %d} %d\n",
	   obj, center.vtX, center.vtY, radius.vtX, radius.vtY, (int )winFn));
  if((winFn == WLZ_WINDOWFN_BLACKMAN) ||
     (winFn == WLZ_WINDOWFN_HAMMING) ||
     (winFn == WLZ_WINDOWFN_HANNING) ||
     (winFn == WLZ_WINDOWFN_PARZEN) ||
     (winFn == WLZ_WINDOWFN_WELCH))
  {
    rVTab = obj->values.r;
    rVData = rVTab->values;
    vType = WlzGreyTableTypeToGreyType(rVTab->type, NULL);
    if((vLimits.xMin = center.vtX - rVTab->kol1) < 0)
    {
      vLimits.xMin = 0;
    }
    if((vLimits.yMin = center.vtY - rVTab->line1) < 0)
    {
      vLimits.yMin = 0;
    }
    vLimits.xMax = rVTab->kol1 + rVTab->width - 1 - center.vtX;
    vLimits.yMax = rVTab->lastln - center.vtY;
    pos.vtY = 0;
    quadOff[1].vtY = quadOff[0].vtY = vLimits.yMin * rVTab->width;
    while((quadOff[0].vtY >= 0) || (quadOff[1].vtY >= 0))
    {
      pos.vtX = 0;
      quadOff[1].vtX = quadOff[0].vtX = vLimits.xMin;
      while((quadOff[0].vtX >= 0) || (quadOff[1].vtX >= 0))
      {
	tD0 = (double )(pos.vtX) / (double )(radius.vtX);
	tD1 = (double )(pos.vtY) / (double )(radius.vtY);
	if(((tD0 * tD0) + (tD1 *  tD1)) >= 1.0)
	{
	  tD0 = 0.0;
	}
	else
	{
	  if(pos.vtX == 0)
	  {
	    distE = radius.vtY;
	    distP = pos.vtY;
	  }
	  else if(pos.vtY == 0)
	  {
	    distE = radius.vtX;
	    distP = pos.vtX;
	  }
	  else
	  {
	    tD0 = (double )(pos.vtY) / (double )(pos.vtX);
	    tD0 *= tD0;
	    tD1 = (double )(radius.vtY * radius.vtY) +
		  ((double )(radius.vtX * radius.vtX) * tD0);
	    tD0 +=  1.0;
	    distE = (double )(radius.vtX * radius.vtY) * sqrt(tD0 / tD1);
	    distP = sqrt((double )(pos.vtX * pos.vtX) +
			 (double )(pos.vtY * pos.vtY));
	  }
	  switch(winFn)
	  {
	    case WLZ_WINDOWFN_BLACKMAN:
	      tD1 = WLZ_M_PI * distP / distE;
	      tD0 = 0.42 + (0.50 * cos(tD1)) + (0.08 * cos(2.0 * tD1));
	      break;
	    case WLZ_WINDOWFN_HAMMING:
	      tD1 = WLZ_M_PI * distP / distE;
	      tD0 = 0.54 + (0.46 * cos(tD1));
	      break;
	    case WLZ_WINDOWFN_HANNING:
	      tD1 = WLZ_M_PI * distP / distE;
	      tD0 = 0.50 + (0.50 * cos(tD1));
	      break;
	    case WLZ_WINDOWFN_PARZEN:
	      tD0 = 1.0 - (distP / distE);
	      break;
	    case WLZ_WINDOWFN_WELCH:
	      tD0 = distP / distE;
	      tD0 *= tD0;
	      tD0 = 1.0 - tD0;
	      break;
	    default:
              break;
	  }
	}
	switch(vType)
	{
	  case WLZ_GREY_INT:
	    if(quadOff[0].vtY >= 0)
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.inp = rVData.inp + quadOff[0].vtY + quadOff[0].vtX;
		*(tData.inp) = WLZ_NINT(*(tData.inp) * tD0);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.inp = rVData.inp + quadOff[0].vtY + quadOff[1].vtX;
		*(tData.inp) = WLZ_NINT(*(tData.inp) * tD0);
	      }
	    }
	    if((pos.vtY != 0) && (quadOff[1].vtY >= 0))
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.inp = rVData.inp + quadOff[1].vtY + quadOff[0].vtX;
		*(tData.inp) = WLZ_NINT(*(tData.inp) * tD0);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.inp = rVData.inp + quadOff[1].vtY + quadOff[1].vtX;
		*(tData.inp) = WLZ_NINT(*(tData.inp) * tD0);
	      }
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    if(quadOff[0].vtY >= 0)
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.shp = rVData.shp + quadOff[0].vtY + quadOff[0].vtX;
		*(tData.shp) = (short )WLZ_NINT(*(tData.shp) * tD0);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.shp = rVData.shp + quadOff[0].vtY + quadOff[1].vtX;
		*(tData.shp) = (short )WLZ_NINT(*(tData.shp) * tD0);
	      }
	    }
	    if((pos.vtY != 0) && (quadOff[1].vtY >= 0))
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.shp = rVData.shp + quadOff[1].vtY + quadOff[0].vtX;
		*(tData.shp) = (short )WLZ_NINT(*(tData.shp) * tD0);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.shp = rVData.shp + quadOff[1].vtY + quadOff[1].vtX;
		*(tData.shp) = (short )WLZ_NINT(*(tData.shp) * tD0);
	      }
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    if(quadOff[0].vtY >= 0)
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.ubp = rVData.ubp + quadOff[0].vtY + quadOff[0].vtX;
		*(tData.ubp) = (WlzUByte )((*(tData.ubp) * tD0) + 0.5);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.ubp = rVData.ubp + quadOff[0].vtY + quadOff[1].vtX;
		*(tData.ubp) = (WlzUByte )((*(tData.ubp) * tD0) + 0.5);
	      }
	    }
	    if((pos.vtY != 0) && (quadOff[1].vtY >= 0))
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.ubp = rVData.ubp + quadOff[1].vtY + quadOff[0].vtX;
		*(tData.ubp) = (WlzUByte )((*(tData.ubp) * tD0) + 0.5);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.ubp = rVData.ubp + quadOff[1].vtY + quadOff[1].vtX;
		*(tData.ubp) = (WlzUByte )((*(tData.ubp) * tD0) + 0.5);
	      }
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    if(quadOff[0].vtY >= 0)
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		*(rVData.flp + quadOff[0].vtY + quadOff[0].vtX) *= tD0;
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		*(rVData.flp + quadOff[0].vtY + quadOff[1].vtX) *= tD0;
	      }
	    }
	    if((pos.vtY != 0) && (quadOff[1].vtY >= 0))
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		*(rVData.flp + quadOff[1].vtY + quadOff[0].vtX) *= tD0;
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		*(rVData.flp + quadOff[1].vtY + quadOff[1].vtX) *= tD0;
	      }
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    if(quadOff[0].vtY >= 0)
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		*(rVData.dbp + quadOff[0].vtY + quadOff[0].vtX) *= tD0;
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		*(rVData.dbp + quadOff[0].vtY + quadOff[1].vtX) *= tD0;
	      }
	    }
	    if((pos.vtY != 0) && (quadOff[1].vtY >= 0))
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		*(rVData.dbp + quadOff[1].vtY + quadOff[0].vtX) *= tD0;
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		*(rVData.dbp + quadOff[1].vtY + quadOff[1].vtX) *= tD0;
	      }
	    }
	    break;
	  case WLZ_GREY_RGBA: /* RGBA to be done RAB */
	  default:
	    break;
	}
	++pos.vtX;
	--(quadOff[0].vtX);
	if(pos.vtX <= vLimits.xMax)
	{
	  ++(quadOff[1].vtX);
	}
	else
	{
	  quadOff[1].vtX = -1;
	}
      }
      ++pos.vtY;
      quadOff[0].vtY -= rVTab->width;
      if(pos.vtY <= vLimits.yMax)
      {
	quadOff[1].vtY += rVTab->width;
      }
      else
      {
	quadOff[1].vtY = -1;
      }
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzWindowApplyFn FX\n"));
}

/*!
* \return	New object.
* \ingroup	WlzValuesFilters
* \brief	Cuts a rectangular region from the given object and then
*		applies the specified 2D data windowing function.
*		The linkcount of the returned object is zero.
* \param	srcObj			Given source object.
* \param	winFn			Required windowing function.
* \param	org			Window origin with respect to the
*					given given source object.
* \param	rad			Window radius.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzWindow(WlzObject *srcObj, WlzWindowFnType winFn,
			      WlzIVertex2 org, WlzIVertex2 rad,
			      WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*dstObj = NULL;
  WlzGreyType	gType;
  WlzIBox2	cutBox;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzWindow FE %p %d {%d %d} {%d %d} %p\n",
	   srcObj, winFn, org.vtX, org.vtY, rad.vtX, rad.vtY, dstErr));
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	dstObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(srcObj->values.core == NULL)
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else if((rad.vtX <= 0) || (rad.vtY <= 0))
	{
	  dstObj = WlzMakeEmpty(&errNum);
	}
	else
	{
	  cutBox.xMin = org.vtX - rad.vtX;
	  cutBox.yMin = org.vtY - rad.vtY;
	  cutBox.xMax = org.vtX  + rad.vtX;
	  cutBox.yMax = org.vtY  + rad.vtY;
	  gType = WlzGreyTableTypeToGreyType(srcObj->values.core->type,
	  				     &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj = WlzCutObjToBox2D(srcObj, cutBox, gType,
				      0, 0.0, 0.0, &errNum);
	    if(dstObj && (errNum == WLZ_ERR_NONE))
	    {
	      if(dstObj->type == WLZ_2D_DOMAINOBJ)
	      {
		WlzWindowApplyFn(dstObj, org, rad, winFn);
	      }
	    }
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      WlzFreeObj(dstObj);
      dstObj = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzWindow FX %p\n",
	   dstObj));
  return(dstObj);
}

/*!
* \return	Constant name string or NULL if match fails.
* \ingroup	WlzValuesFilters
* \brief	Finds the appropriate name string for the given window
*		function value.
* \param	winFnValue		Window function value.
*/
const char 	*WlzWindowFnName(WlzWindowFnType winFnValue)
{
  const char	*winFnName = NULL;

  switch(winFnValue)
  {
    case WLZ_WINDOWFN_BLACKMAN:
      winFnName = wlzWinFnNameBlackman;
      break;
    case WLZ_WINDOWFN_HAMMING:
      winFnName = wlzWinFnNameHamming;
      break;
    case WLZ_WINDOWFN_HANNING:
      winFnName = wlzWinFnNameHanning;
      break;
    case WLZ_WINDOWFN_PARZEN:
      winFnName = wlzWinFnNameParzen;
      break;
    case WLZ_WINDOWFN_RECTANGLE:
      winFnName = wlzWinFnNameRectangle;
      break;
    case WLZ_WINDOWFN_WELCH:
      winFnName = wlzWinFnNameWelch;
      break;
    case WLZ_WINDOWFN_UNSPECIFIED:
    default:
      winFnName = wlzWinFnNameUnspecified;
      break;
  }
  return(winFnName);
}

/*!
* \return	Window function value. WLZ_WINDOWFN_UNSPECIFIED is returned
*		if no match is found.
* \ingroup	WlzValuesFilters
* \brief	Finds the appropriate value for the given window function
*		name string.
* \param	winFnName		Window function name string.
*/
WlzWindowFnType	WlzWindowFnValue(const char *winFnName)
{
  int		matchVal;
  WlzWindowFnType winFnValue = WLZ_WINDOWFN_UNSPECIFIED;
  
  if(WlzStringMatchValue(&matchVal, winFnName,
			 wlzWinFnNameBlackman,  WLZ_WINDOWFN_BLACKMAN,
			 wlzWinFnNameHamming,   WLZ_WINDOWFN_HAMMING,
			 wlzWinFnNameHanning,   WLZ_WINDOWFN_HANNING,
			 wlzWinFnNameParzen,    WLZ_WINDOWFN_PARZEN,
			 wlzWinFnNameRectangle, WLZ_WINDOWFN_RECTANGLE,
			 wlzWinFnNameWelch,     WLZ_WINDOWFN_WELCH,
			 NULL))
  {
    winFnValue = (WlzWindowFnType )matchVal;
  }
  return(winFnValue);
}
