#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzWindow.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Implements a spatial grey value windowing filter.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

static const char wlzWinFnNameBlackman[]  = "blackman",
		wlzWinFnNameHamming[]     = "hamming",
		wlzWinFnNameHanning[]     = "hanning",
		wlzWinFnNameParzen[]      = "parzen",
		wlzWinFnNameRectangle[]   = "rectangle",
		wlzWinFnNameWelch[]       = "welch",
		wlzWinFnNameUnspecified[] = "unspecified";

/************************************************************************
* Function:	WlzWindowApplyFn					*
* Returns:	void							*
* Purpose:	Applies the specified 2D data windowing function to the	*
*		given woolz object. In 1D the window functions are	*
*		  Blackman  0.42 + 0.50 * cos((2 * PI * n) / (N - 1)) +	*
*		  	    0.08 * cos((4 * PI * n) / (N - 1))		*
*		  Hamming   0.54 + 0.46 * cos((2 * PI * n) / (N - 1))	*
*		  Hanning   0.50 + 0.50 * cos((2 * PI * n) / (N - 1))	*
*		  Parzen    1.0 - |((2 * n) - (N - 1)) / (N + 1)|	*
* 		  Welch     1.0 - (((2 * n) - (N - 1)) / (N + 1))^2	*
*		The 1D window functions are applied along lines from	*
*		center of the elipse to its edge, all data outside of	*
*		the elipse are set to zero. The distance from the 	*
*		center of an elipse (at the origin) to its edge passing	*
*		through soome point p is given by 			*
*		  (Rx * Ry) * sqrt((1 + ((px / py)^2)) /		*
*				   ((Ry^2) + ((Rx * px / py)^2)))	*
*		and the distance from its center to the point p is 	*
*		  sqrt((px^2) + (py^2))					*
*		where 							*
*		  R[xy]	are the [xy] radii of the ellipse		*
*		  p[xy] are the coordinates of the point p.		*
*		The 1D window functions are applied along these lines	*
*		using the above distances.				*
*		This function, which is only called by WlzWindow()	*
*		assumes that all its parameters are valid.		*
* Global refs:	-							*
* Parameters:	WlzObject *obj:		Given woolz object which MUST	*
*					HAVE A A NON-NULL RECTANGULAR	*
*					VALUE TABLE as produced by	*
*					WlzCutObjToBox2D().		*
*		WlzIVertex2 center:	Windows center wrt the given	*
*					woolz object's coordinates.	*
*		WlzIVertex2 radius:	Windows radius.			*
*		WlzWindowFnType winFn:	Required windowing function.	*
************************************************************************/
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
	  ("WlzWindowApplyFn FE 0x%lx {%d %d} {%d %d} %d\n",
	   (unsigned long )obj, center.vtX, center.vtY, radius.vtX, radius.vtY,
	   (int )winFn));
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
		*(tData.shp) = WLZ_NINT(*(tData.shp) * tD0);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.shp = rVData.shp + quadOff[0].vtY + quadOff[1].vtX;
		*(tData.shp) = WLZ_NINT(*(tData.shp) * tD0);
	      }
	    }
	    if((pos.vtY != 0) && (quadOff[1].vtY >= 0))
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.shp = rVData.shp + quadOff[1].vtY + quadOff[0].vtX;
		*(tData.shp) = WLZ_NINT(*(tData.shp) * tD0);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.shp = rVData.shp + quadOff[1].vtY + quadOff[1].vtX;
		*(tData.shp) = WLZ_NINT(*(tData.shp) * tD0);
	      }
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    if(quadOff[0].vtY >= 0)
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.ubp = rVData.ubp + quadOff[0].vtY + quadOff[0].vtX;
		*(tData.ubp) = (unsigned int )((*(tData.ubp) * tD0) + 0.5);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.ubp = rVData.ubp + quadOff[0].vtY + quadOff[1].vtX;
		*(tData.ubp) = (unsigned int )((*(tData.ubp) * tD0) + 0.5);
	      }
	    }
	    if((pos.vtY != 0) && (quadOff[1].vtY >= 0))
	    {
	      if(quadOff[0].vtX >= 0)
	      {
		tData.ubp = rVData.ubp + quadOff[1].vtY + quadOff[0].vtX;
		*(tData.ubp) = (unsigned int )((*(tData.ubp) * tD0) + 0.5);
	      }
	      if((pos.vtX != 0) && (quadOff[1].vtX >= 0))
	      {
		tData.ubp = rVData.ubp + quadOff[1].vtY + quadOff[1].vtX;
		*(tData.ubp) = (unsigned int )((*(tData.ubp) * tD0) + 0.5);
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

/************************************************************************
* Function:	WlzWindow						*
* Returns:	WlzObject *:		New object with window applied	*
*					to it.				*
* Purpose:	Cuts a rectangular region from the given object and	*
*		then applies the specified 2D data windowing function.	*
* Notes:	The linkcount of the returned object is not set.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		WlzWindowFnType winFn:	Required windowing function.	*
*		WlzIVertex2 org:	Window origin wrt the given	*
*					source object's coordinates.	*
*		WlzIVertex2 rad:	Window radius.			*
*		WlzErrorNum *wlzErr:	Destination error pointer, may	*
*					be NULL.			*
************************************************************************/
WlzObject	*WlzWindow(WlzObject *srcObj, WlzWindowFnType winFn,
			      WlzIVertex2 org, WlzIVertex2 rad,
			      WlzErrorNum *wlzErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*dstObj = NULL;
  WlzGreyType	gType;
  WlzIBox2	cutBox;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzWindow FE 0x%lx %d {%d %d} {%d %d} 0x%lx\n",
	   (unsigned long )srcObj, winFn, org.vtX, org.vtY, rad.vtX, rad.vtY,
	   (unsigned long )wlzErr));
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
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzWindow FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}

/************************************************************************
* Function:	WlzWindowFnName						*
* Returns:	const char *:		Constant name string or NULL if	*
*					match fails.			*
* Purpose:	Finds the appropriate name string for the given window	*
*		function value.						*
* Global refs:	const char wlzWinFnNameBlackman[]: Window function name	*
*					strings.			*
*		const char wlzWinFnNameHamming[]:			*
*		const char wlzWinFnNameHanning[]:			*
*		const char wlzWinFnNameParzen[]:			*
*		const char wlzWinFnNameRectangle[]:			*
*		const char wlzWinFnNameWelch[]:				*
*		const char wlzWinFnNameUnspecified[]:			*
* Parameters:	WlzWindowFnType winFnValue: Window function value.	*
************************************************************************/
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

/************************************************************************
* Function:	WlzWindowFnValue					*
* Returns:	WlzWindowFnType:	Window function value, which is	*
*					WLZ_WINDOWFN_UNSPECIFIED if no	*
*					match is found.			*
* Purpose:	Finds the appropriate value for the given window	*
*		function name.						*
* Global refs:	const char wlzWinFnNameBlackman[]: Window function name	*
*					strings.			*
*		const char wlzWinFnNameHamming[]:			*
*		const char wlzWinFnNameHanning[]:			*
*		const char wlzWinFnNameParzen[]:			*
*		const char wlzWinFnNameRectangle[]:			*
*		const char wlzWinFnNameWelch[]:				*
*		const char wlzWinFnNameUnspecified[]:			*
* Parameters:	const char *winFnName:	Window function value.		*
************************************************************************/
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
