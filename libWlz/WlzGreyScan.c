#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyScan_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreyScan.c
* \author       Richard Baldock, Bill Hill
* \date         September 2003
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
* \brief	Object grey value scanning functions.
* \ingroup	WlzValuesUtils
*/

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

/*! 
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief        Initialise interval and grey scanning of standard object
*		in standard direction.
*
* \param    	obj			Object to be scanned.
* \param    	iwsp			Interval scanning workspace.
* \param    	gwsp			Grey value table scanning workspace.
*/
WlzErrorNum 
WlzInitGreyScan(WlzObject	*obj,	/* raster and tranpl of 0 provided */
		WlzIntervalWSpace	*iwsp,
		WlzGreyWSpace	*gwsp)
{
  return(WlzInitGreyRasterScan(obj, iwsp, gwsp, WLZ_RASTERDIR_ILIC, 0));
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief	This must be called when a grey workspace is no longer
* 		required as it frees resources of the workspace. 
* \param	iwsp			Interval scanning workspace.
* \param	gwsp			Grey value table scanning workspace.
*/
WlzErrorNum
WlzEndGreyScan(WlzIntervalWSpace *iwsp, WlzGreyWSpace *gwsp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((iwsp == NULL) || (gwsp == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(gwsp->gtable.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(iwsp->gryptr == gwsp)
  {
    WlzTiledValues *tv;
    WlzTiledValueBuffer *tvb;

    tv = gwsp->gtable.t;
    tvb = gwsp->tvb;
    if((WlzGreyTableIsTiled(tv->type)) && (tvb != NULL))
    {
      if(tvb->valid)
      {
	WlzTiledValueBufferFlush(tvb, tv);
      }
      WlzFreeTiledValueBuffer(tvb);
      gwsp->tvb = NULL;
    }
  }
  return(errNum);
}

/*! 
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief        As WlzInitGreyScan(), but with choice of raster direction
* 		and transplanting,
* \param    	obj			Input object to be scanned.
* \param    	iwsp			Interval scanning workspace.
* \param    	gwsp			Grey value table scanning workspace.
* \param    	raster			Direction for the raster scan.
* \param    	tranpl			Flag to allow overwriting of
* 					grey-values.
*/
WlzErrorNum 
WlzInitGreyRasterScan(WlzObject 	*obj,
		      WlzIntervalWSpace *iwsp,
		      WlzGreyWSpace 	*gwsp,
		      WlzRasterDir 	raster,
		      int 		tranpl)
{
  WlzErrorNum	errNum;

  if((errNum = WlzInitRasterScan(obj,iwsp,raster)) != WLZ_ERR_NONE)
  {
    return(errNum);
  }
  return(WlzInitGreyWSpace(obj,iwsp,gwsp,tranpl));
}

/*! 
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief        Attach grey workspace to interval scanning workspace.
* \param    	obj			Object to be scanned.
* \param    	iwsp			Interval scanning workspace.
* \param    	gwsp			Value table scanning workspace.
* \param    	tranpl			Flag to enable value transplanting.
*/
WlzErrorNum 
WlzInitGreyWSpace(WlzObject 		*obj,
		  WlzIntervalWSpace 	*iwsp,
		  WlzGreyWSpace 	*gwsp,
		  int 			tranpl)
{
  WlzValues 	vdmn;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

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
  else if((iwsp == NULL) || (gwsp == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    (void )memset(gwsp, 0, sizeof(WlzGreyWSpace));
    iwsp->gryptr = gwsp;
    gwsp->gtable = obj->values;
    vdmn = gwsp->gtable;
    gwsp->pixeltype = WlzGreyTableTypeToGreyType(vdmn.v->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gwsp->gdomaintype = WlzGreyTableTypeToTableType(vdmn.v->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(gwsp->gdomaintype)
    {

      case WLZ_GREY_TAB_RAGR:
	/* Synchronise grey table line pointer to interval domain
	   - pointer when incremented points to first line. */
	gwsp->gline = vdmn.v->vtblines + iwsp->linpos - vdmn.v->line1;
	break;

      case WLZ_GREY_TAB_RECT:
	/* Rectangular grey table does not have valueline pointers. */
	break;

      case WLZ_GREY_TAB_INTL:
	/* Synchronise grey table line pointer to interval domain
	   - pointer when incremented points to first line. */
	gwsp->gline = (WlzValueLine *) (vdmn.i->vil + iwsp->linpos -
	                                vdmn.i->line1);
	break;
      case WLZ_GREY_TAB_TILED:
	/* Create a buffer large enough to hold any single line of values. */
	if((gwsp->tvb = WlzMakeTiledValueBuffer(vdmn.t, &errNum)) != NULL)
        {
	  gwsp->tvb->pl = -1;
	}
        break;
      default:
	errNum = WLZ_ERR_VALUES_TYPE;
	break;

    }
    gwsp->intptr = iwsp;
    gwsp->tranpl = tranpl;
  }
  return(errNum);
}

/*! 
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief        Obtain next interval and its grey table.
*
* 		Note - cryptic documentation:
*	 	If tranpl and gvio=0 and scan already started (check line
*	 	number) then the previous grey interval must be collected,
*	 	if tranpl and gvio=1 then must output this interval (unless
*	 	scan finished).
* \return       Woolz error.
* \param    	iwsp			Interval scanning workspace - this
* 					must be initialised.
*/
WlzErrorNum 
WlzNextGreyInterval(WlzIntervalWSpace *iwsp)
{
  WlzGreyWSpace *gwsp = iwsp->gryptr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gwsp->tranpl && (gwsp->gvio == 0) &&
      (iwsp->linpos != (iwsp->linbot - iwsp->lineraster)))
  {
    errNum = WlzGreyInterval(iwsp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzNextInterval(iwsp);
  }
  if((errNum == WLZ_ERR_NONE) && (!gwsp->tranpl || (gwsp->gvio == 1)))
  {
    errNum = WlzGreyInterval(iwsp);
  }
  return(errNum);
}

/*! 
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief        Handle grey table for an interval. This must follow a
*		call to WlzNextInterval().
* \param    	iwsp			Interval scaning workspace.
*/
WlzErrorNum 
WlzGreyInterval(WlzIntervalWSpace *iwsp)
{
  int		i,
  		off;
  WlzGreyP	g;
  WlzValues	vdmn;
  WlzGreyWSpace	*gwsp = iwsp->gryptr;
  WlzValueLine *val;
  WlzValueIntervalLine *vil;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  g.v = NULL;
  vdmn = gwsp->gtable;
  switch(gwsp->gdomaintype)
  {
    case WLZ_GREY_TAB_RAGR:
      if(iwsp->nwlpos)
      {
	/* New line signal, move to next line of greytable. */
	gwsp->gline += iwsp->nwlpos*iwsp->lineraster;
      }
      /* Pointer to grey values. */
      off = iwsp->lftpos - vdmn.v->kol1 - gwsp->gline->vkol1;
      switch(gwsp->pixeltype)
      {
	case WLZ_GREY_INT:
	  g.inp = gwsp->gline->values.inp + off;
	  break;
	case WLZ_GREY_SHORT:
	  g.shp = gwsp->gline->values.shp + off;
	  break;
	case WLZ_GREY_UBYTE:
	  g.ubp = gwsp->gline->values.ubp + off;
	  break;
	case WLZ_GREY_FLOAT:
	  g.flp = gwsp->gline->values.flp + off;
	  break;
	case WLZ_GREY_DOUBLE:
	  g.dbp = gwsp->gline->values.dbp + off;
	  break;
	case WLZ_GREY_RGBA:
	  g.rgbp = gwsp->gline->values.rgbp + off;
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_RECT:
      /* Pointer to grey values. */
      off = ((iwsp->linpos - vdmn.r->line1) * vdmn.r->width) +
            iwsp->lftpos - vdmn.r->kol1;
      switch (gwsp->pixeltype)
      {
	case WLZ_GREY_INT:
	  g.inp = vdmn.r->values.inp + off;
	  break;
	case WLZ_GREY_SHORT:
	  g.shp = vdmn.r->values.shp + off;
	  break;
	case WLZ_GREY_UBYTE:
	  g.ubp = vdmn.r->values.ubp + off;
	  break;
	case WLZ_GREY_FLOAT:
	  g.flp = vdmn.r->values.flp + off;
	  break;
	case WLZ_GREY_DOUBLE:
	  g.dbp = vdmn.r->values.dbp + off;
	  break;
	case WLZ_GREY_RGBA:
	  g.rgbp = vdmn.r->values.rgbp + off;
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_INTL:
      if(iwsp->nwlpos)
      {				
	/* New line signal, move to next line of greytable. */
	gwsp->gline = (WlzValueLine *)(((WlzValueIntervalLine *)gwsp->gline) +
				       iwsp->nwlpos*iwsp->lineraster);
      }
      vil = (WlzValueIntervalLine *)gwsp->gline;
      val = vil->vtbint;
      for(i = 0; i < vil->nintvs; i++)
      {
	if((vdmn.i->kol1 + val->vkol1 <= iwsp->lftpos) &&
	   (vdmn.i->kol1 + val->vlastkl >= iwsp->rgtpos))
	{
	  off = iwsp->lftpos - vdmn.i->kol1 - val->vkol1;
	  switch(gwsp->pixeltype)
	  {
	    case WLZ_GREY_INT:
	      g.inp = val->values.inp + off;
	      break;
	    case WLZ_GREY_SHORT:
	      g.shp = val->values.shp + off;
	      break;
	    case WLZ_GREY_UBYTE:
	      g.ubp = val->values.ubp + off;
	      break;
	    case WLZ_GREY_FLOAT:
	      g.flp = val->values.flp + off;
	      break;
	    case WLZ_GREY_DOUBLE:
	      g.dbp = val->values.dbp + off;
	      break;
	    case WLZ_GREY_RGBA:
	      g.rgbp = val->values.rgbp + off;
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	  if(g.v)
	  {
	    break;
	  }
	}
	val++;
      }
      break;
    case WLZ_GREY_TAB_TILED:
      {
	WlzTiledValues *tv;
        WlzTiledValueBuffer *tvb;

	tv = vdmn.t;
	tvb = gwsp->tvb;
	if(tvb->valid)
	{
	  WlzTiledValueBufferFlush(tvb, tv);
	}
	tvb->pl = (tv->dim == 2)? 0: iwsp->plnpos - tv->plane1;
	tvb->ln = iwsp->linpos - tv->line1;
	off = tvb->kl[0] = iwsp->lftpos - tv->kol1;
	tvb->kl[1] = iwsp->rgtpos - tv->kol1;
	WlzTiledValueBufferFill(tvb, tv);
	switch (gwsp->pixeltype)
	{
	  case WLZ_GREY_INT:
	    g.inp = tvb->lnbuf.inp + off;
	    break;
	  case WLZ_GREY_SHORT:
	    g.shp = tvb->lnbuf.shp + off;
	    break;
	  case WLZ_GREY_UBYTE:
	    g.ubp = tvb->lnbuf.ubp + off;
	    break;
	  case WLZ_GREY_FLOAT:
	    g.flp = tvb->lnbuf.flp + off;
	    break;
	  case WLZ_GREY_DOUBLE:
	    g.dbp = tvb->lnbuf.dbp + off;
	    break;
	  case WLZ_GREY_RGBA:
	    g.rgbp = tvb->lnbuf.rgbp + off;
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
      break;
    default:
      errNum = WLZ_ERR_VALUES_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(gwsp->tranpl == 0)
    {
      /* Direct mode */
      gwsp->u_grintptr = g;
    }
    else
    {
      /* Transplant mode */
      /* 
       * MOD 12/1/92 - in transplant mode, we will always copy to or from
       * an integer array.  This change is made so that routines like
       * seqpar can more easily deal with differing pixel types.
       */
      switch(gwsp->gvio)
      {
	case 0 :		/* Input to object */
	  switch(gwsp->pixeltype)
	  {
	    case WLZ_GREY_INT:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		*g.inp++ = gwsp->u_grintptr.inp[i];
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		*g.shp++ = (short )WLZ_CLAMP(gwsp->u_grintptr.inp[i],
					     SHRT_MIN, SHRT_MAX);
	      }
	      break;
	    case WLZ_GREY_UBYTE:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		*g.ubp++ = (WlzUByte )WLZ_CLAMP(gwsp->u_grintptr.inp[i],
					        0, 255);
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		*g.flp++ = (float )(gwsp->u_grintptr.inp[i]);
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		*g.dbp++ = gwsp->u_grintptr.inp[i];
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		*g.rgbp++ = gwsp->u_grintptr.inp[i];
	      }
	      break;
	    default:
	      /* can't get here because grey-type already checked */
	      break;
	  }
	  break;
	case 1 :		/* Output from object */
	  switch(gwsp->pixeltype)
	  {
	    case WLZ_GREY_INT:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		gwsp->u_grintptr.inp[i] = *g.inp++;
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		gwsp->u_grintptr.inp[i] = *g.shp++;
	      }
	      break;
	    case WLZ_GREY_UBYTE:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		gwsp->u_grintptr.inp[i] = *g.ubp++;
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		gwsp->u_grintptr.inp[i] = (int )(*g.flp++);
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		gwsp->u_grintptr.inp[i] = (int )(*g.dbp++);
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      for (i=0; i<iwsp->colrmn; i++)
	      {
		gwsp->u_grintptr.inp[i] = *g.rgbp++;
	      }
	      break;
	    default:
	      /* can't get here because grey-type already checked */
	      break;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_INT_DATA;
	  break;
      }
    }
  }
  return(errNum);
}
