#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzIntervalDomScan.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for scanning through a Woolz object's
*		interval domain.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

/************************************************************************
*   Function   :WlzInitRasterScan					*
*   Returns    :WlzErrorNum: WLZ_ERR_NONE - success			*
*			WLZ_ERR_OBJECT_NULL, WLZ_ERR_OBJECT_TYPE,	*
*			WLZ_ERR_DOMAIN_NULL, WLZ_ERR_DOMAIN_TYPE.	*
*   Parameters :WlzObject *obj: object to be scanned, must be		*
*		WLZ_2D_DOMAINOBJ or WLZ_EMPTY_OBJ			*
*		WlzIntervalWSpace *iwsp: interval scanning work space.	*
*		WlzRasterDir raster: scanning direction.		*
*   Date       : Mon Oct 14 13:23:22 1996				*
*   Synopsis   :initialise interval scanning				*
************************************************************************/
WlzErrorNum 	WlzInitRasterScan(WlzObject *obj, WlzIntervalWSpace *iwsp,
		  		  WlzRasterDir raster)
{
  WlzIntervalDomain *ivdm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(iwsp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	if(obj->domain.i == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}    
	else
	{
	  ivdm = obj->domain.i;
	  iwsp->objaddr = obj;
	  iwsp->dmntype = ivdm->type;
	  switch(raster)
	  {
	    case WLZ_RASTERDIR_ILIC:
	      iwsp->lineraster = 1;
	      iwsp->colraster = 1;
	      break;
	    case WLZ_RASTERDIR_ILDC:
	      iwsp->lineraster = 1;
	      iwsp->colraster = -1;
	      break;
	    case WLZ_RASTERDIR_DLIC:
	      iwsp->lineraster = -1;
	      iwsp->colraster = 1;
	      break;
	    case WLZ_RASTERDIR_DLDC:
	      iwsp->lineraster = -1;
	      iwsp->colraster = -1;
	      break;
	    default:
	      errNum = WLZ_ERR_RASTERDIR_TYPE;
	      break;
	  }
	  iwsp->linbot = (iwsp->lineraster > 0) ? ivdm->line1 : ivdm->lastln;
	  iwsp->linpos = iwsp->linbot - iwsp->lineraster;
	  iwsp->linrmn = ivdm->lastln - ivdm->line1;
	  switch (iwsp->dmntype)
	  {
	    case WLZ_INTERVALDOMAIN_INTVL:
	      iwsp->intvln = (iwsp->lineraster > 0)?
	      		     ivdm->intvlines: ivdm->intvlines + iwsp->linrmn;
	      /* Pointer when incremented points to first line
	         to force nextinterval to move to first line */
	      iwsp->intvln -= iwsp->lineraster;
	      break;
	    case WLZ_INTERVALDOMAIN_RECT:
	      iwsp->intvln = NULL;
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  iwsp->intrmn = 0;
	  iwsp->colrmn = 0;
	  ++iwsp->linrmn;
	  iwsp->intdmn = ivdm;
	}
	break;
      case WLZ_EMPTY_OBJ:
	iwsp->objaddr = obj;
	iwsp->dmntype = WLZ_INTERVALDOMAIN_INTVL;     /* Used in WlzNextLine */
	iwsp->lineraster = 0;
	iwsp->colraster = 0;
	iwsp->linbot = 0;
	iwsp->linpos = 0;
	iwsp->linrmn = -1;    /* So WlzNextInterval & WlzNextLine do nothing */
	iwsp->intvln = NULL;
	iwsp->intrmn = 0;
	iwsp->colrmn = 0;
	iwsp->intdmn = NULL;
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }
  return(errNum);
}


/************************************************************************
*   Function   : WlzNextInterval					*
*   Returns    : WlzErrorNum:	WLZ_ERR_NONE 	- Next interval		*
*			        WLZ_ERR_EOO 	- No more intervals	*
*				*		- Error
*   Parameters :WlzIntervalWSpace *iwsp: interval scanning work space	*
*   Date       : Mon Oct 14 13:53:49 1996				*
*   Synopsis   :							*
* Get next interval in a standard object.
************************************************************************/
WlzErrorNum	WlzNextInterval(WlzIntervalWSpace *iwsp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(iwsp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(iwsp->linrmn < 0)
  {
    errNum = WLZ_ERR_EOO;
  }
  else
  {
    switch (iwsp->dmntype)
    {
      case WLZ_INTERVALDOMAIN_INTVL:
	iwsp->nwlpos = 0;
	if (iwsp->intrmn-- <= 0)
	{
	  /* move to next line */
	  do /* until non-empty line */
	  {
	    iwsp->linpos += iwsp->lineraster;
	    if(iwsp->linrmn-- <= 0)
	    {
	      errNum = WLZ_ERR_EOO;
	    }
	    else
	    {
	      iwsp->intvln += iwsp->lineraster;	  /* Next line of intervals */
	      iwsp->intrmn = iwsp->intvln->nintvs;
	      iwsp->nwlpos++;	   /* Count of newlines since last interval */
	    }
	  } while((iwsp->intrmn-- == 0) && (errNum == WLZ_ERR_NONE));
	  if(errNum == WLZ_ERR_NONE)
	  {
	    iwsp->intpos = (iwsp->colraster > 0)?
			   iwsp->intvln->intvs:
			   iwsp->intvln->intvs + iwsp->intrmn;
	  }
	}
	else
	{
	  iwsp->intpos += iwsp->colraster;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  iwsp->lftpos = iwsp->intpos->ileft + iwsp->intdmn->kol1;
	  iwsp->colpos = iwsp->lftpos ;
	  iwsp->rgtpos = iwsp->intpos->iright + iwsp->intdmn->kol1;
	}
	break;

      case WLZ_INTERVALDOMAIN_RECT:
	iwsp->linpos += iwsp->lineraster;
	if(iwsp->linrmn-- == 0)
	{
	  errNum = WLZ_ERR_EOO;
	}
	else
	{
	  iwsp->nwlpos = 1;
	  iwsp->intpos = (WlzInterval *) &(iwsp->intdmn->kol1);
	  iwsp->lftpos = iwsp->intdmn->kol1;
	  iwsp->colpos = iwsp->lftpos;
	  iwsp->rgtpos = iwsp->intdmn->lastkl;
	}
	break;
      default:
	/* This can't happen if the workspace has been initialised
	   using WlzInitRasterScan, WlzInitGreyScan or WlzInitGreyRasterScan
	 */
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;

    }
    if(errNum == WLZ_ERR_NONE)
    {
      iwsp->colrmn = iwsp->rgtpos - iwsp->lftpos + 1;
    }
  }
  return(errNum);
}


/************************************************************************
*   Function   : WlzNextLine						*
*   Returns    :							*
*   Parameters :							*
*   Date       : Mon Oct 14 13:56:57 1996				*
*   Synopsis   :							*
* moving to the line which is step away from the present line
* (beware when both nextline() and nextinterval() are used 
* since nextinterval() will locate the first interval in that line)
************************************************************************/
WlzErrorNum	WlzNextLine(WlzIntervalWSpace *iwsp, int step)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(iwsp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(iwsp->dmntype)
    {
      case WLZ_INTERVALDOMAIN_INTVL:
	iwsp->nwlpos = 0;
	iwsp->linpos += step*iwsp->lineraster;
	iwsp->linrmn -= step;
	if (iwsp->linrmn < 0)
	{
	  errNum = WLZ_ERR_EOO;
	}
	else
	{
	  iwsp->intvln += step * iwsp->lineraster; /* Next line of intervals */
	  iwsp->intrmn = iwsp->intvln->nintvs;
	  iwsp->intpos = (iwsp->colraster > 0)?
	  		 iwsp->intvln->intvs:
			 iwsp->intvln->intvs + iwsp->intrmn;
	  iwsp->lftpos = iwsp->intpos->ileft + iwsp->intdmn->kol1;
	  iwsp->colpos = iwsp->lftpos ;
	  iwsp->rgtpos = iwsp->intpos->iright + iwsp->intdmn->kol1;
	}
	break;
      case WLZ_INTERVALDOMAIN_RECT:
	iwsp->linpos += step * iwsp->lineraster;
	iwsp->linrmn -= step;
	if (iwsp->linrmn < 0)
	{
	  errNum = WLZ_ERR_EOO;
	}
	else
	{
	  iwsp->nwlpos = 1;
	  iwsp->intpos = (WlzInterval *) &(iwsp->intdmn->kol1);
	  iwsp->lftpos = iwsp->intdmn->kol1;
	  iwsp->colpos = iwsp->lftpos;
	  iwsp->rgtpos = iwsp->intdmn->lastkl;
	}
	break;
      default:
	/* This can't happen if the workspace has been initialised
	   using WlzInitRasterScan, WlzInitGreyScan or
	   WlzInitGreyRasterScan */
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      iwsp->colrmn = iwsp->rgtpos - iwsp->lftpos + 1;
    }
  }
  return(errNum);
}

/************************************************************************
*   Function   : WlzInitLineScan					*
*   Returns    :WlzErrorNum: WLZ_ERR_NONE, WLZ_ERR_OBJECT_NULL,		*
*		WLZ_ERR_OBJECT_TYPE, WLZ_ERR_DOMAIN_NULL,		*
		WLZ_ERR_DOMAIN_TYPE					*
*   Parameters :WlzObject *obj: object to be scanned			*
*		WlzIntervalWSpace *iwsp: interval scanning work space.	*
*		int raster: controls direction of scanning - see 	*
*		WlzInitRasterScan					*
*		int scale: same as step in WlzNextLine			*
*		int firstline:	starting line for scanning		*
*   Date       : Mon Oct 14 13:59:04 1996				*
*   Synopsis   :							*
* must call before the first call of WlzNextLine()
************************************************************************/
WlzErrorNum 	WlzInitLineScan(WlzObject *obj, WlzIntervalWSpace *iwsp,
				WlzRasterDir raster, int scale,
				int firstline)
{
  WlzIntervalDomain *ivdm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  if(iwsp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
	if(obj->domain.i == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  ivdm = obj->domain.i;
	  iwsp->objaddr = obj;
	  iwsp->dmntype = ivdm->type;
	  switch(raster)
	  {
	    case WLZ_RASTERDIR_ILIC:
	      iwsp->lineraster = 1;
	      iwsp->colraster = 1;
	      break;
	    case WLZ_RASTERDIR_ILDC:
	      iwsp->lineraster = 1;
	      iwsp->colraster = -1;
	      break;
	    case WLZ_RASTERDIR_DLIC:
	      iwsp->lineraster = -1;
	      iwsp->colraster = 1;
	      break;
	    case WLZ_RASTERDIR_DLDC:
	      iwsp->lineraster = -1;
	      iwsp->colraster = -1;
	      break;
	    default:
	      errNum = WLZ_ERR_RASTERDIR_TYPE;
	      break;
	  }
	  iwsp->linbot = firstline;
	  iwsp->linpos = iwsp->linbot - scale*iwsp->lineraster;
	  iwsp->linrmn = (iwsp->lineraster > 0)?
	   		 ivdm->lastln-firstline: firstline -ivdm->line1;
	  switch(iwsp->dmntype)
	  {
	    case  WLZ_INTERVALDOMAIN_INTVL:
	      iwsp->intvln = ivdm->intvlines + firstline - ivdm->line1;
	      /* Pointer when incremented points to first line
	         to force nextinterval to move to first line */
	      iwsp->intvln -= scale*iwsp->lineraster;
	      break;
	    case WLZ_INTERVALDOMAIN_RECT:
	      iwsp->intvln = NULL;
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  iwsp->intrmn = 0;
	  iwsp->colrmn = 0;
	  iwsp->linrmn++;
	  iwsp->intdmn = ivdm;
	}
	break;
      case WLZ_EMPTY_OBJ:
	iwsp->objaddr = obj;
	iwsp->dmntype = WLZ_INTERVALDOMAIN_INTVL;     /* Used in WlzNextLine */
	iwsp->lineraster = 0;
	iwsp->colraster = 0;
	iwsp->linbot = 0;
	iwsp->linpos = 0;
	iwsp->linrmn = -1;    /* So WlzNextInterval & WlzNextLine do nothing */
	iwsp->intvln = NULL;
	iwsp->intrmn = 0;
	iwsp->colrmn = 0;
	iwsp->intdmn = NULL;
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}
