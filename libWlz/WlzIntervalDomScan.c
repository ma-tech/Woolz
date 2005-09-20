#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzIntervalDomScan.c
* \author       Richard Baldock
* \date         September 2003
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
* \brief	Functions for scanning through an object's interval domain.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/* function:     WlzInitRasterScan    */
/*! 
* \ingroup      WlzDomainOps
* \brief        Initialise raster scanning.
*
* \return       Woolz error.
* \param    obj	Input object too be scanned.
* \param    iwsp	Interval scan workspace.
* \param    raster	Scanning direction.
* \par      Source:
*                WlzIntervalDomScan.c
*/
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

/* function:     WlzNextInterval    */
/*! 
* \ingroup      WlzDomainOps
* \brief        Get next interval in a standard object.
*
* \return       Woolz error.
* \param    iwsp	Interval scanning workspace.
* \par      Source:
*                WlzIntervalDomScan.c
*/
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

/* function:     WlzNextLine    */
/*! 
* \ingroup      WlzDomainOps
* \brief        Moving to the line which is step away from the present line
 (beware when both nextline() and nextinterval() are used 
 since nextinterval() will locate the first interval in that line)
*
* \return       Woolz error.
* \param    iwsp	Interval scanning workspace.
* \param    step	Distance to the required line (1 - next line etc.)
* \par      Source:
*                WlzIntervalDomScan.c
*/
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

/* function:     WlzInitLineScan    */
/*! 
* \ingroup      WlzDomainOps
* \brief        Initialise line scanning, must be called before the first
call of WlzNextLine().
*
* \return       Woolz error.
* \param    obj	Input object to be scanned.
* \param    iwsp	Interval scanning work space.
* \param    raster	Direction for the raster scan.
* \param    scale	Step for line increments.
* \param    firstline	Starting line for scanning.
* \par      Source:
*                WlzIntervalDomScan.c
*/
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
