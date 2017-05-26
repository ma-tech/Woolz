#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzIntervalDomScan_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzIntervalDomScan.c
* \author       Richard Baldock
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
* \brief	Functions for scanning through an object's interval domain.
* \ingroup	WlzDomainOps
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
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
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

/* Classify the possible intersections. Here s (source), d
 *            * (destination) and o (overlap) are used in comment
 *            strings
 *                       * to represent the six possible interval
 *                       intersection cases. */

/*!
* \return	Value indicating which interval workspace should be advanced:
* 		0 for first, 1 for second.
* \ingroup	WlzDomainOps
* \brief	Classifies the possible intersections between the intervals
* 		of the two given interval workspaces.
*
* 		There are six possible cases which are returned as an
* 		intersection code. The destination interval pointer will
* 		have it's values set to cover the intersection where it
* 		exists. Given first interval F, second interval S and
* 		overlaps between the two intervals represented by O and
* 		gaps by - then the intersection codes returned are:
* 		<table border="1">
    		  <tr><td>Code</td> <td>Intersection</td></tr>
  		  <tr><td>0</td><td><tt>S-F</tt></td></tr>
  		  <tr><td>1</td><td><tt>F-S</tt></td></tr>
  		  <tr><td>2</td><td><tt>FOF</tt> or <tt>O</tt></td></tr>
  		  <tr><td>3</td><td><tt>SOS</tt></td></tr>
  		  <tr><td>4</td><td><tt>SOF</tt></td></tr>
  		  <tr><td>5</td><td><tt>FOS</tt></td></tr>
		</table>
* \param	dstItv			Destination pointer for the interval
* 					values, must be valid.
* \param	iWSp0			First interval workspace, must be
* 					valid. 
* \param	iWSp1			Second interval workspace, must be
* 					valid.
* \param	dstIsn			Destination pointer for the
* 					intersection code, may be NULL.
*/
int				WlzIWSpIntersection(
				  WlzInterval *dstItv,
				  WlzIntervalWSpace *iWSp0,
				  WlzIntervalWSpace *iWSp1,
				  int *dstIsn)
{
  int           t0 = 0,
		t1 = 0,
		isn = 0;
  WlzInterval 	itv = {0};

  if((iWSp1->linpos < iWSp0->linpos) ||
     (iWSp1->rgtpos < iWSp0->lftpos))                                 /* S F */
  {
    t0 = 1;
    isn = 0;
  }
  else if((iWSp1->linpos > iWSp0->linpos) ||
          (iWSp1->lftpos > iWSp0->rgtpos))                            /* F S */
  {
    t0 = 0;
    isn = 1;
  }
  else
  {
    t0 = iWSp1->lftpos >= iWSp0->lftpos;
    t1 = iWSp1->rgtpos <= iWSp0->rgtpos;
    if((t0 != 0) && (t1 != 0))                                   /* FOF || O */
    {
      isn = 2;
      itv.ileft = iWSp1->lftpos;
      itv.iright = iWSp1->rgtpos;
    }
    else if((t0 == 0) && (t1 == 0))                                   /* SOS */
    {
      isn = 3;
      itv.ileft = iWSp0->lftpos;
      itv.iright = iWSp0->rgtpos;
    }
    else
    {
      if(t0 == 0)                                                     /* SOF */
      {
	isn = 4;
	itv.ileft = iWSp0->lftpos;
	itv.iright = iWSp1->rgtpos;
      }
      else                                                            /* FOS */
      {
	isn = 5;
	itv.ileft = iWSp1->lftpos;
	itv.iright = iWSp0->rgtpos;
      }
      t0 = !t0;
    }
  }
  *dstItv = itv;
  if(dstIsn)
  {
    *dstIsn = isn;
  }
  return(t0);
}
