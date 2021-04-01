#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFJSON_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzExtFFJSON.c
* \author       Bill Hill
* \date         March 2021
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2021],
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
* \brief	Read and write Woolz objects as JSON encoded Woolz data
* 		structures.
* \ingroup	WlzExtFF
*/

#include <stdio.h>
#include <Wlz.h>
#include <WlzExtFF.h>



static WlzErrorNum		WlzEffWriteJsn2DomObj(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum 		WlzEffWriteJsnIntervalDomain(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum		WlzEffWriteJsnValueTable(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum		WlzEffWriteJsnTiledValueTable(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum 		WlzEffWriteJsnPropertyList(
				  FILE *fP,
				  WlzPropertyList *pList);
static WlzErrorNum		WlzEffWriteJsnGreyValue(
				  FILE *fP,
				  WlzGreyV v,
				  WlzGreyType t,
				  char *trm);
static WlzErrorNum		WlzEffWriteJsnGreyValues(
				  FILE *fP,
				  WlzGreyP p,
				  WlzGreyType t,
				  int n,
				  const char *trm);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given file where the Woolz object
* 		is encoded using JSON.
* \param	fP			Input file stream.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzEffReadObjJsn(
				  FILE *fP,
				  WlzErrorNum *dstErr) 
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    // HACK TODO
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given file where the
* 		Woolz object is encoded using JSON.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum			WlzEffWriteObjJsn(
				  FILE *fP,
				  WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        errNum = WlzEffWriteJsn2DomObj(fP, obj);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz 2D domain object to the given file
* 		where the Woolz object is encoded using JSON.
* \param	fP			Output file stream.
* \param	obj			Given woolz object, known to be
* 					non-null and of type WLZ_2D_DOMAINOBJ.
*/
static WlzErrorNum WlzEffWriteJsn2DomObj(FILE *fP, WlzObject *obj)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(fprintf(fP, "{\n") < 0)
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    const char 	*s0;

    s0 = WlzStringFromObjTypeValue(obj->type, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if(fprintf(fP, "\"type\": \"%s\",\n", s0) < 0)
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffWriteJsnIntervalDomain(fP, obj);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, ",\n") < 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((obj->values.core == NULL) ||
       (WlzGreyTableIsTiled(obj->values.core->type) == 0))
    {
      errNum = WlzEffWriteJsnValueTable(fP, obj);
    }
    else
    {
      errNum = WlzEffWriteJsnTiledValueTable(fP, obj);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, ",\n") < 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffWriteJsnPropertyList(fP, obj->plist);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "}") < 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return       Woolz error number.
* \ingroup      WlzExtFF
* \brief        Writes the given Woolz objects interval domain to the given
* 		file where the Woolz domain is encoded using JSON.
* \param        fP                      Output file stream.
* \param        obj                    	Given woolz object which is known
* 					to have an interval domain.
*/
static WlzErrorNum WlzEffWriteJsnIntervalDomain(FILE *fP, WlzObject *obj)
{
  const char 	*trm[2] = {"", ","};
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(fprintf(fP, "\"domain\": {\n") < 0)
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  if((errNum == WLZ_ERR_NONE) && obj->domain.core)
  {
    WlzIntervalDomain *idom;

    idom = obj->domain.i;
    switch(idom->type)
    {
      case WLZ_INTERVALDOMAIN_INTVL: /* FALLTHROUGH */
      case WLZ_INTERVALDOMAIN_RECT:
        {
	  const char  *s0;

	  s0 = WlzStringFromObjDomainType(obj, NULL);
	  if(fprintf(fP,
	             "\"type\": \"%s\",\n"
		     "\"line1\": %d,\n"
		     "\"lastln\": %d,\n"
		     "\"kol1\": %d,\n"
		     "\"lastkl\": %d,\n"
		     "\"intvlines\": [\n",
		     s0,
		     idom->line1,
		     idom->lastln,
		     idom->kol1,
		     idom->lastkl) < 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	  if((errNum == WLZ_ERR_NONE) &&
	     (idom->type == WLZ_INTERVALDOMAIN_INTVL))
	  {
	    int		i,
	    		nlines;
	    WlzIntervalLine *ivln;

	    ivln = idom->intvlines;
	    nlines = idom->lastln - idom->line1 + 1;
	    for(i = 0; i < nlines; ++i)
	    {
	      int	j;

	      (void )fprintf(fP, "[\n");
	      for(j = 0; j < ivln->nintvs; ++j)
	      {

	        (void )fprintf(fP, "[%d,%d]%s\n",
		               ivln->intvs[j].ileft, ivln->intvs[j].iright,
			       trm[j + 1 < ivln->nintvs]);
	      }
	      (void )fprintf(fP, "]%s\n", trm[i + 1 < nlines]);
	    }
	    
	  }
	  if(fprintf(fP, "]") < 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "}") < 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief        Writes the given Woolz objects (2D) value table to the
* 		given file where the Woolz domain is encoded using JSON.
* \param        fP                      Output file stream.
* \param        obj                    	Given woolz object which is known
* 					to have a value table or NULL values.
*/
static WlzErrorNum		WlzEffWriteJsnValueTable(
				  FILE *fP,
				  WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fprintf(fP, "\"values\": {\n") < 0)
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  if((errNum == WLZ_ERR_NONE) && obj->values.core)
  {
    const char 	*s0;
    WlzGreyType	gType;
    WlzPixelV	bgd;

    gType = WlzGreyTableTypeToGreyType(obj->values.core->type, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      s0 = WlzStringFromGreyType(gType, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(fprintf(fP, "\"grey_type\": \"%s\"\n,", s0) < 0)
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      bgd = WlzGetBackground(obj, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzValueConvertPixel(&bgd, bgd, gType);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )fprintf(fP, "\"background\": ");
      errNum = WlzEffWriteJsnGreyValue(fP, bgd.v, gType, ",\n");
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )fprintf(fP, "\"values\": [");
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzIntervalWSpace     iWSp;
      WlzGreyWSpace         gWSp;
      const char 	    *trm[2] = {"", ","};

      if((errNum = WlzInitGreyScan(obj, &iWSp, &gWSp)) == WLZ_ERR_NONE)
      {

        while((errNum == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
        {
	  int		moreVal;

	  moreVal = (iWSp.linrmn > 0) || (iWSp.intrmn > 0);
	  errNum = WlzEffWriteJsnGreyValues(fP, gWSp.u_grintptr, gType,
	  				    iWSp.colrmn, trm[moreVal]);
	}
	(void )WlzEndGreyScan(&iWSp, &gWSp);
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )fprintf(fP, "]");
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "}") < 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  return(errNum);
}

static WlzErrorNum		WlzEffWriteJsnTiledValueTable(
				  FILE *fP,
				  WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  // HACK TODO
  return(errNum);
}

static WlzErrorNum 		WlzEffWriteJsnPropertyList(
 				  FILE *fP,
				  WlzPropertyList *pList)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  // HACK TODO
  if(fprintf(fP, "\"property_list\": {}") < 0)
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes a single grey value to the given file.
* \param	fP		File pointer.
* \param	v		Woolz grey value.
* \param	t		Woolz grey type.
* \param	trm		String to append following value.
*/
static WlzErrorNum		WlzEffWriteJsnGreyValue(
				  FILE *fP,
				  WlzGreyV v,
				  WlzGreyType t,
				  char *trm)
{
  int		c = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(t)
  {
    case WLZ_GREY_INT:
      c = fprintf(fP, "%d%s", v.inv, trm);
      break;
    case WLZ_GREY_SHORT:
      c = fprintf(fP, "%d%s", v.shv, trm);
      break;
    case WLZ_GREY_UBYTE:
      c = fprintf(fP, "%d%s", v.ubv, trm);
      break;
    case WLZ_GREY_FLOAT:
      c = fprintf(fP, "%g%s", v.flv, trm);
      break;
    case WLZ_GREY_DOUBLE:
      c = fprintf(fP, "%lg%s", v.dbv, trm);
      break;
    case WLZ_GREY_RGBA:
      {
	int	rgba[4];

	rgba[0] = WLZ_RGBA_RED_GET(v.rgbv);
	rgba[1] = WLZ_RGBA_GREEN_GET(v.rgbv);
	rgba[2] = WLZ_RGBA_BLUE_GET(v.rgbv);
	rgba[3] = WLZ_RGBA_ALPHA_GET(v.rgbv);
	c = fprintf(fP, "\"0x%02x%02x%02x%02x\"%s",
		    rgba[0], rgba[1], rgba[2], rgba[3], trm);
      }
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if((errNum == WLZ_ERR_NONE) && (c < 0))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes a multiple grey values to the given file.
* \param	fP		File pointer.
* \param	v		Woolz grey pointer for grey values.
* \param	t		Woolz grey type.
* \param	n		Number of grey values to write.
* \param	trm		String to append following values.
*/
static WlzErrorNum		WlzEffWriteJsnGreyValues(
				  FILE *fP,
				  WlzGreyP p,
				  WlzGreyType t,
				  int n,
				  const char *trm)
{
  int  		i,
  		l,
		k;
  const char 	*sep[2] = {0},
  		*lbk[2] = {0};
  const int	vol = 20;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  sep[0] = ",";
  sep[1] = trm;
  lbk[0] = "";
  lbk[1] = "\n";
  for(i = 0; i < n; ++i)
  {
    l = (i == n - 1);
    k = l || ((i % vol) == vol - 1);
    switch(t)
    {
      case WLZ_GREY_INT:
        (void )fprintf(fP, "%d%s%s", p.inp[i], sep[l], lbk[k]);
	break;
      case WLZ_GREY_SHORT:
        (void )fprintf(fP, "%d%s%s", p.shp[i], sep[l], lbk[k]);
	break;
      case WLZ_GREY_UBYTE:
        (void )fprintf(fP, "%d%s%s", p.ubp[i], sep[l], lbk[k]);
	break;
      case WLZ_GREY_FLOAT:
	(void )fprintf(fP, "%g%s%s", p.flp[i], sep[l], lbk[k]);
	break;
      case WLZ_GREY_DOUBLE:
	(void )fprintf(fP, "%lg%s%s", p.dbp[i], sep[l], lbk[k]);
	break;
      case WLZ_GREY_RGBA:
	{
	  int	rgba[4];

	  rgba[0] = WLZ_RGBA_RED_GET(p.rgbp[i]);
	  rgba[1] = WLZ_RGBA_GREEN_GET(p.rgbp[i]);
	  rgba[2] = WLZ_RGBA_BLUE_GET(p.rgbp[i]);
	  rgba[3] = WLZ_RGBA_ALPHA_GET(p.rgbp[i]);
	  (void )fprintf(fP, "\"0x%02x%02x%02x%02x\"%s%s",
		      rgba[0], rgba[1], rgba[2], rgba[3], sep[l], lbk[k]);
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (feof(fP) != 0))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

