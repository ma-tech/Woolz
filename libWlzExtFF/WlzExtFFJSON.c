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
#include <string.h>
#include <limits.h>
#include <float.h>
#include <cJSON.h>
#include <Wlz.h>
#include <WlzExtFF.h>


static char			*WlzEffReadJsnFile(
				  FILE *fP,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzEffParseJsn2DomObj(
				  const cJSON *jObj,
				  WlzErrorNum *dstErr);
static WlzIntervalDomain	*WlzEffParseJsnIntervalDomain(
				  const cJSON *jDom,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzEffParseJsnRGBA(
				  const char *str,
				  WlzUInt *rgba);
static WlzErrorNum		WlzEffParseJsnValueTable(
				  WlzObject *obj,
				  const cJSON *jVal);
static WlzErrorNum		WlzEffParseJsnPropertyList(
				  WlzObject *obj,
				  const cJSON *jPrp);
static WlzErrorNum		WlzEffWriteJsnMagic(
				  FILE *fP);
static WlzErrorNum		WlzEffWriteJsn2DomObj(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum 		WlzEffWriteJsnIntervalDomain(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum		WlzEffWriteJsnValueTable(
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
  char 		*buf = NULL;
  cJSON 	*jRoot = NULL;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    buf = WlzEffReadJsnFile(fP, &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && ((jRoot = cJSON_Parse(buf)) == NULL))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    const cJSON *jMag,
    		*jVer,
		*jObj;

    jMag = cJSON_GetObjectItemCaseSensitive(jRoot, "magic");
    jVer = cJSON_GetObjectItemCaseSensitive(jRoot, "version");
    jObj = cJSON_GetObjectItemCaseSensitive(jRoot, "object");
    if((jMag == NULL) || !cJSON_IsString(jMag) || !(jMag->valuestring) ||
       strcmp("Woolz", jMag->valuestring) ||
       (jVer == NULL) || !cJSON_IsString(jVer) || !(jVer->valuestring) ||
       (jObj == NULL) || !cJSON_IsObject(jObj))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      const cJSON *jType = NULL;
      WlzObjectType wType = WLZ_NULL;

      jType = cJSON_GetObjectItemCaseSensitive(jObj, "type");
      if((jType == NULL) || !cJSON_IsString(jType) || !(jType->valuestring))
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
        wType = WlzStringToObjType(jType->valuestring, NULL);
	switch(wType)
	{
	  case WLZ_EMPTY_OBJ:
	    obj = WlzMakeEmpty(&errNum);
	    break;
	  case WLZ_2D_DOMAINOBJ:
	    obj = WlzEffParseJsn2DomObj(jObj, &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
    }
  }
  AlcFree(buf);
  cJSON_Delete(jRoot);
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
    if(fprintf(fP, "{\n") < 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else
    {
      errNum = WlzEffWriteJsnMagic(fP);
    }
  }
  if(errNum == WLZ_ERR_NONE)
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
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "}\n") < 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes identification and version strings.
* \param	fP			Output file stream.
*/
static WlzErrorNum		WlzEffWriteJsnMagic(
				  FILE *fP)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fprintf(fP,
             "\"magic\": \"Woolz\",\n"
	     "\"version\": \"%s\",\n",
	     WlzVersion()) < 0)
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	The file contents as a null terminated string or NULL on
* 		error.
* \ingroup	WlzExtFF
* \brief	Reads the entire contents of the given file into a character
* 		buffer.	The buffer should be freed using AlcFree(0 when no
* 		longer required.
* \param	fP			Input file stream.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static char			*WlzEffReadJsnFile(
				  FILE *fP,
				  WlzErrorNum *dstErr)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  long		len = 0;
  char 		*buf = NULL;

  if((fseek(fP, 0, SEEK_END) != 0) || ((len = ftell(fP)) <= 0) ||
     (fseek(fP, 0, SEEK_SET) != 0))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if((buf = (char *)AlcMalloc((len + 1) * sizeof(char))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else if(fread(buf, sizeof(char), len, fP) != len)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    buf[len] = '\0';
  }
  else
  {
    AlcFree(buf);
    buf = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(buf);
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
static WlzErrorNum 		WlzEffWriteJsn2DomObj(
				  FILE *fP,
				  WlzObject *obj)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(fprintf(fP, "\"object\": {\n") < 0)
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
    if(obj->values.core)
    {
      errNum = WlzEffWriteJsnValueTable(fP, obj);
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
static WlzErrorNum 		WlzEffWriteJsnIntervalDomain(
				  FILE *fP,
				  WlzObject *obj)
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
	      ++ivln;
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

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes the property list of an object to the given file.
* \todo		Property lists are not yet written to or parsed from JSON
* 		encoded files.
* \param        fP                      Output file stream.
* \param        pList                  	Given woolz property list or NULL
* 					for an empty property list.
*/
static WlzErrorNum 		WlzEffWriteJsnPropertyList(
 				  FILE *fP,
				  WlzPropertyList *pList)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

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

/*!
* \return	New Woolz 2D domain object.
* \ingroup	WlzExtFF
* \brief	Parses a Woolz 2D domain object from the given cJSON object.
* \param	jObj			Given cJSON object which encodes a
* 					Woolz 2D domain object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzEffParseJsn2DomObj(
				  const cJSON *jObj,
				  WlzErrorNum *dstErr)
{
  const cJSON	*jDom,
  		*jVal,
		*jPrp;
  WlzObject	*obj = NULL;
  WlzDomain	dom = {0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  jDom = cJSON_GetObjectItemCaseSensitive(jObj, "domain");
  jVal = cJSON_GetObjectItemCaseSensitive(jObj, "values");
  jPrp = cJSON_GetObjectItemCaseSensitive(jObj, "property_list");
  if((jDom == NULL) || !cJSON_IsObject(jDom) ||
     (jVal == NULL) || !cJSON_IsObject(jVal) ||
     (jPrp == NULL) || !cJSON_IsObject(jPrp))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    if((dom.i = WlzEffParseJsnIntervalDomain(jDom, &errNum)) != NULL)
    {
      WlzValues	val = {0};

      obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom, val, NULL, NULL, &errNum);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffParseJsnValueTable(obj, jVal);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffParseJsnPropertyList(obj, jPrp);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(obj)
    {
      (void )WlzFreeObj(obj);
      obj = NULL;
    }
    else if(dom.core)
    {
      (void )WlzFreeDomain(dom);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Interval domain or NULL on error.
* \ingroup	WlzExtFF
* \brief	Continues to parse a JSON encoded object for it's interval
* 		domain.
* \param	jDom		Pointer to cJSON interval domain.
* \param	dstErr		Destination error pointer, may be NULL.
*/
static WlzIntervalDomain	*WlzEffParseJsnIntervalDomain(
				  const cJSON *jDom,
				  WlzErrorNum *dstErr)
{
  const cJSON 	*jType,
  		*jLn1,
		*jLLn,
		*jKl1,
		*jLKl,
		*jItvLns;
  WlzObjectType wType = WLZ_NULL;
  WlzIntervalDomain *iDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  jType = cJSON_GetObjectItemCaseSensitive(jDom, "type");
  if((jType == NULL) || !cJSON_IsString(jType) || !(jType->valuestring))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    wType = WlzStringToObjDomainType(jType->valuestring, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(wType)
    {
      case WLZ_INTERVALDOMAIN_INTVL: /* FALLTHROUGH */
      case WLZ_INTERVALDOMAIN_RECT:
	jLn1 = cJSON_GetObjectItemCaseSensitive(jDom, "line1");
	jLLn = cJSON_GetObjectItemCaseSensitive(jDom, "lastln");
	jKl1 = cJSON_GetObjectItemCaseSensitive(jDom, "kol1");
	jLKl = cJSON_GetObjectItemCaseSensitive(jDom, "lastkl");
	jItvLns = cJSON_GetObjectItemCaseSensitive(jDom, "intvlines");
	if((jLn1 == NULL) || !cJSON_IsNumber(jLn1) ||
	   (jLLn == NULL) || !cJSON_IsNumber(jLLn) ||
	   (jKl1 == NULL) || !cJSON_IsNumber(jKl1) ||
	   (jLKl == NULL) || !cJSON_IsNumber(jLKl) ||
	   (jItvLns == NULL) || !cJSON_IsArray(jItvLns))
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  int	l1,
	      ll,
	      k1,
	      kl;

	  l1 = ALG_NINT(jLn1->valuedouble);
	  ll = ALG_NINT(jLLn->valuedouble);
	  k1 = ALG_NINT(jKl1->valuedouble);
	  kl = ALG_NINT(jLKl->valuedouble);
	  if((l1 > ll) || (k1 > kl))
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	  else
	  {
	    iDom = WlzMakeIntervalDomain(wType, l1, ll, k1, kl, &errNum);
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (wType == WLZ_INTERVALDOMAIN_INTVL))
	{
	  int		nItv = 0,
			nItvLn;
	  const cJSON	*jItvLn;
	  WlzInterval	*itvs = NULL;

	  nItvLn = cJSON_GetArraySize(jItvLns);
	  if(iDom->lastln - iDom->line1 + 1 != nItvLn)
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	  else
	  {
	    cJSON_ArrayForEach(jItvLn, jItvLns)
	    {
	      if(jItvLn && cJSON_IsArray(jItvLn))
	      {
		nItv += cJSON_GetArraySize(jItvLn);
	      }
	      else
	      {
		errNum = WLZ_ERR_DOMAIN_DATA;
		break;
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if((itvs = (WlzInterval *)
		       AlcMalloc(nItv * sizeof(WlzInterval))) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	i,
		j;

	    i = 0;
	    jItvLn = jItvLns->child;
	    for(j = 0; j < nItvLn; ++j)
	    {
	      int	k,
			nItvInLn;
	      const cJSON *jItv;
	      WlzIntervalLine *itvLn;

	      itvLn = iDom->intvlines + j;
	      nItvInLn = cJSON_GetArraySize(jItvLn);
	      itvLn->nintvs = nItvInLn;
	      itvLn->intvs = itvs + i;
	      jItv = jItvLn->child;
	      for(k = 0; k < nItvInLn; ++k)
	      {
		WlzInterval	*itv;

		itv = itvs + i;
		if(jItv && cJSON_IsArray(jItv) &&
		    (cJSON_GetArraySize(jItv) == 2))
		{
		  const cJSON  *jLft,
			*jRgt;

		  jLft = jItv->child;
		  jRgt = jLft->next;
		  if(cJSON_IsNumber(jLft) && cJSON_IsNumber(jRgt))
		  {
		    itv->ileft  = ALG_NINT(jLft->valuedouble);
		    itv->iright = ALG_NINT(jRgt->valuedouble);
		  }
		  else
		  {
		    errNum = WLZ_ERR_DOMAIN_DATA;
		    break;
		  }
		}
		else
		{
		  errNum = WLZ_ERR_DOMAIN_DATA;
		  break;
		}
		jItv = jItv->next;
		++i;
	      }
	      if(errNum != WLZ_ERR_NONE)
	      {
		break;
	      }
	      jItvLn = jItvLn->next;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    iDom->freeptr = AlcFreeStackPush(iDom->freeptr, (void *)itvs, NULL);
	  }
	  else
	  {
	    AlcFree(itvs);
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(iDom);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Parse JSON encoded value table.
* \param	obj		Woolz object with domain for the object being
* 				parsed.
* \param	jVal		Pointer to cJSON value table.
*/
static WlzErrorNum		WlzEffParseJsnValueTable(
				  WlzObject *obj,
				  const cJSON *jVal)
{
  const cJSON 	*jBgd,
  		*jGType,
  		*jValues;
  WlzValues	val = {0};
  WlzPixelV	bgd = {0};
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  jGType = cJSON_GetObjectItemCaseSensitive(jVal, "grey_type");
  jBgd = cJSON_GetObjectItemCaseSensitive(jVal, "background");
  jValues = cJSON_GetObjectItemCaseSensitive(jVal, "values");
  if((jGType == NULL) || !cJSON_IsString(jGType) || !(jGType->valuestring) ||
     (jBgd == NULL) ||
     (jValues == NULL) || !cJSON_IsArray(jValues))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gType = WlzStringToGreyType(jGType->valuestring, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(gType == WLZ_GREY_RGBA)
    {
      bgd.type = WLZ_GREY_RGBA;
      errNum = WlzEffParseJsnRGBA(jBgd->valuestring, &(bgd.v.rgbv));
    }
    else
    {
      bgd.type = WLZ_GREY_DOUBLE;
      bgd.v.dbv = jBgd->valuedouble;
      (void )WlzValueConvertPixel(&bgd, bgd, gType);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObjectType vType;

    vType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, gType, NULL);
    val.v = WlzNewValueTb(obj, vType, bgd, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzLong	nVal;
    WlzIntervalWSpace     iWSp;
    WlzGreyWSpace         gWSp;

    nVal = cJSON_GetArraySize(jValues);
    obj->values = WlzAssignValues(val, NULL);
    if((errNum = WlzInitGreyScan(obj, &iWSp, &gWSp)) == WLZ_ERR_NONE)
    {
      int	i = 0;
      const cJSON *jV;

      jV = jValues->child;
      while((errNum == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
      {
	int	k;

	if(nVal - (i + iWSp.colrmn) < 0)
	{
	  errNum = WLZ_ERR_VALUES_DATA;
	}
	else
	{
	  double	v;

	  switch(gType)
	  {
	    case WLZ_GREY_INT:
	      for(k = 0; k < iWSp.colrmn; ++k)
	      {
	        v = ALG_CLAMP(jV->valuedouble, INT_MIN, INT_MAX);
		gWSp.u_grintptr.inp[k] = ALG_NINT(v);
		jV = jV->next;
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      for(k = 0; k < iWSp.colrmn; ++k)
	      {
		v = ALG_CLAMP(jV->valuedouble, SHRT_MIN, SHRT_MAX);
		gWSp.u_grintptr.shp[k] = (short )ALG_NINT(v);
		jV = jV->next;
	      }
	      break;
	    case WLZ_GREY_UBYTE:
	      for(k = 0; k < iWSp.colrmn; ++k)
	      {
		v = ALG_CLAMP(jV->valuedouble, 0, 255);
		gWSp.u_grintptr.ubp[k] = (WlzUByte )ALG_NINT(v);
		jV = jV->next;
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      for(k = 0; k < iWSp.colrmn; ++k)
	      {
		v = ALG_CLAMP(jV->valuedouble, -FLT_MAX, FLT_MAX);
		gWSp.u_grintptr.flp[k] = (float )v;
		jV = jV->next;
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      for(k = 0; k < iWSp.colrmn; ++k)
              {
		gWSp.u_grintptr.dbp[k] = jV->valuedouble;
		jV = jV->next;
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      for(k = 0; k < iWSp.colrmn; ++k)
	      {
		if((errNum = WlzEffParseJsnRGBA(jV->valuestring,
				    gWSp.u_grintptr.rgbp + k)) != WLZ_ERR_NONE)
		{
		  break;
		}
		jV = jV->next;
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	  i += k;
	}
      }
      (void )WlzEndGreyScan(&iWSp, &gWSp);
      if(errNum == WLZ_ERR_EOO)
      {
	errNum = WLZ_ERR_NONE;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Parses a Woolz property list from the given cJSON property
* 		list.
* \todo		Property lists are not yet written to or parsed from JSON
* 		encoded files.
* \param	obj
* \param	jPrp
*/
static WlzErrorNum		WlzEffParseJsnPropertyList(
				  WlzObject *obj,
				  const cJSON *jPrp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Parses a RGBA value from string.
* \param	str		Given value string.
* \param	rgba		Destination RGBA colour value pointer.
*/
static WlzErrorNum		WlzEffParseJsnRGBA(
				  const char *str,
				  WlzUInt *rgba)
{
  unsigned int 	vr,
		vg,
		vb,
		va;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(sscanf(str, "0x%02x%02x%02x%02x", &vr, &vg, &vb, &va) < 4)
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  else
  {
    WLZ_RGBA_RGBA_SET(*rgba, vr, vg, vb, va);
  }
  return(errNum);
}

