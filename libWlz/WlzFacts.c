#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzFacts_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzFacts.c
* \author       Bill Hill
* \date         june 2002
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
* \brief	Text description (facts) of Woolz objects.
* \ingroup	WlzDebug
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <Wlz.h>

typedef struct
{
  int		indent;
  int		verbose;
  FILE		*fFile;
  char		*fStr;
  int		fStrMax;
  int		fStrLen;
  int		fStrInc;
} WlzObjFactsData;

static WlzErrorNum		WlzObjFactsObject(
				  WlzObjFactsData *fData,
				  WlzObject *obj);
static WlzErrorNum		WlzObjFactsLinkcount(
				  WlzObjFactsData *fData,
				  int linkcount);
static WlzErrorNum		WlzObjFactsAppend(
				  WlzObjFactsData *fData,
				  const char *fmt,
				  ...);
static WlzErrorNum		WlzObjFactsBackground(
				  WlzObjFactsData *fData,
				   WlzPixelV bgdV);
static WlzErrorNum		WlzObjFactsIntervalDom(
				  WlzObjFactsData *fData,
				   WlzObject *obj,
				   WlzDomain dom);
static WlzErrorNum		WlzObjFactsValueTab(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzValues val);
static WlzErrorNum 		WlzObjFactsProperty(
				  WlzObjFactsData *fData,
				  WlzProperty	property);
static WlzErrorNum		WlzObjFactsPropList(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzPropertyList *pList);
static WlzErrorNum		WlzObjFactsPlaneDom(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzDomain dom);
static WlzErrorNum		WlzObjFactsVoxelTab(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzValues val);
static WlzErrorNum	        WlzObjFactsAffineTrans(
				  WlzObjFactsData *fData,
				   WlzAffineTransform *trans);
static WlzErrorNum		WlzObjFactsPolygon2D(
				  WlzObjFactsData *fData,
				  WlzPolygonDomain *poly);
static WlzErrorNum		WlzObjFactsBoundlist(
				  WlzObjFactsData *fData,
				  WlzBoundList *bList);
static WlzErrorNum		WlzObjFactsHistogram(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzHistogramDomain *hist);
static WlzErrorNum 		WlzObjFactsContour(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzContour *ctr);
static WlzErrorNum 		WlzObjFactsGMModel(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzGMModel *model);
static WlzErrorNum 		WlzObjFactsPixelV(
				  WlzObjFactsData *fData,
				  const char *prefix,
				  WlzPixelV pV);
static WlzErrorNum 		WlzObjFactsGreyV(
				  WlzObjFactsData *fData,
				  const char *prefix,
				  WlzGreyType gType,
				  WlzGreyV gV);
static WlzErrorNum	        WlzObjFactsMeshTrans(
				  WlzObjFactsData *fData,
				  WlzMeshTransform *mtrans);		
static WlzErrorNum 		WlzObjFactsCMesh2D(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzCMesh2D *mesh);
static WlzErrorNum 		WlzObjFactsCMesh3D(
				  WlzObjFactsData *fData,
				  WlzObject *obj,
				  WlzCMesh3D *mesh);
static WlzErrorNum 		WlzObjFactsCMeshTrans(
				  WlzObjFactsData *fData,
				  WlzCMeshTransform *cmt);

/*!
* \return	Woolz error code.
* \ingroup	WlzDebug
* \brief	The external facts interface function.
*		Creates the facts data structure and then produces a
*               text description of the given object.
* \param	obj			The given object.
* \param	factsFile		If non-NULL the text is written to the
* 					given file.
* \param	dstStr			If non-NULL the pointer is set to an
* 					allocated text string, which should be
*					free'd using AlcFree().
* \param	verbose			Verbose output if non-zero.
*/
WlzErrorNum	WlzObjectFacts(WlzObject *obj, FILE *factsFile, char **dstStr,
			       int verbose)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObjFactsData fData;
  const int	incTxtStep = 2048,     /* String increment */
  		initTxtLen = 4096;     /* Initial assumption for string size */

  fData.verbose = verbose;
  fData.indent = 0;
  fData.fFile = factsFile;
  fData.fStr = NULL;
  fData.fStrMax = 0;
  fData.fStrLen = 0;
  fData.fStrInc = incTxtStep;
  if(dstStr)
  {
    if((fData.fStr = (char *)AlcMalloc(initTxtLen * sizeof(char))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      fData.fStrMax = initTxtLen;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzObjFactsObject(&fData, obj);
  }
  if(dstStr)
  {
    *dstStr = fData.fStr;
  }
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzDebug
* \brief	Produces a text description of an object.
* \param	fData			Facts data structure.
* \param	obj			The given object.
*/
static WlzErrorNum WlzObjFactsObject(WlzObjFactsData *fData, WlzObject *obj)
{
  int		idx,
		count;
  WlzCompoundArray *objCA;
  WlzObject	*obj2D = NULL;
  WlzValues	val2D;
  WlzDomain	dom2D;
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  if(obj == NULL)
  {
    errNum = WlzObjFactsAppend(fData, "Object NULL.\n");
  }
  else if((tStr = WlzStringFromObjType(obj, NULL)) == NULL)
  {
    errNum = WlzObjFactsAppend(fData, "Object type invalid.\n");
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Object type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsLinkcount(fData, obj->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(obj->type)
      {
        case WLZ_2D_DOMAINOBJ:
	  errNum = WlzObjFactsIntervalDom(fData, obj, obj->domain);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsValueTab(fData, obj, obj->values);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_3D_DOMAINOBJ:
	  errNum = WlzObjFactsPlaneDom(fData, obj, obj->domain);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(obj->values.core)
	    {
	      errNum = WlzObjFactsVoxelTab(fData, obj, obj->values);
	    }
	    else
	    {
	      ++(fData->indent);
	      errNum = WlzObjFactsAppend(fData, "Values: NULL.\n");
	      --(fData->indent);
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  if((errNum == WLZ_ERR_NONE) && fData->verbose &&
	     (obj->domain.core != NULL))
	  {
	    idx = 0;
	    count = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
            ++(fData->indent);
	    fData->verbose = 0;
	    while((errNum == WLZ_ERR_NONE) && (idx < count))
	    {
	      errNum = WlzObjFactsAppend(fData, "Plane %d.\n",
	      				 obj->domain.p->plane1 + idx);
	      if(errNum == WLZ_ERR_NONE)
	      {
		dom2D = *(obj->domain.p->domains + idx);
		if(obj->values.core && (obj->values.vox->values + idx)->core)
		{
		  val2D = *(obj->values.vox->values + idx);
		}
		else
		{
		  val2D.core = NULL;
		}
		obj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom2D, val2D,
				    NULL, NULL, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		errNum = WlzObjFactsIntervalDom(fData, obj2D, obj2D->domain);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		errNum = WlzObjFactsValueTab(fData, obj2D, obj2D->values);
	      }
	      if(obj2D)
	      {
	        WlzFreeObj(obj2D);
		obj2D = NULL;
	      }
	      ++idx;
	    }
	    fData->verbose = 1;
            --(fData->indent);
	  }
	  break;
	case WLZ_TRANS_OBJ:
	  errNum = WlzObjFactsAffineTrans(fData, obj->domain.t);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsObject(fData, obj->values.obj);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_2D_POLYGON:
	  errNum = WlzObjFactsPolygon2D(fData, obj->domain.poly);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_BOUNDLIST:
	  errNum = WlzObjFactsBoundlist(fData, obj->domain.b);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_CONTOUR:
	  errNum = WlzObjFactsContour(fData, obj, obj->domain.ctr);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_CMESH_2D:
	  errNum = WlzObjFactsCMesh2D(fData, obj, obj->domain.cm2);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_CMESH_3D:
	  errNum = WlzObjFactsCMesh3D(fData, obj, obj->domain.cm3);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_HISTOGRAM:
	  errNum = WlzObjFactsHistogram(fData, obj, obj->domain.hist);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_AFFINE_TRANS:
	  errNum = WlzObjFactsAffineTrans(fData, obj->domain.t);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_3D_WARP_TRANS:   /* FALLTHROUGH */
	case WLZ_CONV_HULL:       /* FALLTHROUGH */
	case WLZ_3D_POLYGON:      /* FALLTHROUGH */
	case WLZ_RECTANGLE:       /* FALLTHROUGH */
	case WLZ_CONVOLVE_INT:    /* FALLTHROUGH */
	case WLZ_CONVOLVE_FLOAT:  /* FALLTHROUGH */
	case WLZ_WARP_TRANS:      /* FALLTHROUGH */
	case WLZ_FMATCHOBJ:       /* FALLTHROUGH */
	case WLZ_TEXT:            /* FALLTHROUGH */
	case WLZ_COMPOUND_ARR_1:  /* FALLTHROUGH */
	case WLZ_COMPOUND_ARR_2:
	  objCA = (WlzCompoundArray *)obj;
          tStr = WlzStringFromObjTypeValue(objCA->otype, &errNum);
	  if(tStr)
	  {
	    errNum = WlzObjFactsAppend(fData, "otype: %s.\n", tStr);
	  }
	  else
	  {
	    (void )WlzObjFactsAppend(fData, "otype: Unknown (%d).\n", tStr);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "n: %d.\n", objCA->n);
	  }
	  idx = 0;
	  while((errNum == WLZ_ERR_NONE) && (idx < objCA->n))
	  {
	    errNum = WlzObjFactsObject(fData, *(objCA->o + idx));
	    ++idx;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, objCA->plist);
	  }
	  break;
	case WLZ_MESH_TRANS:
	  errNum = WlzObjFactsMeshTrans(fData, obj->domain.mt);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_CMESH_TRANS:
	  errNum = WlzObjFactsCMeshTrans(fData, obj->domain.cmt);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsPropList(fData, obj, obj->plist);
	  }
	  break;
	case WLZ_EMPTY_OBJ:
	  /* no more facts available */
	  break;
	case WLZ_COMPOUND_LIST_1: /* FALLTHROUGH */
	case WLZ_COMPOUND_LIST_2: /* FALLTHROUGH */
	case WLZ_PROPERTY_OBJ:    /* FALLTHROUGH */
	default:
	  errNum = WlzObjFactsAppend(fData,
	  		"Facts not implemented for this object type.\n");
	  break;
      }
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzDebug
* \brief	Produces a text description of a linkcount.
* \param	fData			Facts data structure.
* \param	linkcount		The given linkcount.
*/
static WlzErrorNum WlzObjFactsLinkcount(WlzObjFactsData *fData, int linkcount)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", linkcount);
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzDebug
* \brief	Builds a text string using the given format and data, cf
* 		printf().
* \param	fData			Facts data structure.
* \param	fmt		Text format string.
*/
static WlzErrorNum WlzObjFactsAppend(WlzObjFactsData *fData, const char *fmt,
		    		     ...)
{
  va_list	args;
  int		indent,
  		len;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  va_start(args, fmt);
  if(fData->fFile)
  {
    indent = fData->indent;
    while(indent-- > 0)
    {
      (void )fprintf(fData->fFile, "%s", "  ");
    }
    (void )vfprintf(fData->fFile, fmt, args);
  }
  if(fData->fStr)
  {
    if((fData->fStrMax - fData->fStrLen) < fData->fStrInc)
    {
      fData->fStrMax += fData->fStrInc;
      fData->fStr = (char *)AlcRealloc(fData->fStr,
      				       fData->fStrMax * sizeof(char));
      if(fData->fStr == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      len = vsprintf(fData->fStr + fData->fStrLen, fmt, args);
      fData->fStrLen += len;
      *(fData->fStr + fData->fStrLen) = '\0';
    }
  }
  va_end(args);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzDebug
* \brief	Produces a text description of a background pixel value.
* \param	fData			Facts data structure.
* \param	bgdV			Background pixel value.
*/
static WlzErrorNum WlzObjFactsBackground(WlzObjFactsData *fData,
					 WlzPixelV bgdV)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzObjFactsPixelV(fData, "Background", bgdV);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzDebug
* \brief	Produces a text description of a pixel value.
* \param	fData			Facts data structure.
* \param	prefix			Prefix string.
* \param	pix			Pixel value.
*/
static WlzErrorNum WlzObjFactsPixelV(WlzObjFactsData *fData,
					const char *prefix, WlzPixelV pV)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*tStr;

  if(errNum == WLZ_ERR_NONE)
  {
    tStr = WlzStringFromGreyType(pV.type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzObjFactsAppend(fData, "%s type: %s.\n",
    			       prefix, tStr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzObjFactsGreyV(fData, prefix, pV.type, pV.v);
  }
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzDebug
* \brief	Produces a text description of a grey value.
* \param	fData			Facts data structure.
* \param	prefix			Prefix string.
* \param	gType			The grey type.
* \param	gV			Grey value.
*/
static WlzErrorNum WlzObjFactsGreyV(WlzObjFactsData *fData,
				    const char *prefix,
				    WlzGreyType gType, WlzGreyV gV)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  switch(gType)
  {
    case WLZ_GREY_LONG:
      errNum = WlzObjFactsAppend(fData, "%s value: %d.\n",
				 prefix, gV.lnv);
      break;
    case WLZ_GREY_INT:
      errNum = WlzObjFactsAppend(fData, "%s value: %d.\n",
				 prefix, gV.inv);
      break;
    case WLZ_GREY_SHORT:
      errNum = WlzObjFactsAppend(fData, "%s value: %d.\n",
				 prefix, gV.shv);
    case WLZ_GREY_UBYTE:
      errNum = WlzObjFactsAppend(fData, "%s value: %d.\n",
				 prefix, gV.ubv);
      break;
    case WLZ_GREY_FLOAT:
      errNum = WlzObjFactsAppend(fData, "%s value: %g.\n",
				 prefix, gV.flv);
      break;
    case WLZ_GREY_DOUBLE:
      errNum = WlzObjFactsAppend(fData, "%s value: %g.\n",
				 prefix, gV.dbv);
      break;
    case WLZ_GREY_RGBA:
      errNum = WlzObjFactsAppend(fData, "%s value: %d %d %d %d.\n",
				 prefix,
				 gV.rgbv & 0xff,
				 (gV.rgbv >> 8) & 0xff,
				 (gV.rgbv >> 16) & 0xff,
				 (gV.rgbv >> 24) & 0xff);
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzDebug
* \brief	Produces a text description of an interval domain.
* \param	fData			Facts data structure.
* \param	obj			Object to determine domain type.
* \param	dom			Domain union for interval domain.
*/
static WlzErrorNum WlzObjFactsIntervalDom(WlzObjFactsData *fData,
					  WlzObject *obj, WlzDomain dom)
{
  int		tI0,
		itvIdx,
	  	line;
  WlzIntervalLine *itvLn;
  WlzInterval	*itv;
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromObjDomainType(obj, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Domain NULL.\n");
      errNum = WLZ_ERR_NONE;
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Domain type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Domain type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", dom.i->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Line bounds: %d %d.\n",
		                 dom.i->line1, dom.i->lastln);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Column bounds: %d %d.\n",
		                 dom.i->kol1, dom.i->lastkl);
    }
    if((errNum == WLZ_ERR_NONE) && (dom.i->type == WLZ_INTERVALDOMAIN_INTVL))
    {
      tI0 = WlzIntervalCount(dom.i, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzObjFactsAppend(fData, "No. of intervals: %d\n", tI0);
      }
      if((fData->verbose) && (errNum == WLZ_ERR_NONE))
      {
	line = dom.i->line1;
        itvLn = dom.i->intvlines;
	while((line <= dom.i->lastln) && (errNum == WLZ_ERR_NONE))
	{
	  errNum = WlzObjFactsAppend(fData, "Line: %d.\n", line);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Intervals: %d.\n",
	    			       itvLn->nintvs);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    ++(fData->indent);
	    itvIdx = 0;
	    itv = itvLn->intvs;
	    while((itvIdx < itvLn->nintvs) && (errNum == WLZ_ERR_NONE))
	    {
	      errNum = WlzObjFactsAppend(fData, "%8d %8d %8d\n",
	      				 itvIdx, itv->ileft, itv->iright);
	      ++itv;
	      ++itvIdx;
	    }
	    --(fData->indent);
	  }
	  ++line;
	  ++itvLn;
	}
      }
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a value table.
* \param	fData			Facts data structure.
* \param	obj			Object to determine value type.
* \param	val			values union for value table.
*/
static WlzErrorNum WlzObjFactsValueTab(WlzObjFactsData *fData,
				       WlzObject *obj, WlzValues val)
{
  const char	*tStr;
  WlzObjectType	valType;
  WlzGreyType	gType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromObjValuesType(obj, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_VALUES_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Values NULL.\n");
      errNum = WLZ_ERR_NONE;
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Values type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Values type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n",
      				 val.core->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      valType = WlzGreyTableTypeToTableType(val.core->type, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      gType = WlzGreyTableTypeToGreyType(val.core->type, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      tStr = WlzStringFromGreyType(gType, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Grey type: %s.\n", tStr);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(valType)
      {
        case WLZ_GREY_TAB_RAGR:
          errNum = WlzObjFactsBackground(fData, val.v->bckgrnd);
	  break;
	case WLZ_GREY_TAB_RECT:
          errNum = WlzObjFactsBackground(fData, val.r->bckgrnd);
	  break;
	case WLZ_GREY_TAB_INTL:
          errNum = WlzObjFactsBackground(fData, val.i->bckgrnd);
	  break;
        default:
	  errNum = WLZ_ERR_VALUES_TYPE;
	  break;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(valType)
      {
        case WLZ_GREY_TAB_RAGR:
	  errNum = WlzObjFactsAppend(fData, "Values line bounds: %d %d\n",
	  			     val.v->line1, val.v->lastln);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Values column bounds: %d %d\n",
	    			       val.v->kol1,
				       val.v->kol1 + val.v->width - 1);
	  }
	  break;
	case WLZ_GREY_TAB_RECT:
	  errNum = WlzObjFactsAppend(fData, "Values line bounds: %d %d\n",
	  			     val.r->line1, val.r->lastln);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Values column bounds: %d %d\n",
	    			       val.r->kol1,
				       val.r->kol1 + val.r->width - 1);
	  }
	  break;
	case WLZ_GREY_TAB_INTL:
	  errNum = WlzObjFactsAppend(fData, "Values line bounds: %d %d\n",
	  			     val.i->line1, val.i->lastln);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Values column bounds: %d %d\n",
	    			       val.i->kol1,
				       val.i->kol1 + val.i->width - 1);
	  }
	  break;
        default:
	  errNum = WLZ_ERR_VALUES_TYPE;
	  break;
      }
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a property.
* \param	fData			Facts data structure.
* \param	prop			Given property.
*/
static WlzErrorNum WlzObjFactsProperty(WlzObjFactsData *fData,
					WlzProperty prop)
{
  const char	*pStr;
  WlzEMAPProperty *eProp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* No indent because already done by WlzObjFactsPropertyList(). */
  if(prop.core == NULL)
  {
    (void )WlzObjFactsAppend(fData, "Property NULL.\n");
  }
  else
  {
    pStr = WlzStringFromPropertyType(prop, NULL);
    if(pStr)
    {
      (void )WlzObjFactsAppend(fData, "Property type:%s\n", pStr);
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Property type: Unknown (%d).\n",
      			       (int )(prop.core->type));
    }
    ++(fData->indent);
    errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n",
    			       prop.core->linkcount);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(prop.core->type)
      {
	case WLZ_PROPERTY_SIMPLE:
	  errNum = WlzObjFactsAppend(fData, "Size: %d.\n",
				     prop.simple->size);
	  break;
	case WLZ_PROPERTY_EMAP:
	  eProp = prop.emap;
	  pStr = WlzStringFromEMAPPropertyType(eProp, NULL);
	  if(pStr)
	  {
	    errNum = WlzObjFactsAppend(fData,
	    			"EMAP property type: %s.\n", pStr);
	  }
	  else
	  {
	    errNum = WlzObjFactsAppend(fData,
	    			"EMAP property type: Unknown (%d).\n",
	    		      	(int )eProp->type);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Model UID: %s\n",
	      			eProp->modelUID);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Anatomy UID: %s\n",
	      			eProp->anatomyUID);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Target UID: %s\n",
	      			eProp->targetUID);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Target Version: %s\n",
	      			eProp->targetVersion);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	      errNum = WlzObjFactsAppend(fData, "Stage: %s\n",
	      			eProp->stage);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	      errNum = WlzObjFactsAppend(fData, "Sub-stage: %s\n",
	      			eProp->subStage);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Model name: %s\n",
	      			eProp->modelName);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Model version: %s\n",
	      			eProp->version);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Filename: %s\n",
	      			eProp->fileName?eProp->fileName:"NULL");
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Creation time: %s",
	      			ctime((time_t *) &(eProp->creationTime)));
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Creation author: %s\n",
	      			eProp->creationAuthor);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Creation machine name: %s\n",
	      			eProp->creationMachineName);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Modification time: %s",
	      			ctime((time_t *) &(eProp->modificationTime)));
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Modification author: %s\n",
	      			eProp->modificationAuthor);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzObjFactsAppend(fData, "Comment: %s\n",
	      			eProp->comment?eProp->comment:"NULL");
	  }
	  break;
	case WLZ_PROPERTY_NAME:
	  errNum = WlzObjFactsAppend(fData, "Name: %s.\n",
	  		    	prop.name->name? prop.name->name: "NULL");
	  break;
	case WLZ_PROPERTY_GREY:
	  errNum = WlzObjFactsAppend(fData, "Name: %s.\n",
	  		    	prop.greyV->name? prop.greyV->name: "NULL");
	  if(errNum == WLZ_ERR_NONE)
	  {
  	    errNum = WlzObjFactsPixelV(fData, "Grey", prop.greyV->value);
	  }
	  break;
	default:
	  break;
      }
    }
    --(fData->indent);
  }
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a property list.
* \param	fData			Facts data structure.
* \param	dummy			Object, unused.
* \param	pList			Given property list.
*/
static WlzErrorNum WlzObjFactsPropList(WlzObjFactsData *fData,
				       WlzObject *dummy,
				       WlzPropertyList *pList)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  if((pList == NULL) || (pList->list == NULL))
  {
    (void )WlzObjFactsAppend(fData, "Property list NULL.\n");
  }
  else
  {
    AlcDLPItem	*item;
    (void )WlzObjFactsAppend(fData, "Property list:\n");
    item = pList->list->head;
    do {
      WlzProperty	property;
      if( item ){
	if( item->entry ){
	  property.core = (WlzCoreProperty *) item->entry;
	  errNum = WlzObjFactsProperty(fData, property);
	}
	item = item->next;
      }
    } while(item != pList->list->head);
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a plane domain.
* \param	fData			Facts data structure.
* \param	obj			Object to determine domain type.
* \param	dom			Domain union for plane domain.
*/
static WlzErrorNum WlzObjFactsPlaneDom(WlzObjFactsData *fData,
				       WlzObject *obj, WlzDomain dom)
{
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromObjDomainType(obj, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Domain NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Domain type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Domain type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", dom.i->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Domain plane bounds: %d %d\n",
				 dom.p->plane1, dom.p->lastpl);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Domain line bounds: %d %d\n",
				 dom.p->line1, dom.p->lastln);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Domain column bounds: %d %d\n",
				 dom.p->kol1, dom.p->lastkl);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "VoxelSize: %g %g %g\n",
				 dom.p->voxel_size[0],
				 dom.p->voxel_size[1],
				 dom.p->voxel_size[2]);
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a voxel table.
* \param	fData			Facts data structure.
* \param	obj			Object to determine values type.
* \param	val			Values union for voxel table.
*/
static WlzErrorNum WlzObjFactsVoxelTab(WlzObjFactsData *fData,
				       WlzObject *obj, WlzValues val)
{
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromObjValuesType(obj, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Values NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Values type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Values type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n",
      				 val.core->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsBackground(fData, val.vox->bckgrnd);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Values plane bounds: %d %d\n",
				 val.vox->plane1, val.vox->lastpl);
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of an affine transform.
* \param	fData			Facts data structure.
* \param	trans			Given affine transform.
*/
static WlzErrorNum WlzObjFactsAffineTrans(WlzObjFactsData *fData,
				          WlzAffineTransform *trans)
{
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromTransformType(trans->type, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Transform NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Transform type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Transform type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", trans->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(WlzAffineTransformDimension(trans, NULL))
      {
        case 2:
	  errNum = WlzObjFactsAppend(fData, "mat: %-10g %-10g %-10g\n",
				     trans->mat[0][0], trans->mat[0][1],
				     trans->mat[0][2]);
	  errNum = WlzObjFactsAppend(fData, "     %-10g %-10g %-10g\n",
				     trans->mat[1][0], trans->mat[1][1],
				     trans->mat[1][2]);
	  errNum = WlzObjFactsAppend(fData, "     %-10g %-10g %-10g\n",
				     trans->mat[2][0], trans->mat[2][1],
				     trans->mat[2][2]);
	  break;
	case 3:
	  errNum = WlzObjFactsAppend(fData, "mat: %-10g %-10g %-10g %-10g\n",
			  trans->mat[0][0], trans->mat[0][1],
			  trans->mat[0][2], trans->mat[0][3]);
	  errNum = WlzObjFactsAppend(fData, "     %-10g %-10g %-10g %-10g\n",
			  trans->mat[1][0], trans->mat[1][1],
			  trans->mat[1][2], trans->mat[1][3]);
	  errNum = WlzObjFactsAppend(fData, "     %-10g %-10g %-10g %-10g\n",
			  trans->mat[2][0], trans->mat[2][1],
			  trans->mat[2][2], trans->mat[2][3]);
	  errNum = WlzObjFactsAppend(fData, "     %-10g %-10g %-10g %-10g\n",
			  trans->mat[3][0], trans->mat[3][1],
			  trans->mat[3][2], trans->mat[3][3]);
	  break;
        default:
	  errNum = WLZ_ERR_TRANSFORM_TYPE;
	  break;
      }
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of an 2D polygon domain.
* \param	fData			Facts data structure.
* \param	poly			Given polygon domain.
*/
static WlzErrorNum WlzObjFactsPolygon2D(WlzObjFactsData *fData,
					WlzPolygonDomain *poly)
{
  int		idx;
  WlzIVertex2	*tIVx0;
  WlzFVertex2	*tFVx0;
  WlzDVertex2	*tDVx0;
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  if(poly == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(poly->type)
    {
      case WLZ_POLYGON_INT:
	tStr = "WLZ_POLYGON_INT";
        break;
      case WLZ_POLYGON_FLOAT:
	tStr = "WLZ_POLYGON_FLOAT";
        break;
      case WLZ_POLYGON_DOUBLE:
	tStr = "WLZ_POLYGON_DOUBLE";
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Domain NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Domain type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Domain type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", poly->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "No. verticies: %d.\n",
      			         poly->nvertices);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Max verticies: %d.\n",
      			         poly->maxvertices);
    }
    if(fData->verbose && (errNum == WLZ_ERR_NONE))
    {
      ++(fData->indent);
      switch(poly->type)
      {
        case WLZ_POLYGON_INT:
	  idx = 0;
	  tIVx0 = poly->vtx;
	  while((idx < poly->nvertices) && (errNum == WLZ_ERR_NONE))
	  {
	    errNum = WlzObjFactsAppend(fData, "%8d %8d %8d.\n",
	    			       idx, tIVx0->vtX, tIVx0->vtY);
	    ++idx;
	    ++tIVx0;
	  }
	  break;
        case WLZ_POLYGON_FLOAT:
	  idx = 0;
	  tFVx0 = (WlzFVertex2 *)(poly->vtx);
	  while((idx < poly->nvertices) && (errNum == WLZ_ERR_NONE))
	  {
	    errNum = WlzObjFactsAppend(fData, "%8d %8g %8g.\n",
	    			       idx, tFVx0->vtX, tFVx0->vtY);
	    ++idx;
	    ++tFVx0;
	  }
	  break;
        case WLZ_POLYGON_DOUBLE:
	  idx = 0;
	  tDVx0 = (WlzDVertex2 *)(poly->vtx);
	  while((idx < poly->nvertices) && (errNum == WLZ_ERR_NONE))
	  {
	    errNum = WlzObjFactsAppend(fData, "%8d %8g %8g.\n",
	    			       idx, tDVx0->vtX, tDVx0->vtY);
	    ++idx;
	    ++tDVx0;
	  }
	  break;
        default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
      }
      --(fData->indent);
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a boundary list.
* \param	fData			Facts data structure.
* \param	bList			Given boundary list.
*/
static WlzErrorNum WlzObjFactsBoundlist(WlzObjFactsData *fData,
					WlzBoundList *bList)
{
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  if(bList == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(bList->type)
    {
      case WLZ_BOUNDLIST_PIECE:
	tStr = "WLZ_BOUNDLIST_PIECE";
        break;
      case WLZ_BOUNDLIST_HOLE:
	tStr = "WLZ_BOUNDLIST_HOLE";
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Domain NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Domain type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Domain type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", bList->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "poly:\n");
    }
    if(fData->verbose && (errNum == WLZ_ERR_NONE))
    {
      if(bList->poly)
      {
        errNum = WlzObjFactsPolygon2D(fData, bList->poly);
      }
    }
    ++(fData->indent);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "down:\n");
    }
    if(fData->verbose && (errNum == WLZ_ERR_NONE))
    {
      if(bList->down)
      {
        errNum = WlzObjFactsBoundlist(fData, bList->down);
      }
    }
    --(fData->indent);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "next:\n");
    }
    if(fData->verbose && (errNum == WLZ_ERR_NONE))
    {
      if(bList->next)
      {
        errNum = WlzObjFactsBoundlist(fData, bList->next);
      }
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a histogram domain.
* \param	fData			Facts data structure.
* \param	obj			Object to determine domain type.
* \param	hist			Given histogram domain.
*/
static WlzErrorNum WlzObjFactsHistogram(WlzObjFactsData *fData,
				        WlzObject *obj,
					WlzHistogramDomain *hist)
{
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromObjDomainType(obj, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Domain NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Domain type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Domain type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", hist->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "No. bins: %d\n",
				 hist->nBins);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Origin: %g\n",
				 hist->origin);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Bin size: %g\n",
				 hist->binSize);
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a contour domain.
* \param	fData			Facts data structure.
* \param	obj			Object for type, not used.
* \param	ctr			Given contour domain.
*/
static WlzErrorNum WlzObjFactsContour(WlzObjFactsData *fData,
				      WlzObject *obj,
				      WlzContour *ctr)
{
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromObjDomainType(obj, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Domain NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Domain type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Domain type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", ctr->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsGMModel(fData, obj, ctr->model);
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a geometric model.
* \param	fData			Facts data structure.
* \param	dummy			Object, not used.
* \param	model			Given geometric model.
*/
static WlzErrorNum WlzObjFactsGMModel(WlzObjFactsData *fData,
				      WlzObject *dummy,
				      WlzGMModel *model)
{
  int		idS,
  		dim;
  WlzGMShell	*s;
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  if(model == NULL)
  {
    (void )WlzObjFactsAppend(fData, "Model NULL.\n");
  }
  else
  {
    tStr = WlzStringFromGMModelType(model->type, &errNum);
    if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
    {
      (void )WlzObjFactsAppend(fData, "Model type invalid.\n");
    }
    else
    {
      switch(model->type)
      {
        case WLZ_GMMOD_2I:
	case WLZ_GMMOD_2D:
	case WLZ_GMMOD_2N:
	  dim = 2;
	  break;
	case WLZ_GMMOD_3I:
	case WLZ_GMMOD_3D:
	case WLZ_GMMOD_3N:
	  dim = 3;
	  break;
      }
      errNum = WlzObjFactsAppend(fData, "Model type: %s.\n", tStr);
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzObjFactsAppend(fData,
				   "Linkcount: %d.\n", model->linkcount);
      }
      if((errNum == WLZ_ERR_NONE) && fData->verbose)
      {
        errNum = WlzObjFactsAppend(fData,
				   "nVertex: %d\n",
				   model->res.vertex.numElm);
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nVertexT: %d\n",
				     model->res.vertexT.numElm);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nVertexG: %d\n",
				     model->res.vertexG.numElm);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nDiskT: %d\n",
				     model->res.diskT.numElm);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nEdge: %d\n",
				     model->res.edge.numElm);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nEdgeT: %d\n",
				     model->res.edgeT.numElm);
	}
        if((errNum == WLZ_ERR_NONE) && (dim == 3))
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nFace: %d\n",
				     model->res.face.numElm);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nLoopT: %d\n",
				     model->res.loopT.numElm);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nShell: %d\n",
				     model->res.shell.numElm);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzObjFactsAppend(fData,
	  			     "nShellG: %d\n",
				     model->res.shellG.numElm);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  idS = 0;
	  ++(fData->indent);
	  while((errNum == WLZ_ERR_NONE) && (idS < model->res.shell.numIdx))
	  {
	    if(((s = (WlzGMShell *)
	             AlcVectorItemGet(model->res.shell.vec, idS)) != NULL) &&
	       (s->idx >= 0))
            {
	      switch(s->geo.core->type)
	      {
	        case WLZ_GMELM_SHELL_G2I:
		  errNum = WlzObjFactsAppend(fData,
		                             "% 6d (%d, %d, %d, %d)\n",
					     idS,
					     s->geo.sg2I->bBox.xMin,
					     s->geo.sg2I->bBox.yMin,
					     s->geo.sg2I->bBox.xMax,
					     s->geo.sg2I->bBox.yMax);
		  break;
	        case WLZ_GMELM_SHELL_G2D:
		  errNum = WlzObjFactsAppend(fData,
		                             "% 6d (%g, %g, %g, %g)\n",
					     idS,
					     s->geo.sg2D->bBox.xMin,
					     s->geo.sg2D->bBox.yMin,
					     s->geo.sg2D->bBox.xMax,
					     s->geo.sg2D->bBox.yMax);
		  break;
	        case WLZ_GMELM_SHELL_G3I:
		  errNum = WlzObjFactsAppend(fData,
		                             "% 6d (%d, %d, %d, %d, %d, %d)\n",
					     idS,
					     s->geo.sg3I->bBox.xMin,
					     s->geo.sg3I->bBox.yMin,
					     s->geo.sg3I->bBox.zMin,
					     s->geo.sg3I->bBox.xMax,
					     s->geo.sg3I->bBox.yMax,
					     s->geo.sg3I->bBox.zMax);
		  break;
	        case WLZ_GMELM_SHELL_G3D:
		  errNum = WlzObjFactsAppend(fData,
		                             "% 6d (%g, %g, %g, %g, %g, %g)\n",
					     idS,
					     s->geo.sg3D->bBox.xMin,
					     s->geo.sg3D->bBox.yMin,
					     s->geo.sg3D->bBox.zMin,
					     s->geo.sg3D->bBox.xMax,
					     s->geo.sg3D->bBox.yMax,
					     s->geo.sg3D->bBox.zMax);
		  break;
	        default:
		  break;
	      }
	    }
	    ++idS;
	  }
	  --(fData->indent);
	}
      }
    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of an mesh transform.
* \param	fData			Facts data structure.
* \param	trans			Given mesh transform.
*/
static WlzErrorNum WlzObjFactsMeshTrans(WlzObjFactsData *fData,
				          WlzMeshTransform *trans)
{
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromTransformType(trans->type, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Transform NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Transform type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Transform type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", trans->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Number of Elements: %d.\n", trans->nElem);
      errNum = WlzObjFactsAppend(fData, "Number of Nodes: %d.\n", trans->nNodes);

    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a 2D constrained mesh.
* \param	fData			Facts data structure.
* \param	obj			Object for type.
* \param	cmt			Given constrained mesh transform.
*/
static WlzErrorNum WlzObjFactsCMesh2D(WlzObjFactsData *fData,
				      WlzObject *obj,
				      WlzCMesh2D *mesh)
{
  int		idx,
  		nElm,
  		nNod;
  const char	*tStr;
  WlzCMeshNod2D	*nod;
  WlzCMeshElm2D *elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromObjDomainType(obj, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Domain NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Domain type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Domain type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", mesh->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nElm = mesh->res.elm.numEnt;
      nNod = mesh->res.nod.numEnt;
      errNum = WlzObjFactsAppend(fData, "Number of Elements: %d.\n", nElm);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzObjFactsAppend(fData, "Number of Nodes: %d.\n", nNod);
      }
    }
    if((errNum == WLZ_ERR_NONE) && (fData->verbose != 0))
    {
      errNum = WlzObjFactsAppend(fData, "Nodes:\n");
      if(errNum == WLZ_ERR_NONE)
      {
	idx = 0;
        ++(fData->indent);
	while((errNum == WLZ_ERR_NONE) && (idx < mesh->res.nod.maxEnt))
	{
	  nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idx);
	  if(nod->idx >= 0)
	  {
	    errNum = WlzObjFactsAppend(fData, "% 8d % 8g % 8g\n",
				       nod->idx, nod->pos.vtX, nod->pos.vtY);
	  }
	  ++idx;
	}
        --(fData->indent);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzObjFactsAppend(fData, "Elements:\n");
      }
      if(errNum == WLZ_ERR_NONE)
      {
	idx = 0;
        ++(fData->indent);
	while((errNum == WLZ_ERR_NONE) && (idx < mesh->res.elm.maxEnt))
	{
	  elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idx);
	  if(elm->idx >= 0)
	  {
	    errNum = WlzObjFactsAppend(fData, "% 8d % 8d % 8d % 8d\n",
				       elm->idx,
				       elm->edu[0].nod->idx,
				       elm->edu[1].nod->idx,
				       elm->edu[2].nod->idx);
	  }
	  ++idx;
	}
        --(fData->indent);
      }

    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a 3D constrained mesh.
* \param	fData			Facts data structure.
* \param	obj			Object for type.
* \param	cmt			Given constrained mesh transform.
*/
static WlzErrorNum WlzObjFactsCMesh3D(WlzObjFactsData *fData,
				      WlzObject *obj,
				      WlzCMesh3D *mesh)
{
  int		idx,
  		nElm,
  		nNod;
  const char	*tStr;
  WlzCMeshNod3D	*nod;
  WlzCMeshElm3D *elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromObjDomainType(obj, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Domain NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Domain type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Domain type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", mesh->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nElm = mesh->res.elm.numEnt;
      nNod = mesh->res.nod.numEnt;
      errNum = WlzObjFactsAppend(fData, "Number of Elements: %d.\n", nElm);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzObjFactsAppend(fData, "Number of Nodes: %d.\n", nNod);
      }
    }
    if((errNum == WLZ_ERR_NONE) && (fData->verbose != 0))
    {
      errNum = WlzObjFactsAppend(fData, "Nodes:\n");
      if(errNum == WLZ_ERR_NONE)
      {
	idx = 0;
        ++(fData->indent);
	while((errNum == WLZ_ERR_NONE) && (idx < mesh->res.nod.maxEnt))
	{
	  nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idx);
	  if(nod->idx >= 0)
	  {
	    errNum = WlzObjFactsAppend(fData, "% 8d % 8g % 8g % 8g\n",
				       nod->idx, nod->pos.vtX,
				       nod->pos.vtY, nod->pos.vtZ);
	  }
	  ++idx;
	}
        --(fData->indent);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzObjFactsAppend(fData, "Elements:\n");
      }
      if(errNum == WLZ_ERR_NONE)
      {
	idx = 0;
        ++(fData->indent);
	while((errNum == WLZ_ERR_NONE) && (idx < mesh->res.elm.maxEnt))
	{
	  elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idx);
	  if(elm->idx >= 0)
	  {
	    errNum = WlzObjFactsAppend(fData, "% 8d % 8d % 8d % 8d % 8d\n",
				       elm->idx,
				       elm->face[0].edu[0].nod->idx,
				       elm->face[0].edu[1].nod->idx,
				       elm->face[0].edu[2].nod->idx,
				       elm->face[1].edu[1].nod->idx);
	  }
	  ++idx;
	}
        --(fData->indent);
      }

    }
  }
  --(fData->indent);
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup      WlzDebug
* \brief	Produces a text description of a constrained mesh transform.
* \param	fData			Facts data structure.
* \param	cmt			Given constrained mesh transform.
*/
static WlzErrorNum WlzObjFactsCMeshTrans(WlzObjFactsData *fData,
				         WlzCMeshTransform *cmt)
{
  int		nElm = 0,
  		nNod = 0;
  const char	*tStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  tStr = WlzStringFromTransformType(cmt->type, &errNum);
  if((tStr == NULL) || (errNum != WLZ_ERR_NONE))
  {
    if(errNum == WLZ_ERR_DOMAIN_NULL)
    {
      (void )WlzObjFactsAppend(fData, "Transform NULL.\n");
    }
    else
    {
      (void )WlzObjFactsAppend(fData, "Transform type invalid.\n");
    }
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Transform type: %s.\n", tStr);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", cmt->linkcount);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(cmt->type)
      {
        case WLZ_TRANSFORM_2D_CMESH:
	  nElm = cmt->mesh.m2->res.elm.numEnt;
	  nNod = cmt->mesh.m2->res.nod.numEnt;
	  break;
	case WLZ_TRANSFORM_3D_CMESH:
	  nElm = cmt->mesh.m3->res.elm.numEnt;
	  nNod = cmt->mesh.m3->res.nod.numEnt;
	  break;
	default:
	  errNum = WLZ_ERR_TRANSFORM_TYPE;
	  break;
      }
      errNum = WlzObjFactsAppend(fData, "Number of Elements: %d.\n", nElm);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzObjFactsAppend(fData, "Number of Nodes: %d.\n", nNod);
      }
    }
  }
  --(fData->indent);
  return(errNum);
}
