#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzFacts.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Text description (facts) of a Woolz object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
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

static WlzErrorNum WlzObjFactsObject(WlzObjFactsData *fData, WlzObject *obj),
		WlzObjFactsLinkcount(WlzObjFactsData *fData, int linkcount),
		WlzObjFactsAppend(WlzObjFactsData *fData,
				  const char *fmt, ...),
		WlzObjFactsBackground(WlzObjFactsData *fData,
				     WlzPixelV bgdV),
		WlzObjFactsIntervalDom(WlzObjFactsData *fData,
				       WlzObject *obj, WlzDomain dom),
		WlzObjFactsValueTab(WlzObjFactsData *fData,
				    WlzObject *obj, WlzValues val),
		WlzObjFactsPropList(WlzObjFactsData *fData,
				    WlzObject *obj,
				    WlzSimpleProperty *pList),
		WlzObjFactsPlaneDom(WlzObjFactsData *fData,
				    WlzObject *obj, WlzDomain dom),
		WlzObjFactsVoxelTab(WlzObjFactsData *fData,
				    WlzObject *obj, WlzValues val),
	        WlzObjFactsAffineTrans(WlzObjFactsData *fData,
				       WlzAffineTransform *trans),
		WlzObjFactsPolygon2D(WlzObjFactsData *fData,
				     WlzPolygonDomain *poly),
		WlzObjFactsBoundlist(WlzObjFactsData *fData,
				     WlzBoundList *bList),
		WlzObjFactsHistogram(WlzObjFactsData *fData,
				     WlzObject *obj,
				     WlzHistogramDomain *hist);
		

/************************************************************************
* Function:	WlzObjectFacts						*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	The external facts interface function.			*
*		Creates the facts data structure and then produces a	*
*		text description of the given object.			*
* Global refs:	-							*
* Parameters:	WlzObject *obj:		The given object.		*
*		FILE *factsFile:	If non-NULL the text is		*
*					written to the file stream.	*
*		char **dstStr:		If *dstStr is non-NULL the	*
*					pointer is set to an allocated	*
*					text string, which should be	*
*					free'd using AlcFree().		*
*		int verbose:		Verbose output if non-zero.	*
************************************************************************/
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

/************************************************************************
* Function:	WlzObjFactsObject					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of an object.		*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzObject *obj:		The given object.		*
************************************************************************/
static WlzErrorNum WlzObjFactsObject(WlzObjFactsData *fData, WlzObject *obj)
{
  int		pIdx,
		pCount;
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
	    pIdx = 0;
	    pCount = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
            ++(fData->indent);
	    while((errNum == WLZ_ERR_NONE) && (pIdx < pCount))
	    {
	      errNum = WlzObjFactsAppend(fData, "Plane %d.\n",
	      				 obj->domain.p->plane1 + pIdx);
	      if(errNum == WLZ_ERR_NONE)
	      {
		dom2D = *(obj->domain.p->domains + pIdx);
		if(obj->values.core && (obj->values.vox->values + pIdx)->core)
		{
		  val2D = *(obj->values.vox->values + pIdx);
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
	      ++pIdx;
	    }
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
	  break;
	case WLZ_BOUNDLIST:
	  errNum = WlzObjFactsBoundlist(fData, obj->domain.b);
	  break;
	case WLZ_HISTOGRAM:
	  errNum = WlzObjFactsHistogram(fData, obj, obj->domain.hist);
	  break;
	case WLZ_AFFINE_TRANS:
	  errNum = WlzObjFactsAffineTrans(fData, obj->domain.t);
	  break;
	case WLZ_3D_WARP_TRANS:   /* FALLTHROUGH */
	case WLZ_CONV_HULL:       /* FALLTHROUGH */
	case WLZ_3D_POLYGON:      /* FALLTHROUGH */
	case WLZ_RECTANGLE:       /* FALLTHROUGH */
	case WLZ_VECTOR_INT:      /* FALLTHROUGH */
	case WLZ_VECTOR_FLOAT:    /* FALLTHROUGH */
	case WLZ_POINT_INT:       /* FALLTHROUGH */
	case WLZ_POINT_FLOAT:     /* FALLTHROUGH */
	case WLZ_CONVOLVE_INT:    /* FALLTHROUGH */
	case WLZ_CONVOLVE_FLOAT:  /* FALLTHROUGH */
	case WLZ_WARP_TRANS:      /* FALLTHROUGH */
	case WLZ_FMATCHOBJ:       /* FALLTHROUGH */
	case WLZ_TEXT:            /* FALLTHROUGH */
	case WLZ_COMPOUND_ARR_1:  /* FALLTHROUGH */
	case WLZ_COMPOUND_ARR_2:  /* FALLTHROUGH */
	case WLZ_COMPOUND_LIST_1: /* FALLTHROUGH */
	case WLZ_COMPOUND_LIST_2: /* FALLTHROUGH */
	case WLZ_PROPERTY_OBJ:    /* FALLTHROUGH */
	case WLZ_EMPTY_OBJ:       /* FALLTHROUGH */
	default:
	  break;
      }
    }
  }
  --(fData->indent);
  return(errNum);
}

/************************************************************************
* Function:	WlzObjFactsLinkcount					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of a linkcount.		*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		int linkcount:		The given linkcount.		*
************************************************************************/
static WlzErrorNum WlzObjFactsLinkcount(WlzObjFactsData *fData, int linkcount)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n", linkcount);
  --(fData->indent);
  return(errNum);
}

/************************************************************************
* Function:	WlzObjFactsAppend					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Builds a text string using the given format and data,	*
*		cf printf.						*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		int linkcount:		The given linkcount.		*
************************************************************************/
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

/************************************************************************
* Function:	WlzObjFactsBackground					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of a background pixel	*
*		value.							*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzPixelV bgdV:		Background pixel value.		*
************************************************************************/
static WlzErrorNum WlzObjFactsBackground(WlzObjFactsData *fData,
					 WlzPixelV bgdV)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*tStr;

  if(errNum == WLZ_ERR_NONE)
  {
    tStr = WlzStringFromGreyType(bgdV.type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzObjFactsAppend(fData, "Background grey type: %s.\n", tStr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(bgdV.type)
    {
      case WLZ_GREY_INT:
	errNum = WlzObjFactsAppend(fData, "Background value: %d.\n",
				   bgdV.v.inv);
        break;
      case WLZ_GREY_SHORT:
	errNum = WlzObjFactsAppend(fData, "Background value: %d.\n",
				   bgdV.v.shv);
      case WLZ_GREY_UBYTE:
	errNum = WlzObjFactsAppend(fData, "Background value: %d.\n",
				   bgdV.v.ubv);
        break;
      case WLZ_GREY_FLOAT:
	errNum = WlzObjFactsAppend(fData, "Background value: %g.\n",
				   bgdV.v.flv);
        break;
      case WLZ_GREY_DOUBLE:
	errNum = WlzObjFactsAppend(fData, "Background value: %g.\n",
				   bgdV.v.dbv);
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzObjFactsIntervalDom					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of an interval domain.	*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzObject *obj:		Object to determine domain	*
*					type.				*
*		WlzDomain dom:		Domain union for interval	*
*					domain.				*
************************************************************************/
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

/************************************************************************
* Function:	WlzObjFactsValueTab					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of a value table.		*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzObject *obj:		Object to determine value type.	*
*		WlzValues val:		values union for value table.	*
************************************************************************/
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
      }
    }
  }
  --(fData->indent);
  return(errNum);
}

/************************************************************************
* Function:	WlzObjFactsPropList					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of a property list.		*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzObject *obj:		Object to determine property	*
*					list type.			*
*		WlzSimpleProperty *pList: Given property list.		*
************************************************************************/
static WlzErrorNum WlzObjFactsPropList(WlzObjFactsData *fData,
				       WlzObject *obj,
				       WlzSimpleProperty *pList)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ++(fData->indent);
  if(pList == NULL)
  {
    (void )WlzObjFactsAppend(fData, "Property list NULL.\n");
  }
  else
  {
    errNum = WlzObjFactsAppend(fData, "Linkcount: %d.\n",
			       pList->linkcount);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "Size: %d.\n",
				 pList->size);
    }
  }
  --(fData->indent);
  return(errNum);
}

/************************************************************************
* Function:	WlzObjFactsPlaneDom					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of a plane domain.		*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzObject *obj:		Object to determine domain	*
*					type.				*
*		WlzDomain dom:		Domain union for plane domain.	*
************************************************************************/
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

/************************************************************************
* Function:	WlzObjFactsVoxelTab					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of a voxel table.		*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzObject *obj:		Object to determine values	*
*					type.				*
*		WlzValues val:		Values union for voxel table.	*
************************************************************************/
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

/************************************************************************
* Function:	WlzObjFactsAffineTrans					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of an affine transform.	*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzAffineTransform *trans: Given affine transform.	*
************************************************************************/
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
      errNum = WlzObjFactsAppend(fData, "tx: %g.\n", trans->tx);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "ty: %g.\n", trans->ty);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "tz: %g.\n", trans->tz);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "scale: %g.\n", trans->scale);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "theta: %g.\n", trans->theta);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "phi: %g.\n", trans->phi);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "alpha: %g.\n", trans->alpha);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "psi: %g.\n", trans->psi);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "xsi: %g.\n", trans->xsi);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzObjFactsAppend(fData, "invert: %d.\n", trans->invert);
    }
  }
  --(fData->indent);
  return(errNum);
}

/************************************************************************
* Function:	WlzObjFactsPolygon2D					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of an 2D polygon domain.	*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzPolygonDomain *poly: Given polygon domain.		*
************************************************************************/
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
	    errNum = WlzObjFactsAppend(fData, "%8g %8d %8g.\n",
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
	    errNum = WlzObjFactsAppend(fData, "%8g %8d %8g.\n",
	    			       idx, tDVx0->vtX, tDVx0->vtY);
	    ++idx;
	    ++tDVx0;
	  }
	  break;
      }
      --(fData->indent);
    }
  }
  --(fData->indent);
  return(errNum);
}

/************************************************************************
* Function:	WlzObjFactsBoundlist					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of a boundary list.		*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzBoundList *bList:	Given boundary list.		*
************************************************************************/
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

/************************************************************************
* Function:	WlzObjFactsBoundlist					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Produces a text description of a histogram domain.	*
* Global refs:	-							*
* Parameters:	WlzObjFactsData *fData:	Facts data structure.		*
*		WlzObject *obj:		Object to determine domain	*
*					type.				*
*		WlzHistogramDomain *hist: Given histogram domain.	*
************************************************************************/
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
