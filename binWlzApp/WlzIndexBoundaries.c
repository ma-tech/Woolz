#if defined(__GNUC__)
#ident "MRC HGU $Id:"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id:"
#else static char _WlzIndexBoundaries_c[] = "MRC HGU $Id:";
#endif
#endif
/*!
* \file         WlzIndexBoundaries.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Thu May  1 22:11:59 2008
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par Copyright:
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
* \ingroup      BinWlzApp
* \brief        Calculate the boundaries for a set of threshold levels.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzindexboundaries WlzIndexBoundaries
\par Name
WlzIndexBoundaries - 
\par Synopsis
\verbatim
WlzIndexBoundaries

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-</b></td>
    <td> </td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>

\par Description

\par Examples
\verbatim
\endverbatim

\par See Also
\par Bugs
None known
*/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>

#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s  [-R] [-h] [-v] [<input file>]\n"
	  "\tOptions are:\n"
	  "\t  -R        Output polyline coordinates relative to\n"
	  "\t            the image bounding box.\n"
	  "\t  -h        Help - prints this usage message.\n"
	  "\t  -v        Verbose operation.\n"
	  "",
	  proc_str);
  return;
}

typedef struct _polyListItem {
  WlzPolygonDomain	*poly;
  int			index;
  int			area;
} PolyListItem;

static void freePolyListItem(
  void *item)
{
  PolyListItem	*polyItem = (PolyListItem *) item;

  if( polyItem != NULL ){
    if( polyItem->poly != NULL ){
      WlzFreePolyDmn( polyItem->poly );
      polyItem->poly = NULL;
    }
    AlcFree( polyItem );
  }

  return;
}

static int polyItemCompArea(
  void	*item1,
  void 	*item2)
{
  PolyListItem	*polyItem1 = (PolyListItem *) item1;
  PolyListItem	*polyItem2 = (PolyListItem *) item2;

  if( polyItem1->area > polyItem2->area ){
    return -1;
  }
  if( polyItem1->area < polyItem2->area ){
    return 1;
  }
  return 0;
}

WlzErrorNum addPolylineToList(
  AlcDLPList	*polyList,
  int		index,
  WlzBoundList	*boundList)
{
  PolyListItem	*polyItem;
  WlzObject	*tmpObj;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* check the boundary type */
  switch( boundList->type ){
  case WLZ_BOUNDLIST_PIECE:
    if( boundList->poly ){
      polyItem = AlcMalloc(sizeof(PolyListItem));
      polyItem->poly = WlzAssignPolygonDomain(boundList->poly, &errNum);
      polyItem->index = index;
      if((tmpObj = WlzPolyToObj(boundList->poly, WLZ_SIMPLE_FILL, &errNum))){
	polyItem->area = WlzArea(tmpObj, &errNum);
	WlzFreeObj(tmpObj);
      }
      AlcDLPListEntryAppend(polyList, NULL, (void *) polyItem, freePolyListItem);
    }
    break;

  case WLZ_BOUNDLIST_HOLE:
  default:
    break;
  }

  /* check next and down */
  if( boundList->next ){
    addPolylineToList(polyList, index, boundList->next);
  }
  if( boundList->down ){
    addPolylineToList(polyList, index, boundList->down);
  }

  return errNum;
}

int main(
  int   argc,
  char  **argv)
{
  WlzObject     *obj, *obj1, *obj2, *obj3;
  WlzObject	*boundObj;
  WlzPixelV	minG, maxG;
  WlzIVertex2	*vtxs;
  int		i, index;
  int		relativeFlg;
  int		verboseFlg;
  FILE		*inFile;
  char 		optList[] = "Rhv";
  int		option;
  AlcDLPList	*polyList;
  AlcDLPItem	*listItem;
  PolyListItem	*polyItem;
  WlzIBox3	bBox3;
  AlcErrno	alcErr;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* default values */
  relativeFlg = 0;
  verboseFlg = 0;
  
  /* check input arguments */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'R':
      relativeFlg = 1;
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  if((obj = WlzReadObj(inFile, &errNum)) != NULL) 
  {
    obj = WlzAssignObject(obj, NULL);
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
	/* check for value table */
	if( obj->values.core == NULL ){
	  fprintf(stderr, "%s: object with value table required\n", argv[0]);
	  usage(argv[0]);
	  return( 1 );
	}

	/* check grey-type - require index image with integer values */
	switch( WlzGreyTypeFromObj(obj, &errNum) ){
	case WLZ_GREY_UBYTE:
	case WLZ_GREY_SHORT:
	case WLZ_GREY_INT:
	  break;

	default:
	  fprintf(stderr, "%s: object with integer grey values required\n", argv[0]);
	  usage(argv[0]);
	  return( 1 );
	}

	/* threshold at each level and build list of polygons */
	polyList = AlcDLPListNew(&alcErr);
	WlzGreyRange(obj, &minG, &maxG);
	WlzValueConvertPixel(&minG, minG, WLZ_GREY_INT);
	WlzValueConvertPixel(&maxG, maxG, WLZ_GREY_INT);
	index = (minG.v.inv < 1) ? 1 : minG.v.inv;
	minG.v.inv = index;
	for(;index < maxG.v.inv; index++){
	  minG.v.inv = index;
	  if( obj1 = WlzThreshold(obj, minG, WLZ_THRESH_HIGH, &errNum) ){
	    if( WlzIsEmpty(obj1, &errNum) ){
	      break;
	    }
	    minG.v.inv++;
	    obj2 = WlzThreshold(obj, minG, WLZ_THRESH_HIGH, &errNum);
	    if( (obj2 == NULL) || (WlzIsEmpty(obj2, &errNum)) ){
	      boundObj = WlzObjToBoundary(obj1, 1, &errNum);
	      addPolylineToList(polyList, index, boundObj->domain.b);
	      WlzFreeObj(boundObj);
	    }
	    else {
	      obj3 = WlzDiffDomain(obj1, obj2, &errNum);
	      if( obj3 ){
		if( !WlzIsEmpty(obj3, &errNum) ){
		  boundObj = WlzObjToBoundary(obj3, 1, &errNum);
		  addPolylineToList(polyList, index, boundObj->domain.b);
		  WlzFreeObj(boundObj);
		}
		WlzFreeObj(obj3);
	      }
	    }
	    if( obj2 ){
	      WlzFreeObj(obj2);
	    }
	    WlzFreeObj(obj1);
	  }
	}

	/* get the last polygons */
	minG.v.inv = index;
	if( obj1 = WlzThreshold(obj, minG, WLZ_THRESH_HIGH, &errNum) ){
	  if( !WlzIsEmpty(obj1, &errNum) ){
	    boundObj = WlzObjToBoundary(obj1, 1, &errNum);
	    addPolylineToList(polyList, index, boundObj->domain.b);
	    WlzFreeObj(boundObj);
	  }
	  WlzFreeObj(obj1);
	}

	/* sort polygon list with respect to area */
	AlcDLPListSort(polyList, polyItemCompArea);

	/* write polygon list to stdout */
	if( AlcDLPListCount(polyList, &alcErr) > 0 ){
	  if( relativeFlg ){
	    bBox3 = WlzBoundingBox3I(obj, &errNum);
	  }
	  fprintf(stdout, " \"polylines\":[\n");
	  listItem = polyList->head;
	  do {
	    polyItem = (PolyListItem *) listItem->entry;
	    fprintf(stdout, "  {\"index\": \"%d\", \"coords\": \"", polyItem->index);

	    vtxs = polyItem->poly->vtx;
	    for(i=0; i < polyItem->poly->nvertices; i++, vtxs++){
	      if( relativeFlg ){
		fprintf(stdout, "%d,%d",
			vtxs->vtX - bBox3.xMin,
			vtxs->vtY - bBox3.yMin);
	      }
	      else {
		fprintf(stdout, "%d,%d", vtxs->vtX, vtxs->vtY);
	      }
	      if( i < (polyItem->poly->nvertices - 1) ){
		fprintf(stdout, ",");
	      }
	      else {
		fprintf(stdout, "\"");
	      }
	    }
	    
	    if( listItem->next == polyList->head ){
	      fprintf(stdout, "}\n");
	    }
	    else {
	      fprintf(stdout, "},\n");
	    }

	    listItem = listItem->next;
	  } while( listItem != polyList->head );

	  fprintf(stdout, "  ]\n");
	}
	break;

	/* Only 2D for now */
      case WLZ_3D_DOMAINOBJ:
      default:
        WlzWriteObj(stdout,obj);
        break;
    }
 
    WlzFreeObj(obj);
  }

  return 0;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
