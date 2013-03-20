#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzClassMap_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzClassMap.c
* \author       Richard Baldock
* \date         October 2006
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
* \ingroup      BinWlzApp
* \brief        Generate a class image from a set of probability
*               images, setting colours from a standard list.
*
* \par Binary
* \ref wlzclassmap "WlzClassMap"
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzclassmap WlzClassMap
\par Name
WlzClassMap - 
\par Synopsis
\verbatim
WlzClassMap

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
#include <WlzExtFF.h>

/* insert code here */
/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static int cmpObjSizeA(
  const void *p1,
  const void *p2)
{
  WlzObject	**obj1 = (WlzObject **) p1;
  WlzObject     **obj2 = (WlzObject **) p2;
  int		size1, size2;

  size1 = WlzArea(*obj1, NULL);
  size2 = WlzArea(*obj2, NULL);

  if( size1 > size2 ){
    return 1;
  }
  else if( size1 < size2 ){
    return -1;
  }
  return 0;
}

static int cmpObjSizeB(
  const void *p1,
  const void *p2)
{
  WlzObject	**obj1 = (WlzObject **) p1;
  WlzObject     **obj2 = (WlzObject **) p2;
  int		size1, size2;

  size1 = WlzArea(*obj1, NULL);
  size2 = WlzArea(*obj2, NULL);

  if( size1 < size2 ){
    return 1;
  }
  else if( size1 > size2 ){
    return -1;
  }
  return 0;
}

static WlzUByte colorTable[][3] =
{
  {255,0,0},
  {0,255,0},
  {0,0,255},
  {255,191,0},
  {255,0,191},
  {191,255,0},
  {0,255,191},
  {191,0,255},
  {0,191,255},
  {255,127,64},
  {255,64,127},
  {127,255,64},
  {64,255,127},
  {127,64,255},
  {64,127,255},
  {255,255,0},
  {255,0,255},
  {0,255,255},
  {255,96,96},
  {96,255,96},
  {96,96,255}
};

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-b <background col>] [-B <zero col>] [-d <domainfile>]"
	  "[-i <indexfile>][-n <num>] [-N] [-s] [-S] [-h] [<input file>]\n"
	  "\tGenerate a classification map from a series of input image.\n"
	  "\tEach image is assumed to represent a class \"probability\", scaled\n"
	  "\tbetween 0 and 255. For each pixel the class with the maximum\n"
	  "\tprobability value is selected and the colour set appropriately.\n"
	  "\tNote colour values are 0-255\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -b r,g,b  Colour for background, i.e. outside of the domain:\n"
	  "\t            default (0,0,0)\n"
	  "\t  -B r,g,b  Colour for zero values  default (0,0,0)\n"
	  "\t  -d <file> Optional input domain over which class-map\n"
	  "\t            will be calculated, default - union of input.\n"
	  "\t  -i <file> Create a colour index strip image (jpeg)\n"
	  "\t  -n#       Maximum number of objects -default=100\n"
	  "\t  -N        Normalise the input objects.\n"
	  "\t  -o        Display an outline of the domain\n"
	  "\t  -s        Sort on size - smallest first.\n"
	  "\t  -S        Sort on size - largest first.\n"
	  "\t  -h        Help - prints this usage message\n",
	  proc_str,
	  WlzVersion());
  return;
}

int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj1, *obj2, *obj, **objlist, *rtnObj, *structElem;
  WlzObjectType	type = (WlzObjectType) -1;
  int 		i, n, nmax;
  FILE		*inFile;
  char 		optList[] = "b:B:d:i:n:NosSh";
  int		option;
  char 		*indexFile=NULL;
  /* int	normFlg=0; */
  int		outlineFlg=0;
  int		sortFlg=0;
  WlzPixelV	bckgrnd;
  WlzPixelV	zeroVal;
  int		red, green, blue;
  int		width, height;
  int		*buf;
/*  WlzUByte	**colorTable; */
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyValueWSpace	**gVWSpArr;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* setup the return object type and background */
  bckgrnd.type = WLZ_GREY_RGBA;
  bckgrnd.v.rgbv = 0;
  zeroVal.type = WLZ_GREY_RGBA;
  zeroVal.v.rgbv = 0;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  nmax = 100;
  rtnObj = NULL;
  obj = NULL;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'b':
      if( sscanf(optarg, "%d,%d,%d", &red, &green, &blue) < 3 ){
	fprintf(stderr, "%s: not enough colour values read in, using default\n",
		argv[0]);
      }
      else {
	WLZ_RGBA_RGBA_SET(bckgrnd.v.rgbv, red, green, blue, 255);
      }
      break;

    case 'B':
      if( sscanf(optarg, "%d,%d,%d", &red, &green, &blue) < 3 ){
	fprintf(stderr, "%s: not enough colour values read in, using default\n",
		argv[0]);
      }
      else {
	WLZ_RGBA_RGBA_SET(zeroVal.v.rgbv, red, green, blue, 255);
      }
      break;

      /* read in a target domain over which to capture the occupancy,
	 set grey type to WLZ_GREY_INT and value to zero */
    case 'd':
      if((inFile = fopen(optarg, "rb")) != NULL){
	if( (obj = WlzReadObj(inFile, &errNum)) == NULL ){
	  fprintf(stderr, "%s: Can't read class-map domain file\n", argv[0]);
	  usage(argv[0]);
	}
	fclose(inFile);
      }
      else {
	fprintf(stderr, "%s: Can't open class-map domain file\n", argv[0]);
	usage(argv[0]);
      }
      break;

    case 'i':
      indexFile = AlcStrDup(optarg);
      break;

    case 'n':
      nmax = atoi(optarg);
      if( nmax < 1 ){
	fprintf(stderr, "%s: nmax = %d is invalid\n", argv[0], nmax);
	usage(argv[0]);
	return( 1 );
      }
      break;

    case 'N':
      /* normFlg = 1; */
      break;

    case 'o':
      outlineFlg = 1;
      break;

    case 's':
      sortFlg = 1;
      break;

    case 'S':
      sortFlg = 2;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  /* check for input domain */
  if( obj ){
    width = obj->domain.i->lastkl - obj->domain.i->kol1 + 7;
    height = obj->domain.i->lastln - obj->domain.i->line1 + 7;
    buf = (int *) AlcMalloc(sizeof(WlzUInt) * width * height);
    rtnObj = WlzMakeRect(obj->domain.i->line1-3,
			 obj->domain.i->lastln+3,
			 obj->domain.i->kol1-3,
			 obj->domain.i->lastkl+3,
			 WLZ_GREY_RGBA,
			 buf,
			 bckgrnd,
			 NULL, NULL, &errNum);
    WlzGreySetValue(rtnObj, bckgrnd);
    if( outlineFlg ){
      structElem = WlzMakeCircleObject(2.0, 0.0, 0.0, &errNum);
      obj2 = WlzStructDilation(obj, structElem, &errNum);
      obj2->values = WlzAssignValues(rtnObj->values, &errNum);
      WLZ_RGBA_RGBA_SET(bckgrnd.v.rgbv, 255, 255, 255, 255);
      WlzGreySetValue(obj2, bckgrnd);
      WlzFreeObj(obj2);
      WlzFreeObj(structElem);
    }
    obj->values = WlzAssignValues(rtnObj->values, &errNum);
    WlzGreySetValue(obj, zeroVal);
    WlzFreeObj(obj);
    obj = NULL;
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  /* allocate space for the object pointers */
  if( (objlist = (WlzObject **)
       AlcMalloc(sizeof(WlzObject *) * nmax)) == NULL ){
    (void )fprintf(stderr, "%s: memory allocation failed.\n",
    		   argv[0]);
    return( 1 );
  }

  /* read objects accumulating compatible types */
  n = 0;
  while( ((obj = WlzAssignObject(WlzReadObj(inFile, NULL),
  			         NULL)) != NULL) && (n < nmax) ) {

    if( type == -1 &&
	(obj->type == WLZ_2D_DOMAINOBJ) ){
      type = obj->type;
    }

    if( obj->type == type ){
      objlist[n++] = WlzAssignObject(obj, NULL);
    } else {
      WlzFreeObj( obj );
    }
  }

  if( type == WLZ_EMPTY_OBJ ){
    return( 0 );
  }

  /* sort on size - area or volume */
  if( sortFlg == 1 ){
    /* smallest first */
    qsort(objlist, n, sizeof(WlzObject *), cmpObjSizeA);
  }
  else if ( sortFlg == 2 ){
    /* largest first */
    qsort(objlist, n, sizeof(WlzObject *), cmpObjSizeB);
  }

  /* check for class-map object */
  if( rtnObj ){
    /* check type against return object - must be the same */
    if( type != rtnObj->type ){
      fprintf(stderr, "%s: Class-map domain object type does not match\n"
	      " input domain type\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    /* use the union of input domains as the occupancy object */
    if((obj1 = WlzUnionN(n, objlist, 0, &errNum)) == NULL) {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to perform union (%s).\n",
		     argv[0], errMsg);
      return(1);
    }
    width = obj1->domain.i->lastkl - obj1->domain.i->kol1 + 7;
    height = obj1->domain.i->lastln - obj1->domain.i->line1 + 7;
    buf = (int *) AlcMalloc(sizeof(WlzUInt) * width * height);
    rtnObj = WlzMakeRect(obj1->domain.i->line1-3,
			 obj1->domain.i->lastln+3,
			 obj1->domain.i->kol1-3,
			 obj1->domain.i->lastkl+3,
			 WLZ_GREY_RGBA,
			 buf,
			 bckgrnd,
			 NULL, NULL, &errNum);
    WlzGreySetValue(rtnObj, bckgrnd);
    if( outlineFlg ){
      structElem = WlzMakeCircleObject(2.0, 0.0, 0.0, &errNum);
      obj2 = WlzStructDilation(obj1, structElem, &errNum);
      obj2->values = WlzAssignValues(rtnObj->values, &errNum);
      WLZ_RGBA_RGBA_SET(bckgrnd.v.rgbv, 255, 255, 255, 255);
      WlzGreySetValue(obj2, bckgrnd);
      WlzFreeObj(obj2);
      WlzFreeObj(structElem);
    }
    obj1->values = WlzAssignValues(rtnObj->values, &errNum);
    WlzGreySetValue(obj1, zeroVal);
    WlzFreeObj(obj1);
  }

  /* set up a colour table */
/*  AlcUnchar2Malloc(&colorTable, n, 3);*/
/* use static table until we can work out a decent algorithm */
  

  /* Find class for each pixel and set colour */
  errNum = WlzInitGreyScan(rtnObj, &iwsp, &gwsp);
  gVWSpArr = (WlzGreyValueWSpace **) AlcMalloc(sizeof(WlzGreyValueWSpace *) * n);
  for(i=0; i < n; i++){
    gVWSpArr[i] = WlzGreyValueMakeWSp(objlist[i], &errNum);
  }
  while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
    WlzGreyP		gptr;
    int	k, max, class;
    gptr = gwsp.u_grintptr;
    for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++, gptr.rgbp++){
      max = 0;
      class = -1;
      for(i=0; i < n; i++){
	if( WlzInsideDomain(objlist[i],0,iwsp.linpos, k, &errNum) ){
	  WlzGreyValueGet(gVWSpArr[i], 0, iwsp.linpos, k);
	  if( gVWSpArr[i]->gVal[0].inv > max ){
	    max = gVWSpArr[i]->gVal[0].inv;
	    class = i;
	  }
	}
      }
      if( class > -1 ){
	WLZ_RGBA_RGBA_SET(*gptr.rgbp,
			  colorTable[class][0],
			  colorTable[class][1],
			  colorTable[class][2],
			  255);
      }
    }
  }
  if( errNum == WLZ_ERR_EOO ){
    errNum = WLZ_ERR_NONE;
  }

  

  /* output class-map object object */
  WlzWriteObj(stdout, rtnObj);

  /* check for index object */
  if( indexFile ){
    FILE	*fp;
    WlzUInt	uiVal;
    int		k;
    if((fp = fopen(indexFile, "wb")) != NULL){
      WlzFreeObj(rtnObj);
      bckgrnd.type = WLZ_GREY_RGBA;
      bckgrnd.v.rgbv = 0;
      width = 25;
      height = 25 * n;
      buf = (int *) AlcMalloc(sizeof(WlzUInt) * width * height);
      rtnObj = WlzMakeRect(0, height-1, 0, width-1,
			   WLZ_GREY_RGBA, buf, bckgrnd,
			   NULL, NULL, &errNum);
      for(i=0; i < n; i++){
	WLZ_RGBA_RGBA_SET(uiVal,
			  colorTable[i][0],
			  colorTable[i][1],
			  colorTable[i][2],
			  255);
	for(k=0; k < (width*width); k++, buf++){
	  *((WlzUInt *) buf) = uiVal;
	}
      }
      WlzEffWriteObj(fp, NULL, rtnObj, WLZEFF_FORMAT_JPEG);
      fclose(fp);
    }
  }

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
