#if defined(__GNUC__)
#ident "MRC HGU $Id:"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id:"
#else
static char _WlzApplyTileFunction_c[] = "MRC HGU $Id:";
#endif
#endif
/*!
* \file         WlzApplyTileFunction.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Thu Oct 21 17:06:05 2010
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
* \brief        Apply a function for each tile to the associated image and output a vector of values.
*               
*
* Maintenance log with most recent changes at top of list.
*/
#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(
  char	*str)
{
  fprintf(stderr,
	  "Usage:\n"
	  "%s -f <function> -h -v <image file> <tiles file>\n"
	  "\tRead in an image file (woolz image object) and tiles file (woolz compound object)\n"
	  "\tand apply the selected function to the iage values within each tile. The calculated\n"
	  "\tvalues are output as an ascii csv string to standard output\n"
	  "Arguments:\n"
	  "\t-f#          parameter to determine tile function default 1\n"
	  "\t             = 1 - average\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n"
	  "\n",
	  str);

  return;
}


int main(
  int   argc,
  char  **argv)
{
  FILE		*inFileObj=NULL, *inFileTiles=NULL;
  char 		optList[] = "f:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		verboseFlg=0;
  int		func=1;
  WlzObject	*imageObj = NULL, *tilesObj=NULL, *obj;
  WlzCompoundArray	*tilesCompObj;
  double	tileVal=0.0;
  int		i;
  WlzPixelV	tmplVal;
  WlzGreyType	dstGType;
  double	dstMin, dstMax, dstSum, dstSumSq, dstMean, dstStdDev;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'f':
      func = atoi(optarg);
      switch (func){
      default:
	fprintf(stderr, "%s: invalid tile function option %d, reset to 1\n", argv[0], func);
	func = 1;
	break;
      case 1: /* average */
	break;
      }
	  
    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* must have image and tile file */
  /* check input file/stream */
  if( (argc - optind) < 2 ){
    fprintf(stderr, "%s: not enough arguments\n", argv[0]);
    usage(argv[0]);
    return 1;
  }
  else {
    if( (inFileObj = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open image object file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (inFileTiles = fopen(*(argv+optind+1), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open tiles file %s\n", argv[0], *(argv+optind+1));
      usage(argv[0]);
      return 1;
    }
  }
  
  /* verbose output */
  if( verboseFlg ){
    fprintf(stderr, "%s: tile function value value: %d\n", argv[0], func);
    fprintf(stderr, "%s: image object file: %s\n", argv[0], *(argv+optind));
    fprintf(stderr, "%s: tiles object file: %s\n", argv[0], *(argv+optind+1));
  }

  /* get woolz objects */
  if((imageObj = WlzReadObj(inFileObj, &errNum)) == NULL){
    fprintf(stderr, "%s: can't read image object\n", argv[0]);
    usage(argv[0]);
    return 1;
  }
  imageObj = WlzAssignObject(imageObj, NULL);

  if((tilesObj = WlzReadObj(inFileTiles, &errNum)) == NULL){
    fprintf(stderr, "%s: can't read tiles object\n", argv[0]);
    usage(argv[0]);
    return 1;
  }
  tilesObj = WlzAssignObject(tilesObj, NULL);
  tilesCompObj = (WlzCompoundArray *) tilesObj;

  /* now calculate the tile values */
  switch( tilesObj->type ){
  case WLZ_COMPOUND_ARR_1:
  case WLZ_COMPOUND_ARR_2:
    tmplVal.type = WLZ_GREY_UBYTE;
    tmplVal.v.ubv = 0;
    for(i=0; i < tilesCompObj->n; i++){
      switch( func ){
      case 1: /* average */
	/* get tile object with values */
	if( (obj = WlzGreyTemplate(imageObj, tilesCompObj->o[i], tmplVal, &errNum)) != NULL ){
	  WlzGreyStats(obj, &dstGType, &dstMin, &dstMax, &dstSum, &dstSumSq,
		       &dstMean, &dstStdDev, &errNum);
	  tileVal = dstMean;
	}
	else {
	  fprintf(stderr, "%s: error in tile template\n", argv[0]);
	  tileVal = -1.0;
	}
	break;

      default:
	break;
      }

      if( i > 0 ){
	fprintf(stdout, ", %8.3f", tileVal);
      }
      else{
	fprintf(stdout, "%8.3f", tileVal);
      }
    }
    break;

  default:
    fprintf(stderr, "%s: invalid tiles object type - must be compound array\n", argv[0]);
    return 1;
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
