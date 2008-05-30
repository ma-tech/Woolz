#if defined(__GNUC__)
#ident "MRC HGU $Id:"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id:"
#else static char _WlzDomainToTiles_c[] = "MRC HGU $Id:";
#endif
#endif
/*!
* \file         WlzDomainToTiles.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed May 28 15:55:01 2008
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
* \brief        cut a woolz domain object into tiles as per the standard tiff tile algorithm.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzdomaintotiles WlzDomainToTiles
\par Name
WlzDomainToTiles - 
\par Synopsis
\verbatim
WlzDomainToTiles

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

static void usage(
  char *proc_str)
{
  fprintf(stderr,
	  "Usage:\n"
	  "%s [-o <file>] [-t #,#] [-h] [-v]  [<input file>]\n"
	  "\tRead in a woolz 2D domain object and generate a set of tile\n"
	  "\timages using the same algorithm as for a tiled tiff. The image\n"
	  "\torigin is set to 0,0 and tiles generated row-wise from top to\n"
	  "\tbottom. The output is a copound object with the object order\n"
	  "\tin tile sequence order.\n"
	  "Arguments:\n"
	  "\t-o <file>         write object to given file, default stdout\n"
	  "\t-t xsize,ysize    tile size, default 256x256\n"
	  "\t-h                print this message\n"
	  "\t-v                verbose operation\n"
	  "\n",
	  proc_str);
  return;
}

int main(
  int   argc,
  char  **argv)
{
  FILE		*inFile, *outFile;
  char 		optList[] = "o:t:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  const char	*errMsg;
  int		verboseFlg=0;
  WlzObject	*inObj, *tileObj, *obj1, *obj2, *obj3;
  WlzCompoundArray	*cObj;
  int		tileWidth, tileHeight, width, height;
  int		i, tileIndx, row, col, numTiles;
  WlzIBox2	bbox;
  WlzDomain	domain;
  WlzValues	values;
  WlzPixelV	bckgrnd;

  /* set defaults */
  outFile = stdout;
  tileWidth = 256;
  tileHeight = 256;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'o':
      if( (outFile = fopen(optarg, "wb")) == NULL ){
	fprintf(stderr, "%s: couldn't open output file: %s\n", argv[0], optarg);
	usage(argv[0]);
	return 1;
      }
      break;

    case 't':
      i = sscanf(optarg, "%d,%d", &tileWidth, &tileHeight);
      if( i < 2 ){
	fprintf(stderr, "%s: couldn't read tile size\n", argv[0]);
	usage(argv[0]);
	return 1;
      }
      if((tileWidth < 1) || (tileHeight < 1)){
	fprintf(stderr, "%s: invalid tile size %dx%d\n", argv[0], tileWidth, tileHeight);
	usage(argv[0]);
	return 1;
      }
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* get input stream */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  /* read the first object */
  if((inObj = WlzReadObj(inFile, &errNum)) != NULL){
    switch( inObj->type ){

    case WLZ_2D_DOMAINOBJ:
      /* calculate number of tiles and set up compound object for output */
      bbox = WlzBoundingBox2I(inObj, &errNum);
      width = bbox.xMax - bbox.xMin + 1;
      height = bbox.yMax - bbox.yMin + 1;
      row = width / tileWidth + ((width%tileWidth)?1:0);
      col = height / tileHeight + ((height%tileHeight)?1:0);
      numTiles = row * col;
      cObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, numTiles, NULL,
				  WLZ_2D_DOMAINOBJ, &errNum);

      /* shift object bounding box to origin */
      obj1 = WlzShiftObject(inObj, -bbox.xMin, -bbox.yMin, 0, &errNum);
      WlzFreeObj(inObj);
      inObj = obj1;

      /* create tiles shifting each to origin */
      domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				       0, tileHeight-1, 0, tileWidth-1, &errNum);
      values.core = NULL;
      tileObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL, NULL, &errNum);
      bckgrnd = WlzGetBackground(inObj, &errNum);
      tileIndx = 0;
      for(row=0; row < height; row += tileHeight){
	for(col=0; col < width; col += tileWidth){
	  obj1 = WlzShiftObject(tileObj, col, row, 0, &errNum);
	  obj2 = WlzIntersect2(inObj, obj1, &errNum);
	  obj3 = WlzGreyTemplate(inObj, obj2, bckgrnd, &errNum);
	  cObj->o[tileIndx] = WlzAssignObject(
	    WlzShiftObject(obj3, -col, -row, 0, &errNum), NULL);
	  WlzFreeObj(obj3);
	  WlzFreeObj(obj2);
	  WlzFreeObj(obj1);
	  tileIndx++;
	}
      }
      WlzFreeObj(tileObj);

      break;

    default:
      fprintf(stderr,
	      "%s: invalid object type - only implmented for 2D domain objects.\n",
	      argv[0]);
      return 1;
    }
  }

  /* write object to stdout */
  if( errNum == WLZ_ERR_NONE ){
    WlzWriteObj(outFile, (WlzObject *) cObj);
  }
  else {
    (void) WlzStringFromErrorNum(errNum, &errMsg);
    fprintf(stderr, "%s: failed to generate the tile object: %s\n", argv[0], errMsg);
    return errNum;
  }

  return 0;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
