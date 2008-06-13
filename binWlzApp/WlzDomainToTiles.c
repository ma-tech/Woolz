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
* \brief        Cut a woolz domain object into tiles or chips of equal size. With the tiff-tile option the object is reset to origin (0,0) and each tile has origin (0,0) and the tiles are ordered in the tiff-tile sequence.
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
	  "%s [-o <file>] [-p #] [-t #,#,#] [-T] [-h] [-v]  [<input file>]\n"
	  "\tRead in a woolz 2D domain object and generate a set of tile\n"
	  "\timages. These will be a set of equal sized 2D or 3D tiles which\n"
	  "\tcover the input domain. If the input object has values then these\n"
	  "\twill also be transferred.\n"
	  "\t\n"
	  "\tWith option \"-T\" use the same algorithm as for a tiled tiff. The image\n"
	  "\torigin is set to (0,0) and tiles generated row-wise from top to\n"
	  "\tbottom. The output is a compound object with the object order\n"
	  "\tin tile sequence order.\n"
	  "Arguments:\n"
	  "\t-o <file>         write object to given file, default stdout\n"
	  "\t-p percent     only generate tiles if intersect size is >= percent \n"
	  "\t               of the maximum tile size.\n"
	  "\t-t xsize,ysize,zsize    tile size, default 256x256 (2D), 16x16x16 (3D)\n"
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
  char 		optList[] = "o:p:t:Thv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  const char	*errMsg;
  int		verboseFlg;
  WlzObject	*inObj, *tileObj, *obj1, *obj2, *obj3;
  WlzCompoundArray	*cObj;
  int		tileWidth, tileHeight, tileDepth, width, height, depth;
  int		i, tileIndx, row, col, plane, numTiles;
  int		tileDefaultFlg;
  int		tiffFlg;
  int		percent, minSize;
  WlzIBox2	bbox2;
  WlzIBox3	bbox3;
  WlzDomain	domain;
  WlzValues	values;
  WlzPixelV	bckgrnd;

  /* set defaults */
  outFile = stdout;
  tileDefaultFlg = 1;
  verboseFlg = 0;
  tiffFlg = 0;
  tileWidth = tileHeight = tileDepth = 0;
  percent = 0;

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

    case 'p':
      percent = atoi(optarg);
      if((percent < 0) || (percent > 100)){
	fprintf(stderr, "%s: invalid size percentage: %d\n", argv[0], percent);
	usage(argv[0]);
	return 1;
      }
      break;

    case 't':
      i = sscanf(optarg, "%d,%d,%d", &tileWidth, &tileHeight, &tileDepth);
      if( i < 3 ){
	fprintf(stderr, "%s: couldn't read tile size\n", argv[0]);
	usage(argv[0]);
	return 1;
      }
      if((tileWidth < 1) || (tileHeight < 1) || (tileDepth < 1)){
	fprintf(stderr, "%s: invalid tile size %dx%dx%d\n", argv[0],
		tileWidth, tileHeight, tileDepth);
	usage(argv[0]);
	return 1;
      }
      tileDefaultFlg = 0;
      break;

    case 'T':
      tiffFlg = 1;
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
      /* check tile default */
      if( tileDefaultFlg ){
	tileWidth = tileHeight = 256;
      }
      minSize = tileWidth * tileHeight * percent / 100;

      /* calculate number of tiles and set up compound object for output */
      bbox2 = WlzBoundingBox2I(inObj, &errNum);
      width = bbox2.xMax - bbox2.xMin + 1;
      height = bbox2.yMax - bbox2.yMin + 1;
      row = width / tileWidth + ((width%tileWidth)?1:0);
      col = height / tileHeight + ((height%tileHeight)?1:0);
      numTiles = row * col;
      cObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, numTiles, NULL,
				  WLZ_2D_DOMAINOBJ, &errNum);

      /* create tiles shifting each to origin */
      domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				       0, tileHeight-1, 0, tileWidth-1, &errNum);
      values.core = NULL;
      tileObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL, NULL, &errNum);

      if( inObj->values.core ){
	bckgrnd = WlzGetBackground(inObj, &errNum);
      }
      tileIndx = 0;
      for(row=bbox2.yMin; row <= bbox2.yMax; row += tileHeight){
	for(col=bbox2.xMin; col <= bbox2.xMax; col += tileWidth){
	  obj1 = WlzShiftObject(tileObj, col, row, 0, &errNum);
	  if((obj2 = WlzIntersect2(inObj, obj1, &errNum))){
	    if( WlzArea(obj2, &errNum) > minSize ){
	      if( inObj->values.core ){
		obj3 = WlzGreyTemplate(inObj, obj2, bckgrnd, &errNum);
	      }
	      else {
		obj3 = WlzCopyObject(obj2, &errNum);
	      }
	      if( tiffFlg ){
		cObj->o[tileIndx] = WlzAssignObject(
		  WlzShiftObject(obj3, -col, -row, 0, &errNum), NULL);
		WlzFreeObj(obj3);
	      }
	      else {
		cObj->o[tileIndx] = WlzAssignObject(obj3, &errNum);
	      }
	      tileIndx++;
	    }
	    WlzFreeObj(obj2);
	  }
	  WlzFreeObj(obj1);
	}
      }
      cObj->n = tileIndx;
      WlzFreeObj(tileObj);
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check tile default */
      if( tileDefaultFlg ){
	tileWidth = tileHeight = tileDepth = 16;
      }
      minSize = tileWidth * tileHeight *tileDepth * percent / 100;

      /* calculate number of tiles and set up compound object for output */
      bbox3 = WlzBoundingBox3I(inObj, &errNum);
      width = bbox3.xMax - bbox3.xMin + 1;
      height = bbox3.yMax - bbox3.yMin + 1;
      depth = bbox3.zMax - bbox3.zMin + 1;
      row = width / tileWidth + ((width%tileWidth)?1:0);
      col = height / tileHeight + ((height%tileHeight)?1:0);
      plane = depth / tileDepth + ((depth%tileDepth)?1:0);
      numTiles = plane * row * col;
      cObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, numTiles, NULL,
				  WLZ_3D_DOMAINOBJ, &errNum);

      /* create tiles shifting each to origin if needed */
      domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				    0, tileDepth - 1,
				    0, tileHeight-1,
				    0, tileWidth-1, &errNum);
      values.core = NULL;
      tileObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values, NULL, NULL, &errNum);
      for(plane=0; plane < tileDepth; plane++){
	domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					 0, tileHeight-1, 0, tileWidth-1, &errNum);
	tileObj->domain.p->domains[plane] = WlzAssignDomain(domain, &errNum);
      }

      /* get background if input object has grey-values */
      if( inObj->values.core ){
	bckgrnd = WlzGetBackground(inObj, &errNum);
      }

      tileIndx = 0;
      for(plane=bbox3.zMin; plane <= bbox3.zMax; plane += tileDepth){
	for(row=bbox3.yMin; row <= bbox3.yMax; row += tileHeight){
	  for(col=bbox3.xMin; col <= bbox3.xMax; col += tileWidth){
	    obj1 = WlzShiftObject(tileObj, col, row, plane, &errNum);
	    if((obj2 = WlzIntersect2(inObj, obj1, &errNum))){
	      if( WlzVolume(obj2, &errNum) > minSize ){
		if( inObj->values.core ){
		  obj3 = WlzGreyTemplate(inObj, obj2, bckgrnd, &errNum);
		}
		else {
		  obj3 = WlzCopyObject(obj2, &errNum);
		}
		if( tiffFlg ){
		  cObj->o[tileIndx] = WlzAssignObject(
		    WlzShiftObject(obj3, -col, -row, -plane, &errNum), NULL);
		  WlzFreeObj(obj3);
		}
		else {
		  cObj->o[tileIndx] = WlzAssignObject(obj3, &errNum);
		}
		tileIndx++;
	      }
	      WlzFreeObj(obj2);
	    }
	    WlzFreeObj(obj1);
	  }
	}
      }
      WlzFreeObj(tileObj);
      cObj->n = tileIndx;
      break;

    default:
      fprintf(stderr,
	      "%s: invalid object type - only implmented for 2D & 3D domain objects.\n",
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
