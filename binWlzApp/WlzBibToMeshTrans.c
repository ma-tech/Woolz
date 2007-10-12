#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzBibToMeshTrans_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzApp/WlzBibToMeshTrans.c
* \author       Richard Baldock
* \date         March 2005
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
* \brief        Generates a mesh transform from an MAPaint bibfiile of
*		basis function parameters and tie-points.
*               
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzbibtomeshtrans "WlzBibToMeshTrans"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzbibtomeshtrans WlzBibToMeshTrans
\par Name
WlzBibToMeshTrans - generates a mesh transform from a bibfiile
		    with basis function parameters and tie-points.
\par Synopsis
\verbatim
WlzBibToMeshTrans [-h] [-v] [-b<bibfile>] [-f<obj>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Version operation.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>A bibfile defining the warp parameters e.g. from MAPaint.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>A 2D Woolz object which is used to define the mesh,
        this should be the same as the object to which the mesh
	transform will be applied.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
</table>
\par Description
Creates a mesh transform from an MAPaint warp bibfile.
By default the 2D object used to define the mesh transfrom
is read from the standard input and
the mesh transform is written to the standard output.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzFacts.c "WlzFacts.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref WlzMeshTransformFromCPts "WlzMeshTransformFromCPts(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>
#include <bibFile.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s -b <parameter bibfile>"
	  " [-o <output file>] [-h] [-v]"
	  " <2D object input file>\n"
	  "\tConvert an MAPaint warp bibfile to a mesh transform\n"
	  "\twriting the mesh transform object to standard output.\n"
	  "\tA 2D image is required to define the mesh extent, read\n"
	  "\tfrom standard input by default\n"
	  "\tOptions are:\n"
	  "\t  -b<bibfile>        bibfile defining the warp parameters e.g.\n"
	  "\t                     from MAPaint\n"
	  "\t  -f<2D image file>  2D Woolz image used to define the mesh\n"
	  "\t  -o<output file>    Output filename, default to stdout\n"
	  "\t  -h                 Help - this message\n"
	  "\t  -v                 verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*inObj, *outObj;
  WlzDomain	domain;
  WlzValues	values;
  FILE		*inFP, *outFP, *bibFP;
  char		*outFile, *bibFile;
  char 		optList[] = "b:f:o:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  int		verboseFlg=0;
  WlzMeshTransform *meshTr = NULL;
  WlzFnType		wlzFnType;
  WlzMeshGenMethod	meshMthd;
  int			basisFnPolyOrder = 3;
  int			meshMinDst;
  int	 		meshMaxDst;
  WlzTransformType	affineType;
  WlzDVertex2		*srcVtxs, *dstVtxs;
  int		relFlg=0;
  int		numVtxs, maxNumVtxs;
  char		*errMsg;
    
  /* additional defaults */
  outFile = "-";
  bibFile = NULL;
  numVtxs = 0;
  maxNumVtxs = 0;
  srcVtxs = NULL;
  dstVtxs = NULL;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'b':
      bibFile = optarg;
      break;

    case 'o':
      outFile = optarg;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 0;

    case 'v':
      verboseFlg = 1;
      break;

    }
  }

  /* check input file/stream */
  inFP = stdin;
  if( optind < argc ){
    if( (inFP = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* check output file/stream */
  if(strcmp(outFile, "-"))
  {
    if((outFP = fopen(outFile, "w")) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    outFP = stdout;
  }

  /* check bibfile - get warp function parameters and tie-points */
  if((errNum == WLZ_ERR_NONE) && (bibFile != NULL)){
    if((bibFP = fopen(bibFile, "r")) != NULL){
      bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, bibFP);
      while( bibFileErr == BIBFILE_ER_NONE ){

	/* test for warp function and mesh parameters */
	if( !strncmp(bibfileRecord->name, "WlzWarpTransformParams", 22) ){

	  WlzEffBibParseWarpTransformParamsRecord(bibfileRecord,
						  &wlzFnType,
						  &affineType,
						  &meshMthd,
						  &meshMinDst, &meshMaxDst);
	}

	/* test for tie points */
	if( !strncmp(bibfileRecord->name, "WlzTiePointVtxs", 15) ){
	  int		index;
	  WlzDVertex3	dstVtx;
	  WlzDVertex3	srcVtx;

	  WlzEffBibParseTiePointVtxsRecord(bibfileRecord, &index,
				       &dstVtx, &srcVtx, &relFlg);
	  if( relFlg ){
	    fprintf(stderr, "%s: warning, relative vertex positions, probably an error\n",
		    argv[0]);
	  }
	  if( numVtxs >= maxNumVtxs ){
	    maxNumVtxs += 512;
	    srcVtxs = (WlzDVertex2 *) AlcRealloc(srcVtxs,
						 sizeof(WlzDVertex2)*maxNumVtxs);
	    dstVtxs = (WlzDVertex2 *) AlcRealloc(dstVtxs,
						 sizeof(WlzDVertex2)*maxNumVtxs);
	  }
	  srcVtxs[numVtxs].vtX = srcVtx.vtX;
	  srcVtxs[numVtxs].vtY = srcVtx.vtY;
	  dstVtxs[numVtxs].vtX = dstVtx.vtX;
	  dstVtxs[numVtxs].vtY = dstVtx.vtY;
	  numVtxs++;
	}

	BibFileRecordFree(&bibfileRecord);
	bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, bibFP);
      }
      fclose( bibFP );
      if( bibFileErr == BIBFILE_ER_EOF ){
	bibFileErr = BIBFILE_ER_NONE;
      }
      if( bibFileErr != BIBFILE_ER_NONE ){
	fprintf(stderr, "%s: error reading bibfile\n", argv[0]);
  	return 1;
      }
    }
    else {
      fprintf(stderr, "%s: can't open parameter bibfile %s\n", argv[0],
	      bibFile);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: warp parameters and tie-points bibfile required\n",
	    argv[0]);
    usage(argv[0]);
    return 1;
  }
  
  /* read objects and section if possible */
  while((errNum == WLZ_ERR_NONE) &&
        ((inObj = WlzReadObj(inFP, &errNum)) != NULL))
  {
    switch( inObj->type )
    {
    case WLZ_2D_DOMAINOBJ:
      /* now create the mesh transform */
      if((meshTr = WlzMeshTransformFromCPts(inObj, wlzFnType, basisFnPolyOrder,
					    numVtxs, srcVtxs, numVtxs, dstVtxs,
					    meshMthd, meshMinDst, meshMaxDst,
					    &errNum)) != NULL){
	/* write out transform and free */
	domain.mt = meshTr;
	values.core = NULL;
	outObj = WlzMakeMain(WLZ_MESH_TRANS, domain, values, NULL, NULL, &errNum);
	WlzWriteObj(outFP, outObj);
	WlzFreeObj(outObj);
      }
      else {
	/* something wrong */
	fprintf(stderr,
		"%s: failed to generate a mesh transform, error code: %d\n",
		argv[0], errNum);
      }
      break;

    default:
      fprintf(stderr, "%s: 2D input object required\n", argv[0]);
      usage(argv[0]);
      return 1;
    }

    WlzFreeObj(inObj);
  }
  if(errNum == WLZ_ERR_READ_EOF)
  {
    errNum = WLZ_ERR_NONE;
  }

  return errNum;
}
 
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
