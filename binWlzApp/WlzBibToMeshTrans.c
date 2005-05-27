#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzBibToMeshTrans.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Mar  9 13:33:37 2005
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzTransform
* \brief        Generate a mesh transform from an MAPaint bibfiile of
 Basis function parameters and tie-points. Also required is a 2D image
 in order to generate a suitable mesh.
*               
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/
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
    if( (inFP = fopen(*(argv+optind), "r")) == NULL ){
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
    if( bibFP = fopen(bibFile, "r") ){
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
      if( meshTr = WlzMeshTransformFromCPts(inObj, wlzFnType, basisFnPolyOrder,
					    numVtxs, srcVtxs, numVtxs, dstVtxs,
					    meshMthd, meshMinDst, meshMaxDst,
					    &errNum) ){
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
 
