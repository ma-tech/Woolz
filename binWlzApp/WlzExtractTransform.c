#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtractTransform_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzApp/WlzExtractTransform.c
* \author       Richard Baldock
* \date         April 1999
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
* \brief	Finds the affine transforms from a reconstruction
* 		bibfile.
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzextracttransform "WlzExtractTransform"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzextracttransform WlzExtractTransform
\par Name
WlzExtractTransform - finds the affine transforms from a reconstruction
                      bibfile.
\par Synopsis
\verbatim
WlzExtractTransform [-h] [-v] [-A] [-R] [-a] [-t] [-T] [-f<filename>]
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
    <td><b>-A</b></td>
    <td>Calculate absolute transform.</td>
  </tr>
  <tr> 
    <td><b>-R</b></td>
    <td>Calculate relative transform, default.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Output and argument string for WlzAffineTransformObj, default.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Output a Woolz affine transform object.</td>
  </tr>
  <tr> 
    <td><b>-T</b></td>
    <td>Output text.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Filename to be found.</td>
  </tr>
</table>
\par Description
WlzExtractTransform finds the transform for the given filename
within the bibfile and outputs either
an argument string for WlzAffineTransformObj,
a Woolz affine transform or
a text description.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzExtractTransform.c "WlzExtractTransform.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include <Wlz.h>
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
	  "Usage:\t%s [-h] [-v] [-A] [-R] [-a] [-t] [-T] [-f<filename>]\n"
	  "                          [input bibfile]\n"
	  "\tFinds the transform for the given filename\n"
	  "\twithin the bib-file and output an argument\n"
	  "\tstring for WlzAffineTransformObj\n"
	  "\tOptions are:\n"
	  "\t  -h                Help - prints this usage message\n"
	  "\t  -v                verbose operation\n"
	  "\t  -A                calculate absolute transform\n"
	  "\t  -R                calculate relative transform (default)\n"
	  "\t  -a                output and argument string for\n"
	  "\t                    WlzAffineTransformObj (default)\n"
	  "\t  -t                output a woolz affine transform object\n"
	  "\t  -T                output text\n"
	  "\t  -f<filename>      filename to be found\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{
  FILE		*inFile;
  char 		optList[] = "ARatTf:h";
  int		option;
  double	tx=0.0, ty=0.0, theta=0.0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*filename=NULL, *recordFilename=NULL;
  char		*errMsg;
  int		numParsedFields=0;
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr=BIBFILE_ER_NONE;
  int		outputOpt=2;
  int		absoluteFlg=0;
  WlzAffineTransformPrim prim;
  WlzAffineTransform	*inTrans, *outTrans, *tmpTrans;
  WlzAffineTransform	*reconTrans;
  WlzObject	*outObj;
  WlzDomain	domain;
  WlzValues	values;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'A':
      absoluteFlg = 1;
      break;

    case 'R':
      absoluteFlg = 0;
      break;

    case 'a':
      outputOpt = 2;
      break;

    case 't':
      outputOpt = 1;
      break;

    case 'T':
      outputOpt = 0;
      break;

    case 'f':
      if( (optarg != NULL) && (strlen(optarg) > 0) ){
	filename = optarg;
      }
      else {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* open the bib file */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* create the first transform */
  inTrans = NULL;
  outTrans = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					   0.0, 0.0, 0.0, 1.0,
					   0.0, 0.0, 0.0, 0.0, 0.0, 0,
					   NULL);
  reconTrans = NULL;

  /* read bibfile records until the filename matches */
  bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
  while( bibFileErr == BIBFILE_ER_NONE ) 
  {
    if( !strncmp(bibfileRecord->name, "Section", 7) ){

      /* parse the record */
      numParsedFields = BibFileFieldParseFmt
	(bibfileRecord->field,
	 (void *) &tx, "%lg", "TransformTx",
	 (void *) &ty, "%lg", "TransformTy",
	 (void *) &theta, "%lg", "TransformTheta",
	 (void *) &recordFilename, "%s", "File",
	 NULL);

      /* make the transform for this record */
      inTrans = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					      tx, ty, 0.0, 1.0,
					      theta, 0.0, 0.0, 0.0, 0.0,
					      0, NULL);

      /* if absolute concatenate the transforms */
      if( absoluteFlg ){
	tmpTrans = WlzAffineTransformProduct(outTrans, inTrans, NULL);
	WlzFreeAffineTransform(inTrans);
	WlzFreeAffineTransform(outTrans);
	outTrans = tmpTrans;
      }
      else {
	WlzFreeAffineTransform(outTrans);
	outTrans = inTrans;
      }
	

      /* test for required file */
      if( filename && recordFilename && strstr(recordFilename, filename) ){
	BibFileRecordFree(&bibfileRecord);
	break;
      }
    }
    else if( !strncmp(bibfileRecord->name, "Reconstruction", 14) ){
      double	scaleX, scaleY, scaleZ;

      /* parse the record */
      numParsedFields = BibFileFieldParseFmt
	(bibfileRecord->field,
	 (void *) &scaleX, "%lg %*lg %*lg", "Scale",
	 (void *) &scaleY, "%*lg %lg %*lg", "Scale",
	 (void *) &scaleZ, "%*lg %*lg %lg", "Scale",
	 NULL);

      /* make the transform for this record */
      reconTrans = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
						 0.0, 0.0, 0.0, scaleX,
						 theta, 0.0, 0.0, 0.0,
						 0.0, 0, NULL);
    }
    BibFileRecordFree(&bibfileRecord);
    bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
  }

  if( outTrans ){
    if( absoluteFlg && reconTrans ){
      tmpTrans = WlzAffineTransformProduct(reconTrans, outTrans, NULL);
      WlzFreeAffineTransform(reconTrans);
      WlzFreeAffineTransform(outTrans);
      outTrans = tmpTrans;
    }

    switch( outputOpt ){
    default:
    case 2:
      errNum = WlzAffineTransformPrimGet(outTrans, &prim);
      if(errNum == WLZ_ERR_NONE)
      {
	fprintf(stdout, " -R -x%g -y%g -a%g -s%g", prim.tx, prim.ty,
		prim.theta, prim.scale);
      }
      break;

    case 1:
      domain.t = outTrans;
      values.core = NULL;
      outObj = WlzMakeMain(WLZ_AFFINE_TRANS, domain, values, NULL, NULL, NULL);
      WlzWriteObj(stdout, outObj);
      WlzFreeObj(outObj);
      break;

    case 0:
      errNum = WlzAffineTransformPrimGet(outTrans, &prim);
      if(errNum == WLZ_ERR_NONE)
      {
	fprintf(stdout,
		"Translation = (%f, %f)\n"
		"Rotation = %g\n"
		"Scale = %g\n",
		prim.tx, prim.ty, prim.theta, prim.scale);
      }
      break;
    }
  }

  return WLZ_ERR_NONE;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
