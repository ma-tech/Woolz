#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        WlzExtractTransform.c
* Date:         April 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Finds the affine transforms from a reconstruction
*		bibfile.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

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
	  "Usage:\t%s [-A] [-R] [-a] [-t] [-T] [-ffilename] [input bibfile]"
	  "[-h] [-v] [<bibfile>]\n"
	  "\tFind the transform for the given filename\n"
	  "\twithin the bib-file and output an argument\n"
	  "\tstring for WlzAffineTransformObj\n"
	  "\tOptions are:\n"
	  "\t  -A                calculate absolute transform\n"
	  "\t  -R                calculate relative transform (default)\n"
	  "\t  -a                output and argument string for\n"
	  "\t                    WlzAffineTransformObj (default)\n"
	  "\t  -t                output a woolz affine transform object\n"
	  "\t  -T                output text\n"
	  "\t  -f<filename>      filename to be found\n"
	  "\t  -h                Help - prints this usage message\n"
	  "\t  -v                verbose operation\n"
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
  BibFileField	*bibfileField;
  BibFileError	bibFileErr=BIBFILE_ER_NONE;
  int		outputOpt=2;
  int		absoluteFlg=0;
  WlzAffineTransform	*inTrans, *outTrans, *tmpTrans;
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
  outTrans = WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE,
					0.0, 0.0, 0.0, 1.0,
					0.0, 0.0, 0.0, 0.0, 0.0, 0, NULL);

  /* read bibfile records until the filename matches */
  bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
  while( bibFileErr == BIBFILE_ER_NONE ) 
  {
    if( strncmp(bibfileRecord->name, "Section", 7) ){
      BibFileRecordFree(&bibfileRecord);
    }
    else {
      /* parse the record */
      numParsedFields = BibFileFieldParseFmt
	(bibfileRecord->field,
	 (void *) &tx, "%lg", "TransformTx",
	 (void *) &ty, "%lg", "TransformTy",
	 (void *) &theta, "%lg", "TransformTheta",
	 (void *) &recordFilename, "%s", "File",
	 NULL);

      /* make the transform for this record */
      inTrans = WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE,
					   tx, ty, 0.0, 1.0,
					   theta, 0.0, 0.0, 0.0, 0.0, 0, NULL);

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

      BibFileRecordFree(&bibfileRecord);
    }
    bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
  }

  if( outTrans ){
    switch( outputOpt ){
    default:
    case 2:
      fprintf(stdout, " -R -x%g -y%g -a%g", outTrans->tx, outTrans->ty,
	      outTrans->theta);
      break;

    case 1:
      domain.t = outTrans;
      values.core = NULL;
      outObj = WlzMakeMain(WLZ_AFFINE_TRANS, domain, values, NULL, NULL, NULL);
      WlzWriteObj(stdout, outObj);
      WlzFreeObj(outObj);
      break;

    case 0:
      fprintf(stdout,
	      "Translation = (%f, %f)\n"
	      "Rotation = %g\n",
	      outTrans->tx, outTrans->ty, outTrans->theta);
      break;
    }
  }

  return WLZ_ERR_NONE;
}
