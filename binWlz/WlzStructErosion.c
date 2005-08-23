#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   Woolz Library					*
*   File       :   WlzStructErosion.c					*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Mon May  7 15:36:24 2001				*
*   $Revision$								*
*   $Name$								*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

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
	  "Usage:\t%s [-s <struct-elem file>] [-h] [<input file>]\n"
	  "\tErode a domain woolz object(s) using given\n"
	  "\tstructuring element writing the new object(s)\n"
	  "\tto standard output\n"
	  "\tOptions are:\n"
	  "\t  -s <file> file holding the structuring\n"
	  "\t            element object. Default is the\n"
	  "\t            first object read in from the input\n"
	  "\t            file or standard input.\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *seObj, *nobj;
  FILE		*inFile;
  FILE		*seFile=NULL;
  char 		optList[] = "hs:";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 's':
      if( (seFile = fopen(optarg, "r")) == NULL ){
	fprintf(stderr, "%s: can't open structuring element file: %s\n",
		argv[0], optarg);
	usage(argv[0]);
	return WLZ_ERR_UNSPECIFIED;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }
  if( !seFile ){
    seFile = inFile;
  }
  
  if( seObj = WlzAssignObject(WlzReadObj(seFile, &errNum), NULL) ){

    /* read objects and threshold if possible */
    while((obj = WlzAssignObject(WlzReadObj(inFile, &errNum), NULL))
	  != NULL) 
    {
      if( nobj = WlzStructErosion(obj, seObj, &errNum) ){
	nobj = WlzAssignObject(nobj, NULL);
	WlzWriteObj(stdout, nobj);
	WlzFreeObj(nobj);
      }
      WlzFreeObj(obj);
    }
    WlzFreeObj(seObj);
  }
  
  /* trap the WLZ_ERR_READ_EOF since this is a legal way of indicating
     the end of objects in a file */
  if( errNum == WLZ_ERR_READ_EOF ){
    errNum = WLZ_ERR_NONE;
  }
  return errNum;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
