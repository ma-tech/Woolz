#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreySetValue.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Sets the grey values of the object to the input value.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
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
	  "Usage:\t%s [-g#] [-h] [-v] [<input mask> [<input obj>]]\n"
	  "\tSet the grey values of the object to the input value.\n"
	  "\tA valuetable will be attached if required.\n"
	  "\tOptions are:\n"
	  "\t  -g#       the new grey value - default 0\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *tmpObj;
  WlzDomain	*domains;
  WlzValues	values, *valuess;
  FILE		*inFile;
  char 		optList[] = "g:hv";
  int		option;
  WlzPixelV	greyVal;
  WlzPixelV	bckgrnd;
  WlzObjectType	type;
  int		verboseFlg=0, p;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  greyVal.type = WLZ_GREY_FLOAT;
  greyVal.v.flv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'g':
      greyVal.v.flv = atof(optarg);
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

  /* check for read from file */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* set up type and background */
  type = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_UBYTE, NULL);
  bckgrnd.type = WLZ_GREY_UBYTE;
  bckgrnd.v.ubv = 0;

  /* read objects and set value if possible */
  while(((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) &&
        (errNum == WLZ_ERR_NONE))
  {
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
      if( obj->values.core == NULL ){
	values.v = WlzNewValueTb(obj, type, bckgrnd, &errNum);
	obj->values = WlzAssignValues(values, NULL);
      }
      errNum = WlzGreySetValue(obj, greyVal);
      break;

    case WLZ_3D_DOMAINOBJ:
      switch( obj->domain.p->type ){
      case WLZ_PLANEDOMAIN_DOMAIN:
	if( obj->values.core == NULL ){
	  if( values.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					       obj->domain.p->plane1,
					       obj->domain.p->lastpl,
					       bckgrnd, NULL, &errNum) ){
	    obj->values = WlzAssignValues(values, NULL);
	    domains = obj->domain.p->domains;
	    valuess = obj->values.vox->values;
	    for(p=0; p < (obj->domain.p->lastpl - obj->domain.p->plane1 + 1);
		p++, domains++, valuess++){
	      if( (*domains).core ){
		values.core = NULL;
		tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, values,
				     NULL, NULL, NULL);
		values.v = WlzNewValueTb(tmpObj, type, bckgrnd, &errNum);
		*valuess = WlzAssignValues(values, NULL);
		WlzFreeObj(tmpObj);
	      }
	    }
	  }
	}
	errNum = WlzGreySetValue(obj, greyVal);
	break;

      default:
	break;
      }
      break;

    default:
      break;
    }

    if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write object (%s).\n",
		     argv[0], errMsg);
      return(1);
    }
    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
