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
*   File       :   WlzTransposePlanes.c					*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Mon Jul 17 17:38:25 2000				*
*   $Revision$								*
*   $Name$								*
*   Synopsis    : Transpose the planes of an object with respect to a	*
*		given range (by default the first and last planes).	*
*		This is a reflection operation in 3D about the center	*
*		of the given range					*
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

static WlzObject *WlzTransposePlanes(WlzObject	*obj,
				     int	offset,
				     WlzErrorNum	*dstErr);

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-o#] [-h] [-v] [<input file>]\n"
	  "\tTranspose planes in the 3D object with respect\n"
	  "\tto the offset so that p' = offset - p\n"
	  "\tThis is equivalent to reflection about offset/2\n"
	  "\tOptions are:\n"
	  "\t  -o#       offset(default plane1+lastpl)\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  int		option;
  char 		optList[] = "ho:v";
  int		offset;
  int		offsetFlg;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		verboseFlg=0;
    
  /* set the defaults */
  offsetFlg = 0;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'o':
      offset = atoi(optarg);
      offsetFlg = 1;
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_PARAM_TYPE;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_ERR_PARAM_TYPE;
    }
  }

  /* read objects and select planes if possible */
  while((obj = WlzReadObj(inFile, &errNum)) != NULL) 
  {
    if( errNum == WLZ_ERR_NONE ){
      switch( obj->type )
      {
      case WLZ_3D_DOMAINOBJ:
	if( obj->domain.p->type == WLZ_EMPTY_OBJ ){
	  errNum = WlzWriteObj(stdout, obj);
	  break;
	}
	/* check the transpose parameters */
	if( !offsetFlg ){
	  offset = obj->domain.p->plane1 + obj->domain.p->lastpl;
	}

	if( nobj = WlzTransposePlanes(obj, offset, &errNum) ){
	  errNum = WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
	break;

      default:
	errNum = WlzWriteObj(stdout, obj);
	break;
      }
    }

    WlzFreeObj(obj);
  }

  return errNum;
}

WlzObject *WlzTransposePlanes(
  WlzObject	*obj,
  int		offset,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  WlzDomain	domain;
  WlzValues	values;
  int		indx, indx1, p;

  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	 errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else {
	switch( obj->domain.p->type ){
	case WLZ_EMPTY_OBJ:
	case WLZ_EMPTY_DOMAIN:
	  return WlzMakeEmpty(dstErr);

	default:
	  domain.p = WlzMakePlaneDomain(obj->domain.p->type,
					offset-obj->domain.p->lastpl,
					offset-obj->domain.p->plane1,
					obj->domain.p->line1,
					obj->domain.p->lastln,
					obj->domain.p->kol1,
					obj->domain.p->lastkl,
					&errNum);
	  if( errNum == WLZ_ERR_NONE ){
	    domain.p->voxel_size[0] = obj->domain.p->voxel_size[0];
	    domain.p->voxel_size[1] = obj->domain.p->voxel_size[1];
	    domain.p->voxel_size[2] = obj->domain.p->voxel_size[2];
	    if( obj->values.core ){
	      values.vox =
		WlzMakeVoxelValueTb(obj->values.vox->type,
				    offset-obj->domain.p->lastpl,
				    offset-obj->domain.p->plane1,
				    obj->values.vox->bckgrnd,
				    NULL, &errNum);
	    }
	    else {
	      values.core = NULL;
	    }
	  }
	  break;
	}
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* make links */
  if( errNum == WLZ_ERR_NONE ){
    indx = obj->domain.p->lastpl - obj->domain.p->plane1;
    indx1 = 0;
    for(p=domain.p->plane1; p <= domain.p->lastpl; p++, indx--, indx1++){
      if( obj->domain.p->domains[indx].core ){
	domain.p->domains[indx1] =
	  WlzAssignDomain(obj->domain.p->domains[indx], &errNum);
      }
      else {
	domain.p->domains[indx1].core = NULL;
      }

      if( values.core ){
	if( obj->values.vox->values[indx].core ){
	  values.vox->values[indx1] =
	    WlzAssignValues(obj->values.vox->values[indx], &errNum);
	}
	else {
	  values.vox->values[indx1].core = NULL;
	}
      }
    }
  }

  /* create object */
  if( errNum == WLZ_ERR_NONE ){
    rtnObj = WlzMakeMain(obj->type, domain, values,
			 NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
