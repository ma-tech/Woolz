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
*   File       :   WlzVtxToObj.c					*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Thu Nov  9 13:46:16 2000				*
*   $Revision$								*
*   $Name$								*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
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
	  "%s\n"
	  "\tBuild a woolz domain object from a vertex\n"
	  "\n",
	  str);

  return;
}

WlzObject *WlzCoordsToObject(
  WlzObjectType	type,
  double	x,
  double	y,
  double	z,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL, *tmpObj;
  WlzDomain	domain;
  WlzValues	values;
  int		k, l, p;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check type */
  switch( type ){

  case WLZ_2D_DOMAINOBJ:
    k = x;
    l = y;
    if( domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					 l, l, k, k, &errNum) ){
      values.core = NULL;
      rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			   NULL, NULL, &errNum);
    }
    break;

  case WLZ_3D_DOMAINOBJ:
    if( tmpObj = WlzCoordsToObject(WLZ_2D_DOMAINOBJ,
				   x, y, z, &errNum) ){
      k = x;
      l = y;
      p = z;
      if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					p, p, l, l, k, k, &errNum) ){
	values.core = NULL;
	domain.p->domains[0] = WlzAssignDomain(tmpObj->domain, NULL);
	rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			     NULL, NULL, &errNum);
      }
      WlzFreeObj(tmpObj);
    }
    break;

  default:
    errNum = WLZ_ERR_OBJECT_TYPE;
    break;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

int main(
  int   argc,
  char  **argv)
{
  
  WlzObject     	*obj;
  WlzDVertex3		vtx;
  char 	optList[] = "hv";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* parse the command line */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* read vertex from stdin and write the objects to stdout */
  while( scanf("%lg %lg %lg", &vtx.vtX, &vtx.vtY, &vtx.vtZ) != EOF ){
    obj = WlzCoordsToObject(WLZ_2D_DOMAINOBJ, vtx.vtX, vtx.vtY, vtx.vtZ,
			    &errNum);
    WlzWriteObj(stdout, obj);
    WlzFreeObj(obj);
    obj = WlzCoordsToObject(WLZ_3D_DOMAINOBJ, vtx.vtX, vtx.vtY, vtx.vtZ,
			    &errNum);
    WlzWriteObj(stdout, obj);
    WlzFreeObj(obj);
  }

  return 0;
}
