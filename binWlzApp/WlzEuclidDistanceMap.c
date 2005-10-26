#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzEuclidDistanceMap.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Tue Oct 18 17:36:31 2005
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
* \ingroup      BinWlzApp
* \brief        Calculate the Euclidean distance from the input
 vertex to all parts of the domain.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

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
	  "%s -p x,y,z  -h -v \n"
	  "\tRead in domain object from stdin and calculate the Euclidean\n"
	  "\tdistance from the input vertex to each point in the domain,\n"
	  "\twriting to stdout.\n"
	  "Arguments:\n"
	  "\t-p#,#,#    input vertex, default (0,0,0)\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n"
	  "\n",
	  str);

  return;
}

int main(
  int   argc,
  char  **argv)
{
  FILE		*inFile=stdin;
  char 		optList[] = "p:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*errMsg;
  int		verboseFlg=0;
  int		type;
  WlzObject	*obj, *obj1;
  int		i;
  double	x, y, z, val;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;

  /* read the argument list and check for an input file */
  opterr = 0;
  x = y = z = 0.0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'p':
      if( sscanf(optarg, "%lf,%lf,%lf", &x, &y, &z) != 3 ){
	fprintf(stderr, "%s: not enough values for vertex\n", argv[0]);
	usage(argv[0]);
	return 1;
      }
      break;

    case 't':
      type = atoi(optarg);
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

  /* get objects from stdin */
  if( obj = WlzReadObj(inFile, &errNum) ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      break;

    case WLZ_3D_DOMAINOBJ: /* to be done */
    default:
      fprintf(stderr, "%s: invalid object type: %d\n", argv[0],
	      obj->type);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: can't read object\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* scan through object setting distance values
     use integer values */
  if( obj1 = WlzConvertPix(obj, WLZ_GREY_INT, &errNum) ){
    WlzInitGreyScan(obj1, &iwsp, &gwsp);
    while( WlzNextGreyInterval(&iwsp) == 0 ){

      gptr = gwsp.u_grintptr;
      switch (gwsp.pixeltype) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
	  val = (iwsp.colpos + i - x)*(iwsp.colpos + i - x) \
	    + (iwsp.linpos - y) * (iwsp.linpos - y);
	  val = sqrt(val);
	  *gptr.inp = WLZ_NINT(val);
	}
	break;

      default:
	/* something gone badly wrong - quit */
	fprintf(stderr, "%s: something wrong, oops\n", argv[0]);
	return 1;
      }
    }
    
    WlzWriteObj(stdout, obj1);
    WlzFreeObj(obj1);
  }

  WlzFreeObj(obj);
  return 0;
}
