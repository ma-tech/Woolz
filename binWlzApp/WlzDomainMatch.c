#pragma ident "MRC HGU $Id$"
/*!
\ingroup      BinWlzApp
\defgroup     wlzdomainmatch WlzDomainMatch
\par Name
WlzDomainMatch - 
\par Synopsis
\verbatim
WlzDomainMatch

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
\author       richard <Richard.Baldock@hgu.mrc.ac.uk>
\date         Fri Sep  9 08:34:25 2005
\version      MRC HGU $Id$
              $Revision$
              $Name$
\par Copyright:
             1994-2003 Medical Research Council, UK.
              All rights reserved.
\par Address:
              MRC Human Genetics Unit,
              Western General Hospital,
              Edinburgh, EH4 2XU, UK.
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
	  "%s -d <delta> -t <type> -h -v \n"
	  "\tRead in domains from stdin and calculate the match value\n"
	  "\taccording to type, writing to stdout.\n"
	  "Arguments:\n"
	  "\t-d#        delta value (default 0.01), must be < 1\n"
	  "\t-t#        type parameter to determine match function default 1\n"
	  "\t             = 1 - Area(intersection)/Area(union)\n"
	  "\t             = 2 - if Area(d1) > Area(d2) as type=1 else inverse\n"
	  "\t             = 3 - Area(intersection)/Area(d1)\n"
	  "\t             = 4 - Area(intersection)/Area(d2)\n"
	  "\t             = 5 - Comparative match between two targets given by\n"
	  "\t               the ratio of type 4 matchs to each domain. For this\n"
	  "\t               option a 3rd object is required in the input stream.\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n"
	  "\n",
	  str);

  return;
}

static int WlzSize(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  int	size;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( obj ){
    switch(obj->type){
    case WLZ_2D_DOMAINOBJ:
      size = WlzArea(obj, &errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
      size = WlzVolume(obj, &errNum);
      break;

    case WLZ_EMPTY_OBJ:
      size = 0;
      break;

    default:
      size = -1;
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    size = -1;
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return size;
}

int main(
  int   argc,
  char  **argv)
{
  FILE		*inFile=stdin;
  char 		optList[] = "t:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*errMsg;
  int		verboseFlg=0;
  int		type=1;
  WlzObject	*obj, *obj1, *obj2, *obj3;
  double	matchVal=0.0;
  double	s1, s2, s3, s4;
  double	delta=0.01;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'd':
      delta = atof(optarg);
      if((delta >= 1.0) || (delta <= 1.0e-10)){
	delta = 0.01;
	fprintf(stderr, "%s: invalid delta, reset to 0.01", argv[0]);
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
  if( obj1 = WlzReadObj(inFile, &errNum) ){
    switch( obj1->type ){
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_EMPTY_OBJ:
      break;

    default:
      fprintf(stderr, "%s: invalid object type: %d\n", argv[0],
	      obj->type);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: can't read first object\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  if( obj2 = WlzReadObj(inFile, &errNum) ){
    if( (obj2->type != obj1->type) && (obj2->type != WLZ_EMPTY_OBJ) ){
      fprintf(stderr, "%s: objects must be same type\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: can't read second object\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* this can fail silently but if there is an object it must
     be valid */
  if( obj3 = WlzReadObj(inFile, &errNum) ){
    if( (obj3->type != obj1->type) && (obj3->type != WLZ_EMPTY_OBJ) ){
      fprintf(stderr, "%s: objects must be same type\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }

  /* now calculate a match value */
  switch( type ){
  case 1:
    if( obj = WlzIntersect2(obj1, obj2, &errNum) ){
      s1 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s1 = 0;
    }
    if( obj = WlzUnion2(obj1, obj2, &errNum) ){
      s2 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s2 = 1;
    }
    matchVal = s1 / s2;
    break;

  case 2:
    if( obj = WlzIntersect2(obj1, obj2, &errNum) ){
      s1 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s1 = 0;
    }
    if( obj = WlzUnion2(obj1, obj2, &errNum) ){
      s2 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s2 = 1;
    }
    matchVal = s1 / s2;
    if( type == 2 ){
      s1 = WlzSize(obj1, &errNum);
      s2 = WlzSize(obj2, &errNum);
      if( s2 > s1 ){
	if( matchVal == 0.0 ){
	  matchVal = 10.0;
	}
	else {
	  matchVal = 1.0 / matchVal;
	  matchVal = WLZ_MIN(matchVal, 10.0);
	}
      }
    }
    break;

  case 3:
    if( obj = WlzIntersect2(obj1, obj2, &errNum) ){
      s1 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s1 = 0;
    }
    s2 = WlzSize(obj1, &errNum);
    matchVal = 0.0;
    if( s2 > 0 ){
      matchVal = s1 / s2;
    }
    break;

  case 4:
    if( obj = WlzIntersect2(obj1, obj2, &errNum) ){
      s1 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s1 = 0;
    }
    s2 = WlzSize(obj2, &errNum);
    matchVal = 0.0;
    if( s2 > 0 ){
      matchVal = s1 / s2;
    }
    break;

  case 5:
    /* this is a coparative measure designed to give a value of 1
       to a random pattern and between zero and infinite for
       matches to one or the other. For analysis for clustering
       probably better to use the logarithm. Zero and infinite are
       delta and 1/delta. */
    /* must be a third object */
    if( obj3 == NULL ){
      fprintf(stderr, "%s: for match option 5 3 input object required\n",
	      argv[0]);
      return 1;
    }
    s1 = WlzSize(obj2, &errNum);
    if( obj = WlzIntersect2(obj1, obj2, &errNum) ){
      s2 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s2 = 0.0;
    }
    s3 = WlzSize(obj3, &errNum);
    if( obj = WlzIntersect2(obj1, obj3, &errNum) ){
      s4 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s4 = 0.0;
    }
    if((s1 < 0.0) || (s2 < 0.0) || (s3 < 0.0) || (s4 < 0.0)){
      /* just fail */
      fprintf(stderr, "%s: something gone wrong, negative size.\n",
	      argv[0]);
      return 1;
    }
    /* calculating (s2/s1) * (s3/s4) */
    /* if the denominator non-zero then simple formula */
    if( (s1*s4) > 0.5 ){
      matchVal = (s2/s1) * (s3/s4);
    }
    /* special cases */
    else if( s1 < 0.5 ){
      /* if area(d1) is zero then so also is the intersection
	 with d1.
	 if in addition the intersection with d2 is non-zero then the
	 ratio is 1.0 (neutral)
	 if the intersection with d2 is zero then if area(2) is
	 zero the ration is 1.0 (neutral) else the ratio is
	 infinite.
      */
      matchVal = 1.0;
      if( s4 < 0.5 ){
	if( s3 > 0.5 ){
	  matchVal = 1.0/delta;
	}
      }
    }
    else if( s4 < 0.5 ){
      /* intersection with d1 is non-zero, if the intersection
	 with d2 is zero then if area(d2) is zero the ratio is 1.0
	 else infinite */
      matchVal = 1.0;
      if( s3 < 0.5 ){
	matchVal = 1.0/delta;
      }
    }
    matchVal = WLZ_MAX(matchVal, delta);
    matchVal = WLZ_MIN(matchVal, 1.0/delta);
    break;

  default:
    fprintf(stderr, "%s: invalid match type\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* print value */
  fprintf(stdout, "%f\n", matchVal);
  return 0;
}
