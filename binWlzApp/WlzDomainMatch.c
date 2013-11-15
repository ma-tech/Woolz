#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDomainMatch_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzDomainMatch.c
* \author       Richard Baldock
* \date         December 2005
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
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
* \ingroup      BinWlzApp
* \brief        Calculate a similarity value for a pair or set of woolz
* 		domains. This has been developed for the purposes of
* 		analysing the EMAGE gene-expression data.
*               
* \par Binary
* \ref wlzdomainmatch "WlzDomainMatch"
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzdomainmatch WlzDomainMatch
\par Name
WlzDomainMatch - calculate a match values between to domains or domain
 sets according to type.
\par Synopsis
\verbatim
WlzDomainMatch -d <delta> -t <type> -m <matrix-file> -s <scale> -h -v 

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-d</b></td>
    <td>delta value (default 0.01), must be < 1 </td>
  </tr>
  <tr>
    <td><b>-m \<file\></b></td>
    <td>input the name of a file containing the mixing
 and contribution matrices - csv format.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>type parameter to determine match function (default 1), values: </td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td>= 1 - Area(intersection)/Area(union)</td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td>= 2 - if Area(d1) > Area(d2) as type=1 else inverse</td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td>= 3 - Area(intersection)/Area(d1)</td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td>= 4 - Area(intersection)/Area(d2)</td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td>= 5 - Comparative match between two targets given by
                       the ratio of type 4 matchs to each domain. For this
                       option a 3rd object is required in the input stream.</td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td>= 6 - use the input mixing and contributing matrices. The matrix dimensions
 must match the number of categories in the category image data. Note the matrices
 do not need to be square.</td>
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
Read in domains or index images from stdin, calculate the match value
according to type and write the calculated value to stdout. The match functions are:
\f[
V_1 = \frac{S(d_1 \wedge d_2)}{S(d_1 \vee d_2)}
\f]
\f[
V_2 = \left \{
\begin{array}{r@{\quad if \quad}l}
V_1  &  S(d_1) \ge S(d_2) \\
\frac{1}{V_1} &  S(d_1) < S(d_2)
\end{array} \right.
\f]
\f[
V_3 = \frac{S(d_1 \wedge d_2)}{S(d_1)}
\f]
\f[
V_4 = \frac{S(d_1 \wedge d_2)}{S(d_2)}
\f]
\f[
V_5 =  \left \{
\begin{array}{c@{\quad : \quad}l}
1.0 & S(d_2) = 0 \quad\mbox{or}\quad  S(d_3) = 0, \\
1 / \delta & S(d_2) = 0 \quad\mbox{and}\quad  S(d_1 \wedge d_3) = 0, \\
\delta & S(d_3) = 0 \quad\mbox{and}\quad  S(d_1 \wedge d_2) = 0, \\
\frac{S(d_1 \wedge d_2)}{S(d_2)} \times \frac{S(d_3)}{S(d_1 \wedge d_3)}  &  \mbox{otherwise}. \\
\end{array} \right.
\f]
\f{eqnarray*}
V_6  & = & \frac{\sum_{Pixels} M_{ll'}}{A_{contrib}} \quad\mbox{where}\\
A_{contrib} & = & \sum_{Pixels} \left \{ \begin{array}{c@{\quad \mbox{if} \quad}l}
1 & C_{ll'} \ne 0,\\
0 & C_{ll'} = 0,
\end{array}\right.
\f}
In these formulae the \f$S()\f$ is the size of the domain - volume or area depending on the nature of the image. \f$l, l'\f$ are the pixel values of the two input category images.

\par Examples
\verbatim
\endverbatim

\par See Also
\par Bugs
None known
*/
#ifndef DOXYGEN_SHOULD_SKIP_THIS

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
	  "Usage: "
	  "%s -d <delta> -t <type> -m <matrix-file> -h -v \n"
	  "\tRead in domains from stdin and calculate the match value\n"
	  "\taccording to type, writing to stdout.\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t-d#        delta value (default 0.01), must be < 1\n"
	  "\t-m<file>   input the name of a file containing the mixing\n"
	  "\t           and contrib matrices - csv format\n"
	  "\t-t#        type parameter to determine match function default 1\n"
	  "\t             = 1 - Area(intersection)/Area(union)\n"
	  "\t             = 2 - if Area(d1) > Area(d2) as type=1 else inverse\n"
	  "\t             = 3 - Area(intersection)/Area(d1)\n"
	  "\t             = 4 - Area(intersection)/Area(d2)\n"
	  "\t             = 5 - Comparative match between two targets given by\n"
	  "\t               the ratio of type 4 matchs to each domain. For this\n"
	  "\t               option a 3rd object is required in the input stream.\n"
	  "\t-s #       scale factor to multiply the data values - useful for heatmap display.\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n",
	  str,
	  WlzVersion());

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

double WlzMixtureValue(
  WlzObject	*obj1,
  WlzObject	*obj2,
  int		numCatRows,
  int		numCatCols,
  double	**mixing,
  double	**contrib,
  WlzErrorNum	*dstErr)
{
  double	val = 0.0, con = 0.0;
  WlzObject	*tmpObj, *tmpObj1, *tmpObj2;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzGreyValueWSpace	*gVWSp;
  WlzPixelV		minP, maxP;
  int		i;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* objects must be same type with grey-values */
  if( (obj1 == NULL) || (obj2 == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj1->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj2->type != obj1->type ){
	errNum = WLZ_ERR_OBJECT_TYPE;
      }
      else {
	if((obj1->values.core == NULL) || (obj2->values.core == NULL)){
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else {
	  /* convert to int and check range */
	  if((tmpObj1 = WlzConvertPix(obj1, WLZ_GREY_INT, &errNum)) != NULL){
	    errNum = WlzGreyRange(tmpObj1, &minP, &maxP);
	    if( errNum == WLZ_ERR_NONE ){
	      if((minP.v.inv < 1) || (minP.v.inv > numCatRows) ||
		 (maxP.v.inv < 1) || (maxP.v.inv > numCatRows)){
		errNum = WLZ_ERR_OBJECT_DATA;
		WlzFreeObj(tmpObj1);
	      }
	    }
	    else {
	      WlzFreeObj(tmpObj1);
	    }
	  }
	  if((errNum == WLZ_ERR_NONE) &&
	     (tmpObj2 = WlzConvertPix(obj2, WLZ_GREY_INT, &errNum)) ){
	    errNum = WlzGreyRange(tmpObj2, &minP, &maxP);
	    if( errNum == WLZ_ERR_NONE ){
	      if((minP.v.inv < 1) || (minP.v.inv > numCatCols) ||
		 (maxP.v.inv < 1) || (maxP.v.inv > numCatCols)){
		errNum = WLZ_ERR_OBJECT_DATA;
		WlzFreeObj(tmpObj1);
		WlzFreeObj(tmpObj2);
	      }
	    }
	    else {
	      WlzFreeObj(tmpObj1);
	      WlzFreeObj(tmpObj2);
	    }
	  }
	}
      }
      break;

    case WLZ_3D_DOMAINOBJ:
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* get the intersection region */
  if( errNum == WLZ_ERR_NONE ){
    if((tmpObj = WlzIntersect2(tmpObj1, tmpObj2, &errNum)) != NULL){
      tmpObj->values = WlzAssignValues(tmpObj1->values, &errNum);
    }
    else {
      WlzFreeObj(tmpObj1);
      WlzFreeObj(tmpObj2);
    }
  }

  /* now calculate the mixture value */
  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(tmpObj, &iwsp, &gwsp);
    gVWSp = WlzGreyValueMakeWSp(tmpObj2, &errNum);
    val = 0.0;
    con = 0.0;
    while((errNum == WLZ_ERR_NONE) &&
	  (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){

      gptr = gwsp.u_grintptr;
      for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
	WlzGreyValueGet(gVWSp, 0, iwsp.linpos, iwsp.colpos+i);
	val += mixing[(*gptr.inp)-1][gVWSp->gVal[0].inv-1];
	con += contrib[(*gptr.inp)-1][gVWSp->gVal[0].inv-1];
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
    WlzGreyValueFreeWSp(gVWSp);
    WlzFreeObj(tmpObj);
    WlzFreeObj(tmpObj1);
    WlzFreeObj(tmpObj2);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  if( con > 0.0 ){
    return val/con;
  }
  else {
    return 0.0;
  }
}

int main(
  int   argc,
  char  **argv)
{
  FILE		*inFile;
  char 		optList[] = "d:m:s:t:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		verboseFlg=0;
  int		type=1;
  WlzObject	*obj = NULL, *obj1, *obj2, *obj3;
  double	matchVal=0.0;
  double	s1, s2, s3, s4;
  double	delta=0.01;
  double	**mixing=NULL, **contrib=NULL;
  double	scaleFactor = 1.0;
  int		k, l;
  int		numCatRows=-1, numCatCols=-1;

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

    case 'm':
      if((inFile = fopen(optarg, "r")) != NULL){
	if( fscanf(inFile, "%d, %d", &numCatCols, &numCatRows) < 2 ){
	  fprintf(stderr, "%s: can't read mixing matrix dimensions\n", argv[0]);
	  usage(argv[0]);
	  return 1;
	}
	AlcDouble2Malloc(&mixing, numCatRows, numCatCols);
	AlcDouble2Malloc(&contrib, numCatRows, numCatCols);
	for(l=0; l < numCatRows; l++){
	  for(k=0; k < numCatCols; k++){
	    if( fscanf(inFile, "%lg,", &(mixing[l][k])) < 1 ){
	      fprintf(stderr, "%s: can't read mixing matrix\n", argv[0]);
	      usage(argv[0]);
	      return 1;
	    }
	  }
	}
	for(l=0; l < numCatRows; l++){
	  for(k=0; k < numCatCols; k++){
	    if( fscanf(inFile, "%lg,", &(contrib[l][k])) < 1 ){
	      fprintf(stderr, "%s: can't read contributing matrix\n", argv[0]);
	      usage(argv[0]);
	      return 1;
	    }
	  }
	}
      }
      else {
	fprintf(stderr, "%s: can't open matrix file\n", argv[0]);
	usage(argv[0]);
	return 1;
      }
      break; 

    case 's':
      scaleFactor = (double) atof(optarg);
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

  /* verbose output */
  if( verboseFlg ){
    fprintf(stderr, "%s: parameter values:\n", argv[0]);
    fprintf(stderr, "\ttype = %d, delta = %f\n", type, delta);
    if( type == 6 ){
      fprintf(stderr, "\t mixing matrix:\n");
      for(l=0; l < numCatRows; l++){
	for(k=0; k < numCatCols; k++){
	  fprintf(stderr, "%f, ", mixing[l][k]);
	}
	fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
      fprintf(stderr, "\t contributing matrix:\n");
      for(l=0; l < numCatRows; l++){
	for(k=0; k < numCatCols; k++){
	  fprintf(stderr, "%f, ", contrib[l][k]);
	}
	fprintf(stderr, "\n");
      }
    }
  }

  /* get objects from stdin */
  inFile = stdin;
  if((obj1 = WlzReadObj(inFile, &errNum)) != NULL){
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

  if((obj2 = WlzReadObj(inFile, &errNum)) != NULL){
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
  if((obj3 = WlzReadObj(inFile, &errNum)) != NULL){
    if( (obj3->type != obj1->type) && (obj3->type != WLZ_EMPTY_OBJ) ){
      fprintf(stderr, "%s: objects must be same type\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }

  /* now calculate a match value */
  switch( type ){
  case 1:
    if((obj = WlzIntersect2(obj1, obj2, &errNum)) != NULL){
      s1 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s1 = 0;
    }
    if((obj = WlzUnion2(obj1, obj2, &errNum)) != NULL){
      s2 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s2 = 1;
    }
    matchVal = s1 / s2;
    break;

  case 2:
    if((obj = WlzIntersect2(obj1, obj2, &errNum)) != NULL){
      s1 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s1 = 0;
    }
    if((obj = WlzUnion2(obj1, obj2, &errNum)) != NULL){
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
    if((obj = WlzIntersect2(obj1, obj2, &errNum)) != NULL){
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
    if((obj = WlzIntersect2(obj1, obj2, &errNum)) != NULL){
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
    if((obj = WlzIntersect2(obj1, obj2, &errNum)) != NULL){
      s2 = WlzSize(obj, &errNum);
      WlzFreeObj(obj);
    }
    else {
      s2 = 0.0;
    }
    s3 = WlzSize(obj3, &errNum);
    if((obj = WlzIntersect2(obj1, obj3, &errNum)) != NULL){
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
    if( s2 > 0.0 ){
      if( s4 > 0.0 ){
	matchVal = (s2/s1) * (s3/s4);
      }
      else {
	matchVal = s2 / s1;
      }
    }
    else {
      if( s4 > 0.0 ){
	matchVal = s3/s4;
      }
      else {
	matchVal = 1.0;
      }
    }
    matchVal = WLZ_MAX(matchVal, delta);
    matchVal = WLZ_MIN(matchVal, 1.0/delta);
    break;

  case 6:
    /* this requires a mixing and contributing matrix and the images
       read in must have grey-values set to the right categories */
    if((numCatRows == -1) || (numCatCols == -1) ||
       (mixing == NULL) || (contrib == NULL)){
      fprintf(stderr, "%s: bad matrix data\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
    matchVal = WlzMixtureValue(obj1, obj2, numCatRows, numCatCols,
			       mixing, contrib, &errNum);
    break;

  default:
    fprintf(stderr, "%s: invalid match type\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* print value */
  fprintf(stdout, "%f\n", scaleFactor * matchVal);
  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
