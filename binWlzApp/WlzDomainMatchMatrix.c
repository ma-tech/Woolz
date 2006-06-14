#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzDomainMatchMatric_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzDomainMatchMatrix.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Jun 14 07:50:43 BST 2006
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par Copyright:
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
* \ingroup      BinWlzApp
* \brief        Calculate a similarity value for a pair or set of woolz domains.
 This has been developed for the purposes of analysing the EMAGE gene-expression data.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzdomainmatchmatrix WlzDomainMatchMatrix
\par Name
WlzDomainMatchMatrix - calculate match values between to domain
 sets according to type.
\par Synopsis
\verbatim
WlzDomainMatchMatrix -d <delta> -t <type> -m <matrix-file> -h -v <rows> <cols>

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
                       two objects per column are required in the input stream.</td>
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
	  "Usage:\n"
	  "%s -d <delta> -t <type> -m <matrix-file> -h -v  <rows> <cols>\n"
	  "\tRead in domains from stdin and calculate the match value\n"
	  "\tmatrix according to type, writing to stdout. The row domains\n"
	  "\tare read firrst followed by the column domains, there must be\n"
	  "\tsufficient for the match type.\n"
	  "Arguments:\n"
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
	  "\t               option two objects are required per column in the input.\n"
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

double WlzMixtureValue(
  WlzObject	*obj1,
  WlzObject	*obj2,
  int		numCatRows,
  int		numCatCols,
  double	**mixing,
  double	**contrib,
  WlzErrorNum	*dstErr)
{
  double	val, con;
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
	  if( tmpObj1 = WlzConvertPix(obj1, WLZ_GREY_INT, &errNum) ){
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
    if( tmpObj = WlzIntersect2(tmpObj1, tmpObj2, &errNum) ){
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
  char 		optList[] = "d:m:t:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*errMsg;
  int		verboseFlg=0;
  int		type=1;
  int		numRows=0, numCols=0;
  WlzObject	*obj, *obj1, *obj2, *obj3;
  WlzObject	**rowDoms, **colDoms;
  WlzObjectType	objType;
  double	matchVal=0.0;
  double	s1, s2, s3, s4;
  double	delta=0.01;
  double	**mixing=NULL, **contrib=NULL;
  int		i, j, k, l;
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
      if( inFile = fopen(optarg, "r") ){
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

  /* there must be two more arguments */
  if( (argc - optind) < 2 ){
    fprintf(stderr, "%s: not enough arguments\n", argv[0]);
    usage(argv[0]);
    return 1;
  }
  numRows = atoi(*(argv+optind));
  numCols = atoi(*(argv+optind+1));
  if((numRows <= 0) || (numCols <= 0)){
    fprintf(stderr, "%s: both number of rows and columns must be > 0\n", argv[0]);
    usage(argv[0]);
    return 1;
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
  rowDoms = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * numRows);
  for(i=0; i < numRows; i++){
    if( rowDoms[i] =  WlzReadObj(inFile, &errNum) ){
    if((i == 0) || (objType == WLZ_EMPTY_OBJ)){
      objType = rowDoms[i]->type;
    }
    if((rowDoms[i]->type != objType) &&
       (rowDoms[i]->type != WLZ_EMPTY_OBJ)){
      fprintf(stderr, "%s: invalid object type: %d\n", argv[0],
	      rowDoms[i]->type);
      usage(argv[0]);
      return 1;
    }
    }
    else {
      fprintf(stderr, "%s: not enough row objects\n", argv[0]);
    }
  }

  /* now the column domains */
  if( type == 5 ){
    colDoms = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * numCols * 2);
    for(i=0; i < numCols*2; i++){
      if( colDoms[i] =  WlzReadObj(inFile, &errNum) ){
	if((colDoms[i]->type != objType) &&
	   (colDoms[i]->type != WLZ_EMPTY_OBJ)){
	  fprintf(stderr, "%s: invalid object type: %d\n", argv[0],
		  colDoms[i]->type);
	  usage(argv[0]);
	  return 1;
	}
      }
      else {
	fprintf(stderr, "%s: not enough row objects\n", argv[0]);
      }
    }
  }
  else {
    colDoms = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * numCols);
    for(i=0; i < numCols; i++){
      if( colDoms[i] =  WlzReadObj(inFile, &errNum) ){
	if((colDoms[i]->type != objType) &&
	   (colDoms[i]->type != WLZ_EMPTY_OBJ)){
	  fprintf(stderr, "%s: invalid object type: %d\n", argv[0],
		  colDoms[i]->type);
	  usage(argv[0]);
	  return 1;
	}
      }
      else {
	fprintf(stderr, "%s: not enough row objects\n", argv[0]);
      }
    }
  }

  /* now calculate the match values */
  for(i=0; i < numRows; i ++){
    if( verboseFlg ){
      fprintf(stderr, "%s: start row %d\n", argv[0], i+1);
    }
    obj1 = rowDoms[i];
    for(j=0; j < numCols; j++){
      if( type == 5 ){
	obj2 = colDoms[j*2];
	obj3 = colDoms[j*2 + 1];
      }
      else {
	obj2 = colDoms[j];
      }
      
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
      fprintf(stdout, "%f", matchVal);
      if( (numCols - j) > 1 ){
	fprintf(stdout, "\t");
      }
      if( verboseFlg ){
	fprintf(stderr, ".");
      }
    }
    fprintf(stdout, "\n");
    if( verboseFlg ){
      fprintf(stderr, "\n%s: completed row %d\n", argv[0], i+1);
    }
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
