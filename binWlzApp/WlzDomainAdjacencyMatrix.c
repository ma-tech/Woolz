#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDomainAdjacencyMatrix_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzDomainAdjacencyMatrix.c
* \author       Richard Baldock
* \date         November 2012
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
* \brief        Calculate the adjacency profiles for a given domain
* 		with respect to a set of domains. Each profile is
* 		output as a row in a matrix. The first row in the
*		matrix indicates the adjacency distance.
*               
* \par Binary
* \ref wlzdomainmatchmatrix "WlzDomainAdjacencyMatrix"
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzdomainadjacencymatrix WlzDomainAdjacencyMatrix
\par Name
WlzDomainAdjacencyMatrix - calculate adjacency profiles for a domain
 with respect to a domain set.
\par Synopsis
\verbatim
WlzDomainAdjacencyMatrix [-a min,max,step] [-f <domain file>] [-n] [-t] [-T][-h] [-v] <test domain file>

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>min and max adjacency distance and calculation step (default -10, 100, 1)</td>
  </tr>
  <tr>
    <td><b>-f \<file\></b></td>
    <td>input the name of a file with the target domain set</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>normalise values by target volume</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Output progress ticks</td>
  </tr>
  <tr>
    <td><b>-T</b></td>
    <td>Output tab-separated values (TSV) default is CSV</td>
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
Read in domains from the input file or from stdin, calculate the adjacency profile
with respect to the test domain and output as a matrix of values. The first row are
the adjacency distances used to calculate each column.
The adjacency value for distance \f$r\f$ is defined as:
\f[
Adjacency (r) = Volume( Dilation(test domain, r) \wedge (target domain) )
\f]

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
	  "%s [-a min,max,step] [-f <domain file>] [-n <norm data file>] [-t] [-T] [-h] [-v]  <test domain file>\n"
	  "\tRead in domains from the input file or from stdin, calculate the \n"
	  "\tadjacency profile with repect to the test domain and output as a \n"
	  "\tmatrix of values. The first row is the adjacency distances used to\n"
	  "\tcalculate each column.\n"
	  "Woolz Version: %s\n"
	  "Arguments:\n"
	  "\t-a#,#,#    min and max values of the adjacency distance and step (default -10,100,1)\n"
	  "\t-f <file>  input the name of a file containing the target domain set.\n"
	  "\t-d <file>  calculate derivative and output to file given\n"
	  "\t-n <file>  normalise by target volume (so max val = 1.0) and output to file given\n"
	  "\t-t         output progress ticks\n"
	  "\t-T         output TSV, default CSV\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n"
	  "\n",
	  str, WlzVersion());

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
  FILE		*infileDomains = NULL, *infile = NULL, *outfile = NULL;
  FILE		*normOutfile = NULL, *derivOutfile = NULL;
  char 		optList[] = "a:d:f:n:tThv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		verboseFlg=0, tickFlg=0;
  int		normaliseFlg=0, derivFlg=0, TSVFlg=0;
  char		separatorChar=',';
  WlzObject	**domains;
  int		numDomains;
  WlzObject	*obj, *obj1, *testDomain;
  WlzObjectType	objType;
  int		dim=0;
  WlzObject	**structElements, **refDomains;
  int		i, j, radius;
  int		rMin=-10, rMax=100, rStep=1;
  int       	adjVal;
  double	normalisedAdjVal = 0.0;

  /* read the argument list and check for an input file */
  infile  = stdin;
  infileDomains = stdin;
  outfile = stdout;
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'a':
      switch( sscanf(optarg, "%d,%d,%d", &rMin, &rMax, &rStep) ){
      case 0: /* set default values */
	fprintf(stderr, "%s: setting default range and step: %d,%d,%d\n",
		argv[0], rMin, rMax, rStep);
	break;
      case 1:
	fprintf(stderr, "%s: setting default range max and step: %d,%d\n",
		argv[0], rMax, rStep);
	break;
      case 2:
	fprintf(stderr, "%s: setting default range step: %d\n",
		argv[0], rStep);
	break;
      default:
      case 3:
	break;
      }
      break;

    case 'd':
      derivFlg = 1;
      if((derivOutfile = fopen(optarg, "r")) == NULL){
	fprintf(stderr, "%s: can't open output file for derivative values %s\n", argv[0], optarg);
	usage(argv[0]);
	return 1;
      }
      break;

    case 'f':
      if((infileDomains = fopen(optarg, "r")) == NULL){
	fprintf(stderr, "%s: can't open domains file %s\n", argv[0], optarg);
	usage(argv[0]);
	return 1;
      }
      break; 

    case 'n':
      normaliseFlg = 1;
      if((normOutfile = fopen(optarg, "w")) == NULL){
	fprintf(stderr, "%s: can't open output file for normalised values %s\n", argv[0], optarg);
	usage(argv[0]);
	return 1;
      }
      break;

    case 't':
      tickFlg = 1;
      break;

    case 'T':
      TSVFlg = 1;
      separatorChar = '\t';
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

  /* Check for test domain file */
  if( (argc - optind) > 0 ){
    if((infile = fopen(*(argv + optind), "r")) == NULL){
      fprintf(stderr, "%s: can't open domains file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* verbose output */
  if( verboseFlg ){
    fprintf(stderr, "%s: calculation of adjacency profiles\n", argv[0]);
    fprintf(stderr, "adjacency range: %d to %d step %d\n", rMin, rMax, rStep);
  }

  /* get test object from infile */
  if( (testDomain = WlzReadObj(infile, &errNum)) == NULL ){
    fprintf(stderr, "%s: can't read test domain object.\n", argv[0]);
    usage(argv[0]);
    return 1;
  }
  switch( testDomain->type ){
  case WLZ_2D_DOMAINOBJ:
    objType = testDomain->type;
    dim = 2;
    break;

  case WLZ_3D_DOMAINOBJ:
    objType = testDomain->type;
    dim = 3;
    break;

  case WLZ_TRANS_OBJ: /* some types to support later */
  case WLZ_2D_POLYGON:
  case WLZ_3D_POLYGON:
  default:
    usage(argv[0]);
    return 1;
  }

  /* read the domain set - must all have same type as test domain or empty */
  domains = NULL;
  numDomains = 0;
  while( (obj = WlzReadObj(infileDomains, &errNum)) ){
    if( (obj->type != objType) && (obj->type != WLZ_EMPTY_OBJ) ){
      fprintf(stderr, "%s: found domain of invalid type\n", argv[0]);
      errNum = WlzFreeObj(obj);
      continue;
    }
    if( (numDomains % 10) == 0 ){
      domains = (WlzObject **) realloc(domains, sizeof(WlzObject *) * (numDomains+10));
    }
    domains[numDomains] = WlzAssignObject(obj, &errNum);
    numDomains++;
  }
  if( verboseFlg ){
    fprintf(stderr, " %s: num domains read: %d\n", argv[0], numDomains);
  }

  /* now calculate the adjacency profiles */
  /* each row is the profile of the test domain against one from the list */
  /* the first row is the set of radii for which the intersection volume is calculated */
  /* output radii and calculate structuring element */
  structElements = (WlzObject **) malloc(sizeof(WlzObject *) * ((rMax - rMin)/rStep + 2));
  refDomains = (WlzObject **) malloc(sizeof(WlzObject *) * ((rMax - rMin)/rStep + 2));
  if( tickFlg ) {
    fprintf(stderr, "Make structuring elements and dilate testdomain:");
  }
  for( radius=rMin, i=0; radius <= rMax; radius += rStep, i++){

    /* output the radius */
    if( (radius + rStep) > rMax ){
      fprintf(outfile, "%d\n", radius);
    } else {
      fprintf(outfile, "%d%c", radius, separatorChar);
    }
    if( normaliseFlg ){
      if( (radius + rStep) > rMax ){
	fprintf(normOutfile, "%d\n", radius);
      } else {
	fprintf(normOutfile, "%d%c", radius, separatorChar);
      }
    }

    /* create the structuring element for the erosion/dilation and dilate the reference domain */
    if( radius < 0 ){
      switch( dim ){
      case 2:
	obj = WlzMakeCircleObject((double) -radius, 0, 0, &errNum);
	break;
      case 3:
	obj = WlzMakeSphereObject(WLZ_3D_DOMAINOBJ, (double) -radius, 0, 0, 0, &errNum);
	break;
      }
      structElements[i] = WlzAssignObject(obj, &errNum);
      if( (obj = WlzStructErosion( testDomain, structElements[i], &errNum)) ){
	refDomains[i] = WlzAssignObject(obj, &errNum);
      }
      else {
	refDomains[i] = NULL;
      }
    }
    else if (radius > 0 ){
      switch( dim ){
      case 2:
	obj = WlzMakeCircleObject((double) radius, 0, 0, &errNum);
	break;
      case 3:
	obj = WlzMakeSphereObject(WLZ_3D_DOMAINOBJ, (double) radius, 0, 0, 0, &errNum);
	break;
      }
      structElements[i] = WlzAssignObject(obj, &errNum);
      if( (obj = WlzStructDilation( testDomain, structElements[i], &errNum)) ){
	refDomains[i] = WlzAssignObject(obj, &errNum);
      }
      else {
	refDomains[i] = NULL;
      }
    }
    else {
      structElements[i] = NULL;
      refDomains[i] = WlzAssignObject(testDomain, &errNum);
    }

    /* tickFlg to give feedback of progress */
    if( tickFlg ) {
      fprintf(stderr, ".");
    }
  }
  if( tickFlg ){
    fprintf(stderr, "done\n");
  }

  /* now calculate intersection volume */
  for( j=0; j < numDomains; j++){
    if( tickFlg ) {
      fprintf(stderr, "Calc Adjacencies index=%d:", j);
    }
    for( radius=rMin, i=0; radius <= rMax; radius += rStep, i++){
      adjVal = 0;
      if( refDomains[i] != NULL ){
	if( (obj1 = WlzIntersect2(refDomains[i], domains[j], &errNum)) ){
	  adjVal = WlzSize(obj1, &errNum);
	  errNum = WlzFreeObj(obj1);
	}
	if( normaliseFlg ){
	  normalisedAdjVal = ((double) adjVal) / ((double) WlzSize(domains[j], &errNum));
	}
      }	

      if( verboseFlg ){
	fprintf(stderr, "%s: radius=%d, domainIndx=%d, adjVal=%d\n", argv[0], radius, j, adjVal);
	if( normaliseFlg ){
	  fprintf(stderr, "%s: radius=%d, domainIndx=%d, normalised adjVal=%.6f\n", argv[0], radius, j, normalisedAdjVal);
	}
      }

      if ( (radius + rStep) > rMax ){
	fprintf(outfile, "%d\n", adjVal);
      } else {
	fprintf(outfile, "%d%c", adjVal, separatorChar);
      }
      if( normaliseFlg ){
	if ( (radius + rStep) > rMax ){
	  fprintf(normOutfile, "%.6f\n", normalisedAdjVal);
	} else {
	  fprintf(normOutfile, "%.6f%c", normalisedAdjVal, separatorChar);
	}
      }
      if( tickFlg ){
	fprintf(stderr, ".");
      }

    }

    if( tickFlg ){
      fprintf(stderr, "completed index\n");
    }

    if( verboseFlg ){
      fprintf(stderr, "%s: completed row %d\n\n", argv[0], i+1);
    }
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
