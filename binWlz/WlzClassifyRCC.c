#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzClassifyRCC_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzClassifyRCC.c
* \author       Bill Hill
* \date         March 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Classifies an ordered pair of domain objects using a
* 		region connected calculus.
* \ingroup	BinWlz
* \par Binary
* \ref wlzclassifyrcc   "WlzClassifyRCC"
*/

/*!
\ingroup BinWlz
\defgroup wlzclassifyrcc   WlzClassifyRCC
\par Name
WlzClassifyRCC - Classifies the given domain objects using the region
calculus defined in Woolz.
\par Synopsis
\verbatim
WlzClassifyRCC [-o<output file>] [-a] [-e] [-h] [<in object 0>] [<in object 1>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>Don't include adjacency classifications.</td>
  </tr>
  <tr>
    <td><b>-e</b></td>
    <td>Don't include enclosure classifications.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Maximum adjacency distance in piels/voxels,
        this should be a resonable distance since domains will be
	dilated using this figure.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints message.</td>
  </tr>
</table>
\par Description
Classifies the given domain objects using the region
calculus defined in Woolz.
\par Examples
The domain of the Woolz object read from the file somites.wlz
is classified with respect to the domain of the Woolz object read
from the file embryonic.wlz.
The classification is that somites is a tangential proper part of
and enclosed by embryonic is written to the standard output.
\verbatim
WlzClassifyRCC somites.wlz embryonic.wlz
WLZ_RCC_TPP|WLZ_RCC_ENC        0.044108|1.000000
\endverbatim

\par File
\ref WlzClassifyRCC.c "WlzClassifyRCC.c"
\par See Also
\ref WlzRCCClass "WlzRCCClass(3)"
\ref WlzRegConCalcRCC "WlzRegConCalcRCC(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		noAdj = 0,
		noEnc = 0,
		maxAdj = 8,
		ok = 1,
		usage = 0;
  double	*nrmVol = NULL;
  WlzRCCClass 	cls = WLZ_RCC_EMPTY;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr = NULL;
  WlzObject	*iObj[2] = {NULL};
  char  	*iObjFileStr[2] = {NULL};
  const char	*errMsg;
  static char	optList[] = "m:o:aeh",
		fileStrDef[] = "-";
  const double	eps = 1.0e-06;

  opterr = 0;
  outFileStr = fileStrDef;
  iObjFileStr[0] = fileStrDef;
  iObjFileStr[1] = fileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
        noAdj = 1;
	break;
      case 'e':
        noEnc = 1;
	break;
      case 'm':
        {
	  double t;

	  if((sscanf(optarg, "%lg", &t) != 1) || (t < -eps) || (t > INT_MAX))
	  {
	    usage = 1;
	  }
	  else
	  {
	    maxAdj = ALG_NINT(t);
	  }
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if((iObjFileStr[0] == NULL) || (*iObjFileStr[0] == '\0') ||
     (iObjFileStr[1] == NULL) || (*iObjFileStr[1] == '\0') ||
     (outFileStr      == NULL) || (*outFileStr      == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    int		i = 0;

    while((i < 2) && (optind < argc))
    {
      iObjFileStr[i] = *(argv + optind);
      ++optind;
      ++i;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
    ok = 0;
  }
  if(ok)
  {
    int		i = 0;

    while((errNum == WLZ_ERR_NONE) && (i <= 1))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((iObjFileStr[i] == NULL) ||
	  (*iObjFileStr[i] == '\0') ||
	  ((fP = (strcmp(iObjFileStr[i], "-")?
		  fopen(iObjFileStr[i], "r"): stdin)) == NULL) ||
	  ((iObj[i] = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to read object %d from file %s (%s)\n",
		       *argv, i, iObjFileStr[i], errMsg);
      }
      if(fP && strcmp(iObjFileStr[i], "-"))
      {
	(void )fclose(fP);
      }
      ++i;
    }
  }
  if(ok)
  {
    cls = WlzRegConCalcRCC(iObj[0], iObj[1], noEnc, noAdj, maxAdj,
    			   NULL, &nrmVol, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute classification (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?  fopen(outFileStr, "w"):
	      				 stdout)) != NULL))
    {
      int	i,
      		first = 1;
      unsigned int m;

      errNum = WLZ_ERR_NONE;
      for(i = 0; i <= WLZ_RCCIDX_CNT; ++i)
      {
	m = 1<<i;
        if(m & cls)
	{
	  const char *str;

	  str = WlzStringFromRCC((WlzRCCClass )m, NULL);
	  if(str == NULL)
	  {
	    errNum = WLZ_ERR_PARAM_TYPE;
	  }
	  else
	  {
	    if(fprintf(fP, "%s%s", (first)? "": "|", str) < 0)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    first = 0;
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
      if((errNum == WLZ_ERR_NONE) && first)
      {
        if(fprintf(fP, "WLZ_RCC_EMPTY") < 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        first = 1;
	for(i = 0; i <= WLZ_RCCIDX_CNT; ++i)
	{
	  m = 1<<i;
	  if(m & cls)
	  {
	    if(first)
	    {
	      first = 0;
	      if(fprintf(fP, "\t%g", nrmVol[i]) < 0)
	      {
		errNum = WLZ_ERR_WRITE_INCOMPLETE;
	      }
	    }
	    else
	    {
	      if(fprintf(fP, "|%g", nrmVol[i]) < 0)
	      {
	        errNum = WLZ_ERR_WRITE_INCOMPLETE;
	      }
	    }
	    if(errNum != WLZ_ERR_NONE)
	    {
	      break;
	    }
	  }
	}
	if((errNum == WLZ_ERR_NONE) && first)
	{
	  if(fprintf(fP, "0.0") < 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if(fprintf(fP, "\n") < 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  AlcFree(nrmVol);
  WlzFreeObj(iObj[0]);
  WlzFreeObj(iObj[1]);
  if(usage)
  {

    fprintf(stderr,
	    "Usage: %s"
	    "  [-a] [-e] [-h] [-m<max adj dist>] [-o<out file>]\n"
	    "\t\t[<in object 0>] [<in object 1>]\n"
	    "Version: %s\n"
	    "Options:\n"
	    "  -a        Don't include adjacency classifications.\n"
	    "  -e        Don't include enclosure classifications.\n"
	    "  -h        Help, prints this usage message.\n"
	    "  -m        Maximum adjacency distance in piels/voxels, this\n"
	    "            should be a resonable distance since domains will\n"
	    "            be dilated using this figure (set to %d)\n"
	    "  -o        Output file name.\n"
	    "Classifies the given domain objects using the region connected\n"
	    "calculus defined in Woolz:\n"
	    "  WLZ_RCC_EMPTY  - One or both of the domains is EMPTY.\n"
	    "  WLZ_RCC_DC     - DisConnected.\n"
	    "  WLZ_RCC_EC     - Edge Connected.\n"
	    "  WLZ_RCC_EQ     - EQual.\n"
	    "  WLZ_RCC_PO     - Partial Overlap.\n"
	    "  WLZ_RCC_TPP    - Tangential Proper Part.\n"
	    "  WLZ_RCC_TPPI   - Tangential Proper Part Inverse.\n"
	    "  WLZ_RCC_NTPP   - Non-Tangential Proper Part.\n"
	    "  WLZ_RCC_NTPPI  - Non-Tangential Proper Part Inverse.\n"
	    "  WLZ_RCC_TSUR   - Tangential SURrounds.\n"
	    "  WLZ_RCC_TSURI  - Tangential SURrounds Inverse.\n"
	    "  WLZ_RCC_NTSUR  - Non-Tangential SURrounds.\n"
	    "  WLZ_RCC_NTSURI - Non-Tangential SURrounds Inverse.\n"
	    "  WLZ_RCC_ENC    - ENCloses.\n"
	    "  WLZ_RCC_ENCI   - ENCloses Inverse.\n"
	    "  WLZ_RCC_ADJ    - ADJacent.\n"
	    "  WLZ_RCC_ADJI   - ADJacent Inverse.\n"
	    "See WlzRegConCalcRCC(3) for more information.\n"
	    "Example:\n%s somites.wlz embryonic.wlz\n"
	    "giving output\n"
	    "WLZ_RCC_TPP|WLZ_RCC_ENC        0.044108|1.000000\n"
	    "which should be interpreted as the classification that the\n"
	    "domain somites is a tangential proper part of the domain\n"
	    "embryonic and the domain somites is enclosed by the domain\n"
	    "embryonic. In this case the scores are the intersection of the\n"
	    "two domains divided by their union, followed by the proportion\n"
	    "of somites which is within the convex hull of embryonic.\n",
	    *argv,
	    WlzVersion(),
	    maxAdj,
	    *argv);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
