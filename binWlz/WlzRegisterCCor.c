#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzRegisterCCor_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzRegisterCCor.c
* \author       Bill Hill
* \date         August 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
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
* \brief	Registers a pair of 2D domain objects with grey values
* 		using frequency domain cross-correlation.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzregisterccor "WlzRegisterCCor"
*/

/*!
\ingroup BinWlz
\defgroup wlzregisterccor WlzRegisterCCor
\par Name
WlzRegisterCCor - registers a pair of 2D domain objects with grey values
                  using frequency domain cross-correlation.
\par Synopsis
\verbatim
WlzRegisterCCor [-h] [-v] [-o<out obj>] [-i <init tr>]
                [-t] [-r] [-R #] [-T #,#]
		[<in obj 0>] [<in obj 1>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name for affine transform.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Initial affine transform object.</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Find the rigid body (aka registration) transform, default.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Find the translation only transform.</td>
  </tr>
  <tr> 
    <td><b>-R</b></td>
    <td>Maximum rotation.</td>
  </tr>
  <tr> 
    <td><b>-T</b></td>
    <td>Maximum translation.</td>
  </tr>
</table>
\par Description
Attempts to register two objects using an frequency domain
cross-correlation algorithm.  The two objects must be 2D spatial
domain objects with grey values.
It is important to give sensible maximum translation and rotation
values. If they are to small the itteration probably will not converge
to the true maximum cross-correlation values, but if the values are
too large the algorithm will take both more time and memory.
The input objects are read from stdin and values are written to stdout
unless the filenames are given.
\par Examples
\verbatim
WlzRegisterCCor -o out-tr.wlz -R 10 -T 20,30 -r in0.wlz in1.wlz
\endverbatim
A rigid body affine transform is found by registering in1.wlz to
in0.wlz using a cross-correlation based registration algorithm. The
affine transform is then written to out-tr.wlz. The maximum rotation
in a single itteration is 10 degrees and the maximum translation in
a single itteration is 20 columns and 30 rows.
\par File
\ref WlzRegisterCCor.c "WlzRegisterCCor.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzregistericp "WlzRegisterICP(1)"
\ref WlzRegCCorObjs "WlzRegCCorObjs(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
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
  int		idx,
		option,
		cConv = 0,
		verbose = 0,
		ok = 1,
		usage = 0;
  double 	maxRot,
  		cCor;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzTransformType trType = WLZ_TRANSFORM_2D_REG;
  WlzDomain	outDom;
  WlzValues	nullVal;
  WlzObject	*inTrObj = NULL,
  		*outObj = NULL;
  WlzObject	*inObj[2];
  WlzDVertex2 	maxTran;
  FILE		*fP = NULL;
  char 		*inTrObjFileStr = NULL,
  		*outObjFileStr;
  char  	*inObjFileStr[2];
  const char	*errMsg;
  static char	optList[] = "i:o:R:T:hrtv",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outDom.core = NULL;
  nullVal.core = NULL;
  inObj[0] = NULL;
  inObj[1] = NULL;
  maxRot = 45;
  maxTran.vtX = maxTran.vtY = 100.0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr[0] = inObjFileStrDef;
  inObjFileStr[1] = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'i':
        inTrObjFileStr = optarg;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'r':
        trType = WLZ_TRANSFORM_2D_REG;
	break;
      case 't':
        trType = WLZ_TRANSFORM_2D_TRANS;
	break;
      case 'R':
	if(sscanf(optarg, "%lg", &maxRot) != 1)
	{
	  usage = 1;
	}
        break;
      case 'T':
	if(sscanf(optarg, "%lg,%lg", &(maxTran.vtX), (&maxTran.vtY)) != 2)
	{
	  usage = 1;
	}
        break;
      case 'v':
        verbose = 1;
	break;
      case 'h':
      default:
        usage = 1;
	break;
    }
  }
  ok = !usage;
  if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
     (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    idx = 0;
    while((idx < 2) && (optind < argc))
    {
      inObjFileStr[idx] = *(argv + optind);
      ++optind;
      ++idx;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
    ok = 0;
  }
  if(ok && inTrObjFileStr)
  {
    /* Read initial affine transform. */
    if(((fP = (strcmp(inTrObjFileStr, "-")?
               fopen(inTrObjFileStr, "r"): stdin)) == NULL) ||
       ((inTrObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	     "%s: Failed to read initial affine transform from file %s (%s)\n",
	     *argv, errMsg);
    }
    if(fP && strcmp(inTrObjFileStr, "-"))
    {
      fclose(fP);
    }
    if(inTrObj &&
       ((inTrObj->type != WLZ_AFFINE_TRANS) || (inTrObj->domain.core == NULL)))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Initial affine transform object invalid type\n",
		     *argv);
    }
  }
  if(ok)
  {
    /* Read objects. */
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx] = WlzAssignObject(WlzReadObj(fP,
	  					    &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to read object %d from file %s (%s)\n",
		       *argv, idx, inObjFileStr[idx], errMsg);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
      }
      ++idx;
    }
  }
  if(ok)
  {
    /* Check object types. */
    if((inObj[0]->type != WLZ_2D_DOMAINOBJ) ||
       (inObj[0]->type != inObj[1]->type))
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if((inObj[0]->domain.core == NULL) ||
            (inObj[1]->domain.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else if((inObj[0]->values.core == NULL) ||
            (inObj[1]->values.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: Input object(s) not appropriate\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    outDom.t = WlzRegCCorObjs(inObj[0], inObj[1],
			      inTrObj? inTrObj->domain.t: NULL, trType,
			      maxTran, maxRot, 1000,
			      &cConv, &cCor, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to register objects (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    outObj = WlzMakeMain(WLZ_AFFINE_TRANS, outDom, nullVal, NULL, NULL,
    			 &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to make affine transform object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?  fopen(outObjFileStr, "w"):
	      				    stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok && verbose)
  {
    (void )printf("%d %g\n", cConv, cCor);
  }
  if(inObj[0])
  {
    (void )WlzFreeObj(inObj[0]);
  }
  if(inObj[1])
  {
    (void )WlzFreeObj(inObj[1]);
  }
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  else if(outDom.core)
  {
    (void )WlzFreeAffineTransform(outDom.t);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-v] [-o<out obj>] [-i <init tr>]\n"
    "                      [-t] [-r] [-R #] [-T #,#]\n"
    "                      [<in obj 0>] [<in obj 1>]\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -v  Verbose output to the standard output. This consists of two\n"
    "      numbers: The convergence flag and the cross correlation value.\n"
    "      When the algorithm has converged the convergence flag has\n"
    "      non-zero value. Cross correlation values are normalized to the\n"
    "      range [0.0-1.0].\n"
    "  -o  Output file name for affine transform.\n"
    "  -i  Initial affine transform object.\n"
    "  -r  Find the rigid body (aka registration) transform, default.\n"
    "  -t  Find the translation only transform.\n"
    "  -R  Maximum rotation.\n"
    "  -T  Maximum translation.\n"
    "Attempts to register two objects using an frequency domain\n"
    "cross-correlation algorithm.  The two objects must be 2D spatial\n"
    "domain objects with grey values.\n"
    "It is important to give sensible maximum translation and rotation\n"
    "values. If they are to small the itteration probably will not converge\n"
    "to the true maximum cross-correlation values, but if the values are\n"
    "too large the algorithm will take both more time and memory.\n"
    "The input objects are read from stdin and values are written to stdout\n"
    "unless the filenames are given.\n",
    *argv,
    " -o out-tr.wlz -R 10 -T 20,30 -r in0.wlz in1.wlz\n"
    "A rigid body affine transform is found by registering in1.wlz to\n"
    "in0.wlz using a cross-correlation based registration algorithm. The\n"
    "affine transform is then written to out-tr.wlz. The maximum rotation\n"
    "in a single itteration is 10 degrees and the maximum translation in\n"
    "a single itteration is 20 columns and 30 rows.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
