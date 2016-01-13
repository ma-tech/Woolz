#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshTransformObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshTransformObj.c
* \author       Zsolt Husz, Bill Hill
* \date         December 2009
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
* \brief	Transforms an object using a constrained mesh transform
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshtransformobj "WlzCMeshTransformObj"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcmeshtransformobj WlzCMeshTransformObj
\par Name
WlzCMeshTransformObj - transforms an object using a constrained mesh transform.
\par Synopsis
\verbatim
WlzCMeshTransformObj [-h] [-t<input transform>] [-i] [-N]
                     [-x<interpolation value>] [-o<output woolz file>]
                     [-s<number of interpolations> [-b <output body>]
		     [-E] [-e <output extension>]] ] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Transform object.</td>
  </tr>
  <tr> 
    <td><b>-E</b></td>
    <td>Output evaluation times to stderr.</td>
  </tr>
  <tr>
    <td><b>-i</b></td>
    <td>Invert the transform after reading.</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>Use WLZ_INTERPOLATION_NEAREST (default WLZ_INTERPOLATION_LINEAR).</td>
  </tr>
  <tr>
    <td><b>-x</b></td>
    <td>Interpolation value, with 0 <= ivalue <= 1.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object filename.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Number of intermediate interpolations.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file body.</td>
  </tr>

  <tr> 
    <td><b>-e</b></td>
    <td>Output object file extension (default wlz).</td>
  </tr>
</table>
\par Description
Reads a constrained mesh transform object, applies it to an object 
and generates the outputs transformed object. Partly transformed 
object can be generated if interpolation value between 0 and 1 
(0 means no transformation). Alternatively a set of interpolations
covering the full range of transformation is computed if the number
of interpolations are given. For this, the output base filename and
its extension are separately given.
\par Example
\verbatim
WlzCMeshTransformObj -t transform.wlz -i -o out.wlz in.wlz
\endverbatim
Applies the inverse transform read from transform.wlz to the in.wlz object
and writes the result to out.wlz
\par Example
\verbatim
WlzCMeshTransformObj -t transform.wlz -s 5 -b out_ -e wlz in.wlz
\endverbatim
Applies the transform read from transform.wlz to the in.wlz and
generates out_000000.wlz ... out_000005.wlz with interpolation
values 0, 0.2, ... 1.0.

\par File
\ref WlzCMeshTransformObj.c "WlzCMeshTransformObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */


int		main(int argc, char *argv[])
{
  int		objCnt,
  		option,
  		inv = 0,
		nStep = 1,
		timer = 0,
		useStep = 0,
  		ok = 1,
  		usage = 0;
  double        transition = 1.0;
  FILE		*inFP = NULL,
  		*outFP = NULL;
  char		*txFileStr = NULL,
                *inFileStr = NULL,
                *outFileStr = NULL,
                *outObjFileBodyStr = NULL,
                *outObjFileExtStr = NULL,
                outObjFileStr[FILENAME_MAX];
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  const char	*errMsgStr;
  WlzObject	*trObj = NULL,
                *inObj = NULL,
                *outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  struct timeval times[3];
  static char   optList[] = "iho:t:s:x:ELb:e:";
  const char    txFileStrDef[] = "-",
  		inFileStrDef[] = "-",
                outFileStrDef[] = "-",
                outObjFileExtStrDef[] = "wlz";


  opterr = 0;
  txFileStr = (char *)txFileStrDef;
  inFileStr = (char *)inFileStrDef;
  outFileStr = (char *)outFileStrDef;
  outObjFileExtStr = (char *)outObjFileExtStrDef;
  (void )memset(times, 0, 3 * sizeof(struct timeval));
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'b':
        outObjFileBodyStr = optarg;
        break;
      case 'e':
        outObjFileExtStr = optarg;
        break;
      case 'i':
        inv = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 't':
        txFileStr = optarg;
        break;
      case 'E':
        timer = 1;
	break;
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
        break;
      case 's':
        useStep = 1;
	if(sscanf(optarg, "%d", &nStep) != 1)
	{
	  usage = 1;
	}
        break;
      case 'x':
        useStep = 0;
        if(sscanf(optarg, "%lg", &transition) != 1)
	{
	  usage = 1;
	}
        break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      inFileStr = argv[optind];
    }
  }
  ok = usage == 0;
  if(useStep && ok && (transition < 0 || transition > 1.0))
  {
      ok = 0;
      (void )fprintf(stderr,
          "%s: Error: transition must be between 0 and 1.\n",
          *argv);
  }
  if(outObjFileBodyStr && outObjFileExtStr &&
     (strlen(outObjFileExtStr) +
      strlen(outObjFileExtStr)) > (FILENAME_MAX - 24))
  {
    ok = 0;
    (void )fprintf(stderr,
                   "%s: Output filename length exceeded!\n",
                   *argv);
  }

  if(ok)
  {
    if(((inFP = (strcmp(txFileStr, "-")?
	        fopen(txFileStr, "r"): stdin)) == NULL) ||
       ((trObj = WlzAssignObject(WlzReadObj(inFP, &errNum), NULL)) == NULL))
    {
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(trObj->type)
      {
        case WLZ_CMESH_2D:  /* FALLTHROUGH */
        case WLZ_CMESH_2D5: /* FALLTHROUGH */
	case WLZ_CMESH_3D:
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	  "%s: Failed to read constrained mesh transform object from\n"
          "file %s (%s)\n",
	  *argv, txFileStr, errMsgStr);
    }
    if((inFP != NULL) && strcmp(txFileStr, "-"))
    {
      (void )fclose(inFP); inFP = NULL;
    }
  }
  if(ok && inv)
  {
    WlzObject *nTrObj;

    nTrObj = WlzCMeshTransformInvert(trObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzFreeObj(trObj);
      trObj = WlzAssignObject(nTrObj, NULL);
    }
    else
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to invert transform (%s).\n",
                     argv[0], errMsgStr);
    }
  }
  inFP = NULL;
  if(ok)
  {
      if((inFP = (strcmp(inFileStr, "-")?
                 fopen(inFileStr, "r"): stdin)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
  }
  if(ok)
  {
    objCnt = 0;
    while(errNum == WLZ_ERR_NONE)
    {
      errNum = WLZ_ERR_READ_EOF;
      if((inObj = WlzAssignObject(WlzReadObj(inFP, &errNum), NULL)) == NULL)
      {
	if(objCnt > 0)
	{
	  errNum = WLZ_ERR_EOO;
	}
	else
	{
	  (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	  (void )fprintf(stderr,
	      "%s: Failed to open input object file %s (%s).\n",
	      *argv, inFileStr, errMsgStr);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(useStep)
	{
	  int step = FILENAME_MAX;

	  step = 0;
	  while((step < nStep) && (errNum == WLZ_ERR_NONE))
	  {
	    WlzObject *stepTrObj = NULL;

	    (void )sprintf(outObjFileStr, "%s%04d%04d.%s",
			   outObjFileBodyStr, objCnt, step, outObjFileExtStr);
	    stepTrObj = WlzAssignObject(
			WlzCopyScaleCMeshValue(
			  (double )step / (double )(nStep -1),
			  trObj, &errNum), NULL);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      outObj = WlzCMeshTransformObj(inObj, stepTrObj, interp, &errNum);
	      if(errNum != WLZ_ERR_NONE)
	      {
		(void )WlzStringFromErrorNum(errNum, &errMsgStr);
		(void )fprintf(stderr,
		    "%s: Failed to apply transform (%s).\n",
		    *argv, errMsgStr);
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(outObjFileBodyStr)
	      {
		if((outFP = fopen(outObjFileStr, "w")) == NULL)
		{
		  (void )fprintf(stderr,
		      "%s: Failed to open output file %s.\n",
		      *argv, outObjFileStr);
		}
	      }
	      else
	      {
		outFP = stdout;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if((errNum = WlzWriteObj(outFP, outObj)) != WLZ_ERR_NONE)
	      {
		(void )WlzStringFromErrorNum(errNum, &errMsgStr);
		(void )fprintf(stderr,
		    "%s: Failed to write output object (%s).\n",
		    *argv, errMsgStr);
	      }
	    }
	    WlzFreeObj(stepTrObj);
	    WlzFreeObj(outObj);
	    ++step;
	  }
	  if(errNum != WLZ_ERR_NONE && ok)
	  {
	    (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	    (void )fprintf(stderr,
		"%s: Failed to create object sequence (%s).\n",
		*argv, errMsgStr);
	  }
	}
	else
	{
	  errNum = WlzScaleCMeshValue(transition, trObj);
	  if(errNum == WLZ_ERR_NONE)
	  {
            gettimeofday(times + 0, NULL);
	    outObj = WlzCMeshTransformObj(inObj, trObj, interp, &errNum);
            gettimeofday(times + 1, NULL);
	    if(errNum != WLZ_ERR_NONE)
	    {
	      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	      (void )fprintf(stderr,
		  "%s: Failed to apply transform (%s).\n",
		  *argv, errMsgStr);
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WLZ_ERR_WRITE_EOF;
	    if(((outFP  = (strcmp(outFileStr, "-")?
		      fopen(outFileStr, "w"):
		      stdout)) == NULL) ||
		((errNum = WlzWriteObj(outFP , outObj)) != WLZ_ERR_NONE))
	    {
	      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	      (void )fprintf(stderr,
		  "%s: Failed to write output object %s (%s).\n",
		  *argv, outFileStr, errMsgStr);
	    }
	    if(outFP && strcmp(outFileStr, "-"))
	    {
	      (void )fclose(outFP);
	    }
	  }
	  (void )WlzFreeObj(outObj);
	}
      }
      (void )WlzFreeObj(inObj);
      ++objCnt;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(errNum == WLZ_ERR_EOO)
      {
	errNum = WLZ_ERR_NONE;
      }
      else
      {
	ok = 0;
      }
    }
  }
  (void )WlzFreeObj(trObj);
  if((inFP != NULL) && strcmp(inFileStr, "-"))
  {
    (void )fclose(inFP); inFP = NULL;
  }
  if(timer)
  {
    ALC_TIMERSUB(times + 1, times + 0, times + 2);
    (void )fprintf(stderr,
    "%s: Elapsed time to transform the object was: %g seconds\n",
    *argv,
    times[2].tv_sec + (0.000001 * times[2].tv_usec));
  }
  if(usage)
  {
      fprintf(stderr,
            "Usage: %s [-h] [-t<input transfrom>] [-i] [-N] \n"
            "        [-x<interpolation value> [-o<output woolz file>]] | \n"
            "        [-s<number of interpolations> [-b <output body>] \n"
            "        [-e <output extension>]] ] [<input object>]\n"
            "Reads a constrained mesh transform object, applies it\n"
            "to objects generating output transformed objects.\n"
            "Partly transformed object can be generated if interpolation\n"
            "values between 0 and 1 are given(0 means no transformation)\n"
            "Alternatively a set of interpolations covering the full range\n"
	    "of transformation is computed if the number of interpolations\n"
	    "is given. For this, the output base filename and its extension\n"
	    "are must be supplied.\n"
	    "Version: %s\n"
            "Options are:\n"
            "  -h  Help, prints this usage message.\n"
            "  -t  Transform object.\n"
            "  -i  Invert the transform after reading.\n"
            "  -E  Output execution time to stderr.\n"
            "  -L  Use WLZ_INTERPOLATION_LINEAR (instead of the default\n"
	    "      WLZ_INTERPOLATION_NEAREST).\n"
            "  -x  Interpolation value, with 0 <= ivalue <= 1.\n"
            "  -o  Output object file\n"
            "  -s  Number of intermediate interpolations.\n"
            "  -b  Output object file body\n"
            "  -e  Output object file extension (default wlz)\n"
            "Parameters -x (and the auxiliary -o) respectively\n"
            "-s (and the auxiliary -b and -e) are mutually exclusive.\n"
	    "Example:\n"
            "  %s -t transform.wlz -i -o out.wlz in.wlz\n"
            "Applies the inverse transform read from transform.wlz to\n"
            "the in.wlz object and writes the result to out.wlz\n"
            "  %s -t transform.wlz -s 5 -b out_ -e wlz in.wlz\n"
            "Applies the transform read from transform.wlz to\n"
            "the in.wlz and generates out_000000.wlz ... out_000005.wlz\n"
            "with interpolation values 0, 0.2, ... 1.0\n.",
            argv[0],
	    WlzVersion(),
	    argv[0],
	    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
