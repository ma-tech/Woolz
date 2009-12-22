#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshTransformObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzCMeshTransformObj.c
* \author       Zsolt Husz
* \date         December 2009
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2009 Medical research Council, UK.
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
* \todo         -
* \bug          None known.
*
* \ref 		WlzCMeshTransformObj "WlzCMeshTransformObj"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcopyobj WlzCMeshTransformObj
\par Name
WlzCMeshTransformObj - transforms an object using a constrained mesh transform.
\par Synopsis
\verbatim
WlzCMeshTransformObj [-h] [-t<input transform>] [-i] [-N] [-x<interpolation value> [-o<output woolz file>]] |
  [-s<number of interpolations> [-b <output body>] [-e <output extension>]] ] [<input object>]
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
  int		dim = 0,
  		inv = 0,
  		ok = 1,
  		option,
  		usage = 0;
  FILE		*inFP = NULL,
  		*outFP = NULL;
  char		*txFileStr = NULL,
                *inFileStr = NULL,
                *outFileStr = NULL,
                *outObjFileBodyStr = NULL,
                *outObjFileExtStr = NULL,
                outObjFileStr[FILENAME_MAX];

  const char	*errMsgStr;
  WlzObject	*tx = NULL,
                *inObj = NULL,
                *outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "iho:t:s:x:nb:e:";
  const char    txFileStrDef[] = "-",
  		inFileStrDef[] = "-",
                outFileStrDef[] = "-",
                outObjFileExtStrDef[] = "wlz";

  double        transition = 1.0;
  int           stepsno  = 1;
  int           useStep = 0;
  WlzInterpolationType interp = WLZ_INTERPOLATION_LINEAR;

  opterr = 0;
  txFileStr = (char *)txFileStrDef;
  inFileStr = (char *)inFileStrDef;
  outFileStr = (char *)outFileStrDef;
  outObjFileExtStr = (char *)outObjFileExtStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'i':
        inv = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 't':
        txFileStr = optarg;
        useStep = 0;
        break;
      case 'x':
        transition = strtod(optarg, NULL);
        break;
      case 'n':
        interp = WLZ_INTERPOLATION_NEAREST;
        break;
      case 's':
        stepsno = atoi(optarg);
        useStep = 1;
        break;
      case 'b':
        outObjFileBodyStr = optarg;
        break;
      case 'e':
        outObjFileExtStr = optarg;
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
  if(useStep && ok && (transition<0 || transition>1)) {
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
                   "%s: output filename length exceeded!\n",
                   *argv);
  }

  if(ok)
  {
    if(((inFP = (strcmp(txFileStr, "-")?
	        fopen(txFileStr, "r"): stdin)) == NULL) ||
       ((tx = WlzAssignObject(WlzReadObj(inFP, &errNum), NULL)) == NULL))
    {
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(tx->type)
      {
        case WLZ_CMESH_2D:
	  dim = 2;
	  break;
	case WLZ_CMESH_3D:
	  dim = 3;
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
  if(ok)
  {
      errNum = WLZ_ERR_READ_EOF;
      if(((inFP = (strcmp(inFileStr, "-")?
                  fopen(inFileStr, "r"): stdin)) == NULL) ||
         ((inObj = WlzAssignObject(WlzReadObj(inFP, &errNum), NULL)) == NULL))
      {
        if(errNum == WLZ_ERR_NONE)
        {
          errNum = WLZ_ERR_READ_INCOMPLETE;
        }
      }

      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
        (void )WlzStringFromErrorNum(errNum, &errMsgStr);
        (void )fprintf(stderr,
            "%s: Failed to open input object file %s (%s).\n",
            *argv, inFileStr, errMsgStr);
      }
      if((inFP != NULL) && strcmp(inFileStr, "-"))
      {
        (void )fclose(inFP); inFP = NULL;
      }
  }
  if(ok && inv)
  {
    if((errNum = WlzCMeshTransformInvert(tx)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to invert transform (%s).\n",
                     argv[0],
                     errMsgStr);
    }
  }


  if(ok) {
      if (useStep) {
        int step = FILENAME_MAX;
        WlzObject *trobj = NULL;
        step = 0;
        for  (step = 0; step < stepsno && ok; step++) {

           (void )sprintf(outObjFileStr,
                    "%s%06d.%s", outObjFileBodyStr, step, outObjFileExtStr);

           trobj = WlzCopyScaleCMeshValue((double)step/(double)(stepsno-1), tx, &errNum);

           if (errNum != WLZ_ERR_NONE)
               break;

           trobj = WlzAssignObject(trobj, &errNum);
           if (errNum != WLZ_ERR_NONE)
               break;

           if(errNum == WLZ_ERR_NONE)
              outObj = WlzCMeshTransformObj( inObj, trobj, interp, &errNum);  //volume transform

           if(errNum != WLZ_ERR_NONE)
           {
             ok = 0;
             (void )WlzStringFromErrorNum(errNum, &errMsgStr);
             (void )fprintf(stderr,
                 "%s: Failed to apply transform (%s).\n",
                 *argv, errMsgStr);
             errNum = WLZ_ERR_NONE; // error already handled
           }
           if(outObjFileBodyStr) {
               if((outFP = fopen(outObjFileStr, "w")) == NULL)
               {
                 ok = 0;
                 (void )fprintf(stderr,
                                "%s: failed to open output file %s.\n",
                                *argv, outObjFileStr);
               }
           }
           else
           {
               outFP = stdout;
           }
           if((errNum = WlzWriteObj(outFP, outObj)) != WLZ_ERR_NONE)
           {
             ok = 0;
             (void )WlzStringFromErrorNum(errNum, &errMsgStr);
             (void )fprintf(stderr,
                            "%s: failed to write output object (%s).\n",
                            *argv, errMsgStr);
           }
           WlzFreeObj(trobj);
           WlzFreeObj(outObj);
       }
       if(errNum != WLZ_ERR_NONE && ok)
       {
           ok = 0;
           (void )WlzStringFromErrorNum(errNum, &errMsgStr);
           (void )fprintf(stderr,
               "%s: Failed to create object sequence (%s).\n",
               *argv, errMsgStr);
       }
      } else {
        WlzObject *trobj = tx;
        errNum = WlzScaleCMeshValue (transition, trobj);
        if(errNum == WLZ_ERR_NONE)
        {
            outObj = WlzCMeshTransformObj( inObj, trobj, interp, &errNum);  //volume transform
            if(errNum != WLZ_ERR_NONE)
            {
              ok = 0;
              (void )WlzStringFromErrorNum(errNum, &errMsgStr);
              (void )fprintf(stderr,
                  "%s: Failed to apply transform (%s).\n",
                  *argv, errMsgStr);
            }
         }
          if(ok)
          {
            errNum = WLZ_ERR_WRITE_EOF;
            if(((outFP  = (strcmp(outFileStr, "-")?
                      fopen(outFileStr, "w"):
                      stdout)) == NULL) ||
               ((errNum = WlzWriteObj(outFP , outObj)) != WLZ_ERR_NONE))
            {
              ok = 0;
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
      }
  }
  (void )WlzFreeObj(tx);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
      fprintf(stderr,
            "Usage: %s [-h] [-t<input transfrom>] [-i] [-N] \n"
            "  [-x<interpolation value> [-o<output woolz file>]] | \n"
            "  [-s<number of interpolations> [-b <output body>] \n"
            "    [-e <output extension>]] ] [<input object>]\n"
            "Reads a constrained mesh transform object, applies it\n"
            "to an object and generates the outputs transformed object.\n"
            "Partly transformed object can be generated if\n"
            "interpolation value between 0 and 1 (0 means no \n"
            "transformation). Alternatively a set of interpolations\n"
            "covering the full range of transformation is computed\n"
            "if the number of interpolations are given. For this,.\n"
            "the output base filename and its extension are\n"
            "separately given.\n"
            "Options are:\n"
            "  -h  Help, prints this usage message.\n"
            "  -t  Transform object.\n"
            "  -i  Invert the transform after reading.\n"
            "  -n  Use WLZ_INTERPOLATION_NEAREST (default WLZ_INTERPOLATION_LINEAR).\n"
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
            argv[0], argv[0],argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
