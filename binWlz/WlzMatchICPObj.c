#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzMatchICPObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzMatchICPObj.c
* \author       Bill Hill
* \date         March 2002
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
* \brief	Object matching using ICP based registration to
* 		find corresponding points.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzmatchicpobj "WlzMatchICPObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzmatchicpobj WlzMatchICPObj
\par Name
WlzMatchICPObj - object matching using ICP based registration.
\par Synopsis
\verbatim
WlzMatchICPObj [-b#] [-d#] [-D#] [-E#] [-o#] [-t#] [-g#] [-i#] [-s#]
               [-A#] [-S#] [-F#] [-I#] [-N] [-h]
	       [<input object 0>] [<input object 1>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Shell breaking control value.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Debug output file.</td>
  </tr>
  <tr> 
    <td><b>-D</b></td>
    <td>Debug flag, 0 no debug output or 1 output untransformed
        decomposed source model.</td>
  </tr>
  <tr> 
    <td><b>-E</b></td>
    <td>Tolerance in mean registration metric value.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Initial affine transform.</td>
  </tr>
  <tr> 
    <td><b>-P</b></td>
    <td>Minimum number of simplices per matched shell segment, with a
        pair of correspondence points possibly being generated per
	matched shell segment.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Maximum number of iterations.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Minmum number of simplicies.</td>
  </tr>
  <tr> 
    <td><b>-A</b></td>
    <td>Maximum angle (degrees) from a global transformation.</td>
  </tr>
  <tr> 
    <td><b>-S</b></td>
    <td>Maximum displacement from a global transformed position.</td>
  </tr>
  <tr> 
    <td><b>-F</b></td>
    <td>Maximum deformation from a global transformation.</td>
  </tr>
  <tr> 
    <td><b>-I</b></td>
    <td>Implausibility threshold for rejecting implausible
        correspondence points which should be greater than zero,
although the useful range is probably [0.5-2.5]. Higher
values allow more implausible matches to be returned.</td>
  </tr>
  <tr> 
    <td><b>-N</b></td>
    <td>Number of match points in neighbourhood when checking the
        plausibility of the correspondence points.</td>
  </tr>
</table>
\par Description
Reads a pair of contours and computes a set of correspondence points
using an ICP based matching algorithm. The correspondence points are
ranked by plausibility, with the most plausible first.
The tie points are generated using and algorithm which
continualy matches fragments of the contours (shells),
while these fragments ar decomposed into smaller fragments.
\par Examples
\verbatim
WlzMatchICPObj -I1000 -N2 -A15 -S15.0 -F0.5 -s15 -b1000
               -o out.num pair_A_OPT_c_s50.wlz pair_A_HE_c_s50.wlz
\endverbatim
Reads two contour models from files pair_A_OPT_c_s50.wlz and
pair_A_HE_c_s50.wlz, then generates a list of tie points
and writes them to the file out.num.
The fragments of the contour models (shells)
will decompose until a decomposition would
result in a shell having less than 15 edges.
If any source shell fails to match a target shell using
constrains of 15 degrees maximum rotation,
15 pixels displacement and a maximum deformation of 50%
then the shell will be rejected.
At the end of the decomposition process
each shell will be used to generate
two tie points.
The implausibility threshold is set so high that no ti points will be
rejected as being implausible.
\par File
\ref WlzMatchICPObj.c "WlzMatchICPObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzbasisfntransformobj "WlzBasisFnTransformObj(1)"
\ref WlzMatchICPObjs "WlzMatchICPObjs(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>
#include <string.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           idx,
		brkFlg = INT_MAX,
		dbgFlg = 0,
		nMatch,
		maxItr = 200,
		minSpx = 15,
		minSegSpx = 10,
		matchImpNN = 7,
  		option,
  		ok = 1,
		usage = 0;
  double	delta = 0.01,
  		maxAng = 30.0 * ALG_M_PI / 180.0,
  		maxDeform = 0.5,
  		maxDisp = 25.0,
		matchImpThr = 1.5;
  FILE		*fP = NULL;
  char		*inTrFileStr,
  		*outFileStr,
		*outDbgFileStr;
  char		*inObjFileStr[2];
  WlzObject	*inTrObj;
  WlzObject	*inObj[2];
  WlzVertexP	matchTP,
  		matchSP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "ahrA:b:d:D:E:F:N:o:t:i:I:s:P:S:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  inTrObj = NULL;
  inObj[0] = inObj[1] = NULL;
  outFileStr = (char *)outFileStrDef;
  outDbgFileStr = NULL;
  inTrFileStr = NULL;
  inObjFileStr[0] = inObjFileStr[1] = (char *)inObjFileStrDef;
  matchTP.v = matchSP.v = NULL;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        if((sscanf(optarg, "%d", &brkFlg) != 1) || (brkFlg < 0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'd':
	dbgFlg = 1;
        outDbgFileStr = optarg;
	break;
      case 'A':
        if((sscanf(optarg, "%lg", &maxAng) != 1) ||
	   (maxAng < 0.0) || (maxAng > 180.0))
	{
	  usage = 1;
	  ok = 0;
	}
	else
	{
	  maxAng *= ALG_M_PI / 180.0;
	}
	break;
      case 'D':
        if((sscanf(optarg, "%d", &dbgFlg) != 1) ||
	   (dbgFlg < 0) || (dbgFlg > 2))
        {
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'E':
        if((sscanf(optarg, "%lg", &delta) != 1) || (delta < 0.0))
        {
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'F':
        if((sscanf(optarg, "%lg", &maxDeform) != 1) || (maxDeform < 0.0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'N':
        if((sscanf(optarg, "%d", &matchImpNN) != 1) || (matchImpNN < 1))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 't':
        inTrFileStr = optarg;
	break;
      case 'i':
        if((sscanf(optarg, "%d", &maxItr) != 1) || (maxItr <= 0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'I':
        if((sscanf(optarg, "%lg", &matchImpThr) != 1) || (matchImpThr < 0.0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'P':
        if(sscanf(optarg, "%d", &minSegSpx) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
        if((sscanf(optarg, "%d", &minSpx) != 1) || (minSpx <= 10))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'S':
        if((sscanf(optarg, "%lg", &maxDisp) != 1) || (maxDisp <= 0.0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 2) != argc)
      {
        usage = 1;
	ok = 0;
      }
      else
      {
        inObjFileStr[0] = *(argv + optind);
        inObjFileStr[1] = *(argv + optind + 1);
      }
    }
  }
  if(ok)
  {
    idx = 0;
    while(ok && (idx < 2))
    {
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx] = WlzAssignObject(WlzReadObj(fP, &errNum),
	  				 NULL)) == NULL) ||
	  (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to read object from file %s\n",
		       *argv, inObjFileStr[idx]);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
        fP = NULL;
      }
      ++idx;
    }
  }
  if(ok && inTrFileStr && (*inTrFileStr != '\0'))
  {
    if(((fP = (strcmp(inTrFileStr, "-")?
    	       fopen(inTrFileStr, "r"): stdin)) == NULL) ||
       ((inTrObj = WlzAssignObject(WlzReadObj(fP, &errNum),
       				   NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s Failed to read the initial affine transform from\n"
      		     "file %s.\n",
		     *argv, inTrFileStr);
    }
    if(fP && strcmp(inTrFileStr, "-"))
    {
      fclose(fP);
      fP = NULL;
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
    errNum = WlzMatchICPObjs(inObj[0], inObj[1],
    			     (inTrObj)? inTrObj->domain.t: NULL,
			     &nMatch, &matchTP, &matchSP,
			     maxItr, minSpx, minSegSpx, brkFlg,
			     maxDisp, maxAng, maxDeform,
			     matchImpNN, matchImpThr, delta);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to match contours.\n",
      		     argv[0]);
    }
  }
  if(ok)
  {
    if(outDbgFileStr && (dbgFlg > 0))
    {
      errNum = WLZ_ERR_WRITE_EOF;
      if(((fP = (strcmp(outDbgFileStr, "-")?
      	        fopen(outDbgFileStr, "w"):
		stdout)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open debug object "
		       "file %s (%s).\n",
		       *argv, outDbgFileStr, errMsg);
      }
      else if((errNum = WlzWriteObj(fP, inObj[1])) != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write debug object "
		       "to file %s (%s).\n",
		       *argv, outDbgFileStr, errMsg);
      }
      if(fP && strcmp(outDbgFileStr, "-"))
      {
        (void )fclose(fP);
      }
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"):
	      stdout)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to open correspondences "
		     "file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    else
    {
      /* Output the correspondces: Assumed to be WlzDVertex2. */
      for(idx = 0; idx <nMatch; ++idx)
      {
        fprintf(fP, "%g %g %g %g\n",
		(matchSP.d2 + idx)->vtX, (matchSP.d2 + idx)->vtY,
		(matchTP.d2 + idx)->vtX - (matchSP.d2 + idx)->vtX,
		(matchTP.d2 + idx)->vtY - (matchSP.d2 + idx)->vtY);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  AlcFree(matchTP.v);
  AlcFree(matchSP.v);
  (void )WlzFreeObj(inTrObj);
  (void )WlzFreeObj(inObj[0]);
  (void )WlzFreeObj(inObj[1]);
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s",
      *argv,
      " [-b#] [-d#] [-D#] [-E#] [-o#] [-t#] [-g#] [-i#] [-s#]\n"
      "          [-A#] [-S#] [-F#] [-I#] [-N]\n"
      "          [-h] [<input object 0>] [<input object 1>]\n"
      "Options:\n"
      "  -b  Shell breaking control value.\n"
      "  -d  Debug output file.\n"
      "  -D  Debug flag:\n"
      "        0  no debug output.\n"
      "        1  untransformed decomposed source model.\n"
      "  -E  Tolerance in mean registration metric value.\n"
      "  -h  Prints this usage information.\n"
      "  -o  Output file name.\n"
      "  -t  Initial affine transform.\n"
      "  -P  Minimum number of simplices per matched shell segment, with a\n"
      "      pair of correspondence points possibly being generated per\n"
      "      matched shell segment.\n"
      "  -i  Maximum number of iterations.\n"
      "  -s  Minmum number of simplicies.\n"
      "  -A  Maximum angle (degrees) from a global transformation.\n"
      "  -S  Maximum displacement from a global transformed position.\n"
      "  -F  Maximum deformation from a global transformation.\n"
      "  -I  Implausibility threshold for rejecting implausible\n"
      "      correspondence points which should be greater than zero,\n"
      "      although the useful range is probably [0.5-2.5]. Higher\n"
      "      values allow more implausible matches to be returned.\n"
      "  -N  Number of match points in neighbourhood when checking the\n"
      "      plausibility of the correspondence points.\n"
      "Reads a pair of contours and computes a set of correspondence points\n"
      "using an ICP based matching algorithm. The correspondence points are\n"
      "ranked by plausibility, with the most plausible first\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
