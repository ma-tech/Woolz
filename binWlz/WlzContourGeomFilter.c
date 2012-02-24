#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzContourGeomFilter_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzContourGeomFilter.c
* \author       Bill Hill
* \date         September 2002
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
* \brief	Filters a contour model.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcontourgeomfilter "WlzContourGeomFilter"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzcontourgeomfilter WlzContourGeomFilter
\par Name
WlzContourGeomFilter - filters a contour model.
\par Synopsis
\verbatim
WlzContourGeomFilter [-h] [-o<output file>]
        [-l#] [-m#] [-n#]
        [-p#] [-s#] [-P#] [-S#] [<input file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-N</b></td>
    <td>Allow non manifild vertices to be filtered.</td>
  </tr>
  <tr>
    <td><b>-l</b></td>
    <td>Filter parameter \f$\lambda\f$.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Filter parameter \f$\mu\f$.</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>The number of itterations.
        If the number of itterations is not given, then the pass and stop
	band ripple determine the number of itterations.</td>
  </tr>
  <tr>
    <td><b>-p</b></td>
    <td>The band pass frequency parameter (kPB)</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>The band stop frequency parameter (kSB).</td>
  </tr>
  <tr>
    <td><b>-P</b></td>
    <td>The maximum pass band ripple (dPB)</td>
  </tr>
  <tr>
    <td><b>-S</b></td>
    <td>The maximum stop band ripple (dSB)</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
</table>

\par Description
Low pass filters the geometry of a contour object to smooth without
shrinkage.

For simple use to smooth a geometric model from a 2D or 3D contour
don't use lambda and mu. Instead use: Values of kPB in the range
[0.0100-0.1000], values of kSB in the range [1.000-1.9999]. The default
ripple values are probably acceptable for most simple filters.
Simple example, which filters a contour smoothing out voxel edges:
    prompt% WlzContourGeomFilter -o out.wlz -p 0.1 -s 1.5 -n 100 in.wlz
    This example passes the geometric model from the contour object
    in in.wlz through a low pass filter with a transition band starting
    at 0.1 and ending with the stop band at 1.5. A maximum of 100
    itterations are performed. If the output from the filter was not
sufficiently smooth then the ripple values could be decreased.
If you want to design a more complex filter or one that is not a
low pass (smoothing) filter then see either:
Gabriel Taubin <i>Curve and Surface Smoothing without Shrinkage
International Conference on Computer Vision (ICCV'95) Conference
Procedings 827-857 1995</i>
or
Gabriel Taubin <i>A Signal Processing Approach to Fair Surface
Design Computer Graphics, 29, 351-358 1995</i>.
Use the -h flag to see the default parameter values.
\par Examples
\verbatim
WlzContourGeomFilter -o smooth.wlz in.wlz
\endverbatim
Smooths the input contour model in.wlz and writes the smoothed result to
smooth.wlz.

\par File
\ref WlzContourGeomFilter.c "WlzContourGeomFilter.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzGMFilterGeomLPParam "WlzGMFilterGeomLPParam(3)"
\ref WlzGMFilterGeomLPLM "WlzGMFilterGeomLPLM(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <Wlz.h>


typedef enum _WLZCtrFilterFlags
{
  WLZWLZ_CTRFFLAGS_NONE     = 0,
  WLZWLZ_CTRFFLAGS_LAMBDA   = (1 << 0),
  WLZWLZ_CTRFFLAGS_MU       = (1 << 1),
  WLZWLZ_CTRFFLAGS_NITR     = (1 << 2),
  WLZWLZ_CTRFFLAGS_KPB      = (1 << 3),
  WLZWLZ_CTRFFLAGS_KSB      = (1 << 4),
  WLZWLZ_CTRFFLAGS_DPB      = (1 << 5),
  WLZWLZ_CTRFFLAGS_DSB      = (1 << 6)
} WLZCtrFilterFlags;

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		tI0,
  		ok = 1,
  		option,
  		usage = 0,
		nItr = 200,
		nonMan = 0;
  unsigned int	flags = 0;
  double 	lambda,
  		mu,
		kPB = 0.1000,
		kSB = 1.0000,
		dPB = 0.2500,
		dSB = 0.1000;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  static char   optList[] = "hNo:l:m:n:p:s:P:S:";
  const char    inObjFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-",
		setMsg[] =    "set on command line",
		notSetMsg[] = "not set / default  ";

  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'N':
        nonMan = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'l':
	if((sscanf(optarg, "%lg", &lambda) != 1) ||
	   (fabs(lambda) < DBL_EPSILON) || (fabs(lambda) > 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	flags |= WLZWLZ_CTRFFLAGS_LAMBDA;
	break;
      case 'm':
	if((sscanf(optarg, "%lg", &mu) != 1) ||
	   (fabs(mu) < DBL_EPSILON) || (fabs(mu) > 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	flags |= WLZWLZ_CTRFFLAGS_MU;
	break;
      case 'n':
	if((sscanf(optarg, "%d", &nItr) != 1) || (nItr < 0))
	{
	  usage = 1;
	  ok = 0;
	}
	flags |= WLZWLZ_CTRFFLAGS_NITR;
	break;
      case 'p':
	if((sscanf(optarg, "%lg", &kPB) != 1) ||
	   (kPB < DBL_EPSILON) || (kPB > 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	flags |= WLZWLZ_CTRFFLAGS_KPB;
	break;
      case 's':
	if((sscanf(optarg, "%lg", &kSB) != 1) ||
	   (kSB < 1.0) || (kSB > 2.0))
	{
	  usage = 1;
	  ok = 0;
	}
	flags |= WLZWLZ_CTRFFLAGS_KSB;
	break;
      case 'P':
	if((sscanf(optarg, "%lg", &dPB) != 1) ||
	   (dPB < DBL_EPSILON) || (dPB > 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	flags |= WLZWLZ_CTRFFLAGS_DPB;
	break;
      case 'S':
	if((sscanf(optarg, "%lg", &dSB) != 1) ||
	   (dSB < DBL_EPSILON) || (dSB > 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	flags |= WLZWLZ_CTRFFLAGS_DSB;
	break;
      case 'h':
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    tI0 = WLZWLZ_CTRFFLAGS_LAMBDA | WLZWLZ_CTRFFLAGS_MU;
    if((flags & tI0) != tI0)
    {
      errNum = WlzGMFilterGeomLPParam(&lambda, &mu, &tI0,
      				      kPB, kSB, dPB, dSB);
    }
    if((flags & WLZWLZ_CTRFFLAGS_NITR) == 0)
    {
      nItr = tI0;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      "%s: Failed to compute filter parameters lambda and mu from the\n"
      "             band pass and band stop frequency parameters (%s).\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if(obj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    if(obj->type != WLZ_CONTOUR)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if(obj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      errNum = WlzGMFilterGeomLPLM(obj->domain.ctr->model, lambda, mu, nItr,
      				   nonMan);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to filter the geometric model. (%s).\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outObjFileStr, "-")?
	     fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outObjFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output file>]\n"
    "        [-l#] [-m#] [-n#]\n"
    "        [-p#] [-s#] [-P#] [-S#] [<input file>]\n"
    "  -h  Output this usage message.\n"
    "  -N  Allow non manifild vertices to be filtered.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -l  Filter parameter lambda.\n"
    "  -m  Filter parameter mu.\n"
    "  -n  The number of itterations. If the number of itterations\n"
    "      is not given, then the pass and stop band ripple determine\n"
    "      the number of itterations.\n"
    "  -p  The band pass frequency parameter (kPB).\n"
    "  -s  The band stop frequency parameter (kSB).\n"
    "  -P  The maximum pass band ripple (dPB).\n"
    "  -S  The maximum stop band ripple (dSB).\n"
    "For simple use to smooth a geometric model from a 2D or 3D contour\n"
    "don't use lambda and mu. Instead use: Values of kPB in the range\n"
    "[0.0100-0.1000], values of kSB in the range [1.000-1.9999]. The default\n"
    "ripple values are probably acceptable for most simple filters.\n"
    "Simple example, which filters a contour smoothing out voxel edges:\n"
    "    prompt%% WlzContourGeomFilter -o out.wlz -p 0.1 -s 1.5 -n 100 in.wlz\n"
    "This example passes the geometric model from the contour object\n"
    "in in.wlz through a low pass filter with a transition band starting\n"
    "at 0.1 and ending with the stop band at 1.5. A maximum of 100\n"
    "itterations are performed. If the output from the filter was not\n"
    "sufficiently smooth then the ripple values could be decreased.\n"
    "If you want to design a more complex filter or one that is not a\n"
    "low pass (smoothing) filter then see either:\n"
    "    Gabriel Taubin. Curve and Surface Smoothing without Shrinkage.\n"
    "    International Conference on Computer Vision (ICCV'95) Conference\n"
    "    Procedings, 827-857, 1995.\n"
    "or\n"
    "    Gabriel Taubin. A Signal Processing Approach to Fair Surface\n"
    "    Design. Computer Graphics, 29, 351-358, 1995.\n",
    argv[0]);
    (void )fprintf(stderr,
		   "The parameter values are:\n"
		   "    %s kPB =    %g\n"
		   "    %s kSB =    %g\n"
		   "    %s dPB =    %g\n"
		   "    %s dSB =    %g\n"
		   "    %s lambda = %g\n"
		   "    %s mu =     %g\n"
		   "    %s nItr =   %d\n",
		   (flags & WLZWLZ_CTRFFLAGS_KPB)? setMsg: notSetMsg, kPB,
		   (flags & WLZWLZ_CTRFFLAGS_KSB)? setMsg: notSetMsg, kSB,
		   (flags & WLZWLZ_CTRFFLAGS_DPB)? setMsg: notSetMsg, dPB,
		   (flags & WLZWLZ_CTRFFLAGS_DSB)? setMsg: notSetMsg, dSB,
		   (flags & WLZWLZ_CTRFFLAGS_LAMBDA)? setMsg: notSetMsg, lambda,
		   (flags & WLZWLZ_CTRFFLAGS_MU)? setMsg: notSetMsg, mu,
		   (flags & WLZWLZ_CTRFFLAGS_NITR)? setMsg: notSetMsg, nItr);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
