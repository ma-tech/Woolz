#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzShellFilter3D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzShellFilter3D.c
* \author       Bill Hill
* \date         November 2020
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Computes a new Woolz object by passing all planes of a 3D
* 		object through an external (shell) filter.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzshellfilter3d "WlzShellFilter3D"
*/

#include <stdio.h>
#include <float.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <Wlz.h>
#include <WlzExtFF.h>

/* HACK TODO copy in the doc below. */
/*!
\ingroup BinWlz
\defgroup wlzshellfilter3d WlzShellFilter3D
\par Name
WlzShellFilter3D
\par Synopsis
\verbatim
WlzShellFilter3D [-g] [-h] [-L] [-R] [-v] [-a<pit><yaw>[,<roll>]]
                 [-d<min[,max>]] [-e<shell cmd>] [-f<x,y,z>] [-F<fmt>]
		 [-m<mode>] [-p<env par>[,<env par>...]] [-s<scale>]
		 [-t<tmp file>,<rmp file>] [-u<x,y,z>] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-g</b></td>
    <td>Shift planes to the origin before writing to the filter and
        then restore the offset, default false.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Use linear interpolation instead of nearest neighbour when
        cutting planes, default false.</td>
  </tr>
  <tr> 
    <td><b>-R</b></td>
    <td>Angles are in radians not degrees, default false.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose output messages.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>The cutting angles pitch, yaw and optionally roll.
        If roll is defined then the mode is absolute.
	Default 0, 0, 0.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Distance range, default object min, object max.</td>
  </tr>
  <tr> 
    <td><b>-e</b></td>
    <td>Shell filter command to be executed. If not supplied then
        the filter is simply an identity operation.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file, default standard output.</td>
  </tr>
  <tr> 
    <td><b>-F</b></td>
    <td>File format for shell filter command input/output, default wlz.
        The list of file format strings recognized can be seen by using
	-F help.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Fixed point, default 0, 0, 0.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Viewing mode, possible values:
  </tr>
  <tr>
    <td> </td>
    <td>
      <table width="500" border="0">
        <tr> <td><b>Parameter value</b></td> <td><b>Viewing mode</b></td> </tr>
        <tr> <td>0</td> <td>up-is-up, default</td> </tr>
        <tr> <td>1</td> <td>statue</td> </tr>
        <tr> <td>2</td> <td>absolute</td> </tr>
      </table>
    Default 0.
    </td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Comma separated list of additional shell environment parameters
        each with the format variable_name=variable_value.</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Voxel size rescaling flags:
  </tr>
  <tr>
    <td> </td>
    <td>
      <table width="500" border="0">
	<tr> <td>bit 1 set</td><td>use voxel-size rescaling</td> </tr>
	<tr> <td>bit 2 set</td><td>use global scaling</td> </tr>
      </table>
    Default 1.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Scaling, default 1.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Temporary files to use for input and output to and from the shell
        command. These are required whenever a shell command is specified
        and must be regular files can not (for example) the standard input
        or standard output. The file format must be set appropriately
        for the temporary files and shell command. Default not set.</td>
  </tr>
  <tr>
    <td><b>-u</b></td>
    <td>Up vector only used for up-is-up mode. Default 0, 0, -1.</td>
  </tr>
</table>
\par Description
Applies the given shell command to each plane of the given object
where the planes are cut parallel to each other, are consecutive
and are defined by the cutting angles, mode and distance(s).
Ranges of planes may be specified by the fixed point and distance
range. By default the input object is read from stdin and the output
object is written to stdout.
\par Examples
\verbatim
WlzShellFilter3D -o out.wlz -a 90,0 -F tif -t t1.tif,t2.tif \
                 -e 'convert -negate t1.tif t2.tif' in.wlz
\endverbatim
Inverts all grey values of the 3D domain object read from in.wlz
and makes the corresponding y-z planes x-y planes in the 3D output
file out.wlz. Each plane is written to t1.tif, the ImageMagick
convert command is executed (to invert grey values) and then the
file t2.tif is read to form an x-y plane in the 3D output object.
\par File
\ref WlzShellFilter3D.c "WlzShellFilter3D.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref WlzGetSectionFromObject "WlzGetSectionFromObject(3)"
\ref WlzEffReadObj "WlzEffReadObj(3)"
\ref WlzEffWriteObj "WlzEffWriteObj(3)"
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
  int		option,
		mode = 0,
		shift = 0,
		envCnt = 0,
		useRadians = 0,
		voxRescale = 0,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  char		*inStr,
		*outStr,
	 	*cmd = NULL,
		*envBuf = NULL,
		*fFmtStr = "wlz",
		*tmpFile[2] = {0};
  double	scale = 1.0;
  char		**envVar = NULL,
		**envVal = NULL;
  double	angles[3] = {0},
  		distRange[2] = {-DBL_MAX, DBL_MAX};
  WlzEffFormat	fFmt = WLZEFF_FORMAT_WLZ;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzDVertex3	fixed = {0.0, 0.0, 0.0},
  		up = {0.0, 0.0, -1.0};
  WlzThreeDViewStruct *viewStr=NULL;
  WlzDomain	outDom = {0};
  WlzValues	outVal = {0};
  struct timeval times[3];
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  static char   optList[] = "ghLRva:d:e:f:F:m:o:p:r:t:u:";
  const char    *falseTrue[2] = {"false", "true"};
  const char    inStrDef[] = "-",
		outStrDef[] = "-";

  opterr = 0;
  inStr = (char *)inStrDef;
  outStr = (char *)outStrDef;
  while((usage == 0) && (errNum == WLZ_ERR_NONE) &&
        ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'g':
        shift = 1;
	break;
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
        break;
      case 'R':
        useRadians = 1;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'a':
        if(sscanf(optarg, "%lg,%lg,%lg", angles, angles + 1, angles + 2) < 2)
        {
	  usage = 1;
	}
	break;
      case 'd':
	{
	  int	pc;

	  if((pc = sscanf(optarg, "%lg,%lg", distRange, distRange + 1)) < 1)
	  {
	    usage = 1;
	  }
	  else if(pc < 2)
	  {
	    int di;
	    
	    di = (distRange[0] > -DBL_MAX);
	    distRange[di] = distRange[!di];
	  }
	}
	break;
      case 'e':
        cmd = optarg;
        break;
      case 'f':
        if(sscanf(optarg, "%lg,%lg,%lg",
	          &(fixed.vtX), &(fixed.vtY), &(fixed.vtZ)) < 1)
        {
	  usage = 1;
	}
	break;
      case 'F':
	fFmtStr = optarg;
        if((fFmt = WlzEffStringExtToFormat(fFmtStr)) == 0)
        {
	  usage = 1;
	}
	break;
      case 'm':
        if((sscanf(optarg, "%d", &mode) != 1) || (mode < 0) || (mode > 2))
	{
	  usage = 1;
	}
	break;
      case 'o':
        outStr = optarg;
	break;
      case 'p':
        {
	  char	*c = optarg;

          envCnt = 0;
	  while(*c)
	  {
	    envCnt += (*c == ',');
	    ++c;
	  }
	  ++envCnt;
	}
	if(((envVar = (char **)
	              AlcMalloc(sizeof(char *) * 2 * envCnt)) == NULL) ||
	   ((envBuf = (char *)
	              AlcMalloc(sizeof(char) * strlen(optarg))) == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  int	i = 0;
	  char 	*buf,
	  	*pair;

	  buf = envBuf;
	  envVal = envVar + envCnt;
	  pair = WlzStringWhiteSpSkipLeading(strtok(optarg, ","));
	  while(pair)
	  {
	    char *sep;

	    (void)strcpy(buf, pair);
	    sep = strchr(buf, '=');
	    envVar[i] = buf;
	    if(sep)
	    {
	      *sep = '\0';
	      envVal[i] = envVar[i] + strlen(envVar[i]) + 1;
	    }
	    else
	    {
	      envVal[i] = envVar[i] + strlen(envVar[i]);
	    }
	    ++i;
	    buf += strlen(pair) + 1;
	    pair = WlzStringWhiteSpSkipLeading(strtok(NULL, ","));
	  }
	}
	break;
      case 'r':
        if(sscanf(optarg, "%d", &voxRescale) != 1)
	{
	  usage = 1;
	}
	break;
      case 't':
	tmpFile[0] = WlzStringWhiteSpSkipLeading(strtok(optarg, ","));
	tmpFile[1] = WlzStringWhiteSpSkipLeading(strtok(NULL, ","));
	if((tmpFile[0] == NULL) || (strlen(tmpFile[0]) < 1) ||
	   (strcmp(tmpFile[0], "-") == 0) ||
	   (tmpFile[1] == NULL) || (strlen(tmpFile[1]) < 1) ||
	   (strcmp(tmpFile[1], "-") == 0))
	{
	  usage = 1;
	}
	break;
      case 'u':
        if(sscanf(optarg, "%lg,%lg,%lg", &(up.vtX), &(up.vtY), &(up.vtZ)) != 3)
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
  if(cmd && ((tmpFile[0] == NULL)  || (tmpFile[1] == NULL)))
  {
    usage = 1;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    ok = 0;
    (void )WlzStringFromErrorNum(errNum, &errMsgStr);
    (void )fprintf(stderr,
		   "%s Failed to parse command line, %s.\n",
		   argv[0],
		   errMsgStr);
  }
  else
  {
    ok = (usage == 0);
  }
  if(ok)
  {
    int		i;

    if(verbose)
    {
      (void )fprintf(stderr,
		     "%s Environment parameters to pass to shell are:\n",
		     argv[0]);
    }
    for(i = 0; i < envCnt; ++i)
    {
      if(verbose)
      {
        (void )fprintf(stderr, "\t%s=%s\n", envVar[i], envVal[i]);
      }
      (void )setenv(envVar[i], envVal[i], 1);
    }
  }
  if(ok)
  {
    if((inStr == NULL) || (*inStr == '\0') ||
       (outStr == NULL) || (*outStr == '\0'))
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
        inStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    FILE	*fP = NULL;

    if((inStr == NULL) ||
       (*inStr == '\0') ||
       ((fP = (strcmp(inStr, "-")?
              fopen(inStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inStr);
    }
    if(fP && strcmp(inStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    if(inObj->type != WLZ_3D_DOMAINOBJ)
    {
      ok = 0;
      errNum = WLZ_ERR_OBJECT_TYPE;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s 3D domain object required, %s.\n",
      		     argv[0],
		     errMsgStr);
      
    }
  }
  if(ok)
  {
    if((viewStr = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum)) != NULL)
    {
      double  dtr;

      dtr = (useRadians)? 1.0: ALG_M_PI / 180.0;
      viewStr->phi   = dtr * angles[0];
      viewStr->theta = dtr * angles[1];
      viewStr->zeta  = dtr * angles[2];
      viewStr->dist  = 0.0;
      viewStr->fixed = fixed;
      viewStr->up    = up;
      viewStr->view_mode = mode;
      viewStr->scale = scale;
      viewStr->voxelRescaleFlg = voxRescale;
      errNum = WlzInit3DViewStruct(viewStr, inObj);
      if(errNum == WLZ_ERR_NONE)
      {
	distRange[0] = ALG_CLAMP(distRange[0],
	                         viewStr->minvals.vtZ, viewStr->maxvals.vtZ);
	distRange[1] = ALG_CLAMP(distRange[1],
	                         distRange[0], viewStr->maxvals.vtZ);
	if(verbose)
        {
	  (void )fprintf(stderr,
			 "%s: Section parameters are\n"
			 "  fixed = %lg,%lg,%lg\n"
			 "  theta = %lg\n"
			 "  phi = %lg\n"
			 "  zeta = %lg\n"
			 "  dist = %lg\n"
			 "  scale = %lg\n"
			 "  viewMode = %s\n"
			 "  up = %lg,%lg,%lg\n",
			 argv[0],
			 viewStr->fixed.vtX, fixed.vtY, fixed.vtZ,
			 viewStr->theta,
			 viewStr->phi,
			 viewStr->zeta,
			 viewStr->dist,
			 viewStr->scale,
			 WlzStringFromThreeDViewMode(viewStr->view_mode, NULL),
			 viewStr->up.vtX, viewStr->up.vtY, viewStr->up.vtZ);
	  (void )fprintf(stderr,
                         "%s: Distance range is %lg - %lg\n",
			 argv[0],
			 distRange[0], distRange[1]);
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to initialise 3D view structure, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    WlzIBox3	outBB;

    outBB.xMin = (int )floor(viewStr->minvals.vtX),
    outBB.xMax = (int )floor(viewStr->maxvals.vtX),
    outBB.yMin = (int )floor(viewStr->minvals.vtY),
    outBB.yMax = (int )floor(viewStr->maxvals.vtY),
    outBB.zMin = (int )floor(distRange[0]),
    outBB.zMax = (int )floor(distRange[1]),
    outDom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
    			          outBB.zMin, outBB.zMax,
				  outBB.yMin, outBB.yMax,
			          outBB.xMin, outBB.xMax, &errNum);
    if((errNum == WLZ_ERR_NONE) && (inObj->values.core != NULL) &&
       (inObj->values.core->type == WLZ_VOXELVALUETABLE_GREY))
    {
      outVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
    			               outBB.zMin, outBB.zMax,
				       inObj->values.vox->bckgrnd, NULL,
				       &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      outObj = WlzMakeMain(inObj->type, outDom, outVal, NULL, NULL, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(outObj == NULL)
      {
        (void )WlzFreePlaneDomain(outDom.p);
	(void )WlzFreeVoxelValueTb(outVal.vox);
      }
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to create output object, %s.\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    int		i,
    		nPlanes;

    nPlanes = outDom.p->lastpl - outDom.p->plane1 + 1;
    for(i = 0; i < nPlanes; ++i)
    {
      WlzObject	*iObj2 = NULL,
      		*oObj2 = NULL;
      WlzIBox2   iBB2 = {0};

      viewStr->dist = outDom.p->plane1 + i;
      errNum = WlzInit3DViewStruct(viewStr, inObj);
      if(errNum == WLZ_ERR_NONE)
      {
	iObj2 = WlzAssignObject(
	        WlzGetSubSectionFromObject(inObj, NULL, viewStr, interp, NULL,
					   &errNum), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	iBB2 = WlzBoundingBox2I(iObj2, &errNum);
      }
      if((errNum == WLZ_ERR_NONE) && shift)
      {
        WlzObject *tObj;

	tObj = WlzAssignObject(
               WlzShiftObject(iObj2, -(iBB2.xMin), -(iBB2.yMin), 0,
	                      &errNum), NULL);
        (void )WlzFreeObj(iObj2);
	iObj2 = tObj;
      }
      if(cmd == NULL)
      {
        oObj2 = iObj2;
	iObj2 = NULL;
      }
      else
      {
	FILE  *fP = NULL;

	errNum = WLZ_ERR_WRITE_EOF;
	if((fP = fopen(tmpFile[0], "w")) != NULL)
	{
	  errNum = WlzEffWriteObj(fP, tmpFile[0], iObj2, fFmt);
	  if(fclose(fP) && (errNum == WLZ_ERR_NONE))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if(system(cmd))
	  {
	    errNum = WLZ_ERR_UNSPECIFIED;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WLZ_ERR_READ_EOF;
          if((fP = (strcmp(tmpFile[1], "-")?
                   fopen(tmpFile[1], "r"): stdin)) != NULL)
	  {
	    oObj2 = WlzAssignObject(
		    WlzEffReadObj(fP, tmpFile[1], fFmt, 0, 0, 0,
				  &errNum), NULL);
	    if(fclose(fP) && (errNum == WLZ_ERR_NONE))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	  }
	}
      }
      (void )WlzFreeObj(iObj2);
      if((errNum == WLZ_ERR_NONE) && shift)
      {
        WlzObject *tObj;

	tObj = WlzAssignObject(
               WlzShiftObject(oObj2, iBB2.xMin, iBB2.yMin, 0, &errNum), NULL);
        (void )WlzFreeObj(oObj2);
	oObj2 = tObj;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	outDom.p->domains[i] = WlzAssignDomain(oObj2->domain, NULL);
	outVal.vox->values[i] = WlzAssignValues(oObj2->values, NULL);
      }
      (void )WlzFreeObj(oObj2);
    }
  }
  if(ok)
  {
    errNum = WlzStandardPlaneDomain(outDom.p, outVal.vox);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		     "%s Failed to standardize output object domain, %s.\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if(verbose)
    {
      gettimeofday(times + 0, NULL);
    }
    if(verbose)
    {
      gettimeofday(times + 1, NULL);
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
                     "%s: Elapsed time for WlzShellFilter3D() %gus\n",
                     argv[0], (1.0e06 * times[2].tv_sec) + times[2].tv_usec);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to filter object, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    FILE	*fP = NULL;

    if((fP = (strcmp(outStr, "-")?  fopen(outStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outStr);
    }
    else
    {
      errNum = WlzWriteObj(fP, outObj);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to write object to file %s\n",
		       *argv, outStr);
      }
    }
    if(fP && strcmp(outStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  AlcFree(envBuf);
  AlcFree(envVar);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  (void )WlzFree3DViewStruct(viewStr);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-g] [-h] [-L] [-R] [-v] [-a<pit><yaw>[,<roll>]]\n"
    "\t\t[-d<min[,max>]] [-e<shell cmd>] [-f<x,y,z>] [-F<fmt>]\n"
    "\t\t[-m<mode>] [-p<env par>[,<env par>...]] [-s<scale>]\n"
    "\t\t[-t<tmp file>,<rmp file>] [-u<x,y,z>] [<input object>]\n"
    "Applies the given shell command to each plane of the given object\n"
    "where the planes are cut parallel to each other, are consecutive\n"
    "and are defined by the cutting angles, mode and distance(s).\n"
    "Ranges of planes may be specified by the fixed point and distance\n"
    "range. By default the input object is read from stdin and the output\n"
    "object is written to stdout.\n"
    "Version: %s\n"
    "Options:\n"
    "  -g  Shift planes to the origin before writing to the filter and\n"
    "      then restore the offset, set to %s.\n"
    "  -h  Help, prints this usage message.\n"
    "  -R  Angles are in radians not degrees, set to %s.\n"
    "  -v  Verbose output messages.\n"
    "  -a  Cutting angles, set to %g, %g, %g.\n"
    "  -d  Distance range, default object min, object max.\n"
    "  -e  Shell filter command to be executed. If not supplied then the\n"
    "      filter is simply an identity operation.\n"
    "  -o  Output file, default standard output.\n"
    "  -F  File format for shell filter command input/output, set to %s.\n"
    "      The list of file format strings recognized can be seen by using\n"
    "      -F help.\n"
    "  -f  Fixed point, set to %g, %g, %g.\n"
    "  -L  Use linear interpolation instead of nearest neighbour when\n"
    "      cutting planes, set to %s.\n"
    "  -m  Viewing mode, possible values:\n"
    "        0 - up-is-up\n"
    "        1 - statue\n"
    "        2 - absolute\n"
    "      set to %d.\n"
    "  -p  Comma separated list of additional shell environment parameters\n"
    "      each with the format <variable name>=<value>.\n"
    "  -r  Voxel size rescaling flags:\n"
    "        bit 1 set - use voxel-size rescaling,\n"
    "        bit 2 set - enable global scaling,\n"
    "      set to %d\n"
    "  -s  Scaling, set to %g.\n"
    "  -t  Temporary files to use for input and output to and from the shell\n"
    "      command. These are required whenever a shell command is specified\n"
    "      and must be regular files can not (for example) the standard input\n"
    "      or standard output. The file format must be set appropriately\n"
    "      for the temporary files and shell command. Set to %s,%s\n"
    "  -u  Up vector only used for up-is-up mode, set to  %g, %g, %g.\n"
    "Example:\n"
    "  %s -o out.wlz -a 90,0 -F tif -t t1.tif,t2.tif \n"
    "\t\t-e 'convert -negate t1.tif t2.tif' in.wlz\n"
    "Inverts all grey values of the 3D domain object read from in.wlz\n"
    "and makes the corresponding y-z planes x-y planes in the 3D output\n"
    "file out.wlz. Each plane is written to t1.tif, the ImageMagick\n"
    "convert command is executed (to invert grey values) and then the\n"
    "file t2.tif is read to form an x-y plane in the 3D output object.\n",
    argv[0],
    WlzVersion(),
    falseTrue[shift && 1],
    falseTrue[useRadians && 1],
    angles[0], angles[1], angles[2],
    fFmtStr,
    fixed.vtY, fixed.vtY, fixed.vtZ,
    falseTrue[interp != WLZ_INTERPOLATION_NEAREST],
    mode,
    voxRescale,
    scale,
    (tmpFile[0])? tmpFile[0]: "NULL", (tmpFile[1])? tmpFile[1]: "NULL",
    up.vtX, up.vtY, up.vtZ,
    argv[0]);
  }
  if(fFmt == WLZEFF_FORMAT_NONE)
  {
    (void )fprintf(stderr,
    "Recognised file formats are:\n"
    "  Description                                       Extension\n"
    "  ***********                                       *********\n%s\n",
    WlzEffFormatTable(2, 50, 10, NULL));
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
