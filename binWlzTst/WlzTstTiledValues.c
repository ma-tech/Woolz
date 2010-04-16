#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTstTiledValues_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTstTiledValues.c
* \author       Bill Hill
* \date         April 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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
* \brief	Test for creating and accessing objects with tiled values
* 		tables.
* \ingroup	Tst
*/


#include <sys/time.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		option,
  		ok = 1,
		usage = 0,
		dim = 2,
		section = 0,
		timer = 0;
  double	yaw = 0.0,
  		pitch = 0.0,
		roll =  0.0,
		dist = 0.0,
		scale = 1.0;
  WlzPixelV	bgdV;
  WlzGreyType	gType;
  WlzDVertex3	up,
  		fixed;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzThreeDViewMode mode = WLZ_UP_IS_UP_MODE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
  		*tlObj = NULL;
  char		*inFileStr,
		*outFileStr,
  		*secFileStr;
  struct timeval times[3];
  const char	*errMsg;
  const size_t	tlSz = 4096;
  static char	optList[] = "hsto:S:",
  		inFileStrDef[] = "-";

  opterr = 0;
  outFileStr = NULL;
  inFileStr = inFileStrDef;
  secFileStr = inFileStrDef;
  gType = WLZ_GREY_UBYTE;
  bgdV.type = WLZ_GREY_UBYTE;
  bgdV.v.ubv = 0;
  up.vtX = up.vtY = 0.0; up.vtZ = -1.0;
  fixed.vtX = fixed.vtY = fixed.vtZ = 0.0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 's':
        section = 1;
	break;
      case 'S':
        secFileStr = optarg;
	break;
      case 't':
        timer = 1;
	break;
      case 'h': /* FALLTHROUGH */
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
      inFileStr = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    fP = NULL;
    errNum = WLZ_ERR_READ_EOF;
    if(((fP = (strcmp(inFileStr, "-")? fopen(inFileStr, "r"):
                                       stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
		     *argv, inFileStr, errMsg);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    if(inObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if(inObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      switch(inObj->type)
      {
        case WLZ_2D_DOMAINOBJ:
	  dim = 2;
	  break;
        case WLZ_3D_DOMAINOBJ:
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
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: invalid object read from file %s (%s).\n",
		     *argv, inFileStr, errMsg);
    }
  }
  if(ok)
  {
    if(timer)
    {
      gettimeofday(times + 0, NULL); 
    }
    tlObj = WlzMakeTiledValuesFromObj(inObj, tlSz, 1, gType, bgdV, &errNum);
    if(timer)
    {
      gettimeofday(times + 1, NULL); 
      timersub(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
                     "%s: Elapsed time for WlzMakeTiledValuesFromObj() %gs\n",
                     *argv,
		     times[2].tv_sec + (0.000001 * times[2].tv_usec));

    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to create object with tiled values (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok && (outFileStr != NULL))
  {
    if(((fP = (strcmp(outFileStr, "-")? fopen(outFileStr, "w"):
                                        stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, tlObj)) != WLZ_ERR_NONE))
    
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write tiled object to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
  }
  if(ok && section && (dim == 3))
  {
    WlzObject 	*secObj = NULL;
    WlzThreeDViewStruct *view = NULL;

    view = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      view->theta = yaw * WLZ_M_PI / 180.0;
      view->phi = pitch * WLZ_M_PI / 180.0;
      view->zeta = roll * WLZ_M_PI / 180.0;
      view->dist = dist;
      view->fixed = fixed;
      view->up = up;
      view->view_mode = mode;
      view->scale = scale;
      errNum = WlzInit3DViewStruct(view, tlObj);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(timer)
      {
	gettimeofday(times + 0, NULL); 
      }
      secObj = WlzGetSubSectionFromObject(tlObj, NULL, view, interp,
      					  NULL, &errNum);
      if(timer)
      {
	gettimeofday(times + 1, NULL); 
	timersub(times + 1, times + 0, times + 2);
	(void )fprintf(stderr,
		     "%s: Elapsed time for WlzGetSubSectionFromObject() %gs\n",
		     *argv,
		     times[2].tv_sec + (0.000001 * times[2].tv_usec));
      }
      fP = NULL;
      errNum = WLZ_ERR_WRITE_EOF;
      if(((fP = (strcmp(secFileStr, "-")? fopen(secFileStr, "w"):
					  stdout)) == NULL) ||
	 ((errNum = WlzWriteObj(fP, secObj)) != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to write section object to file %s (%s).\n",
		       *argv, secFileStr, errMsg);
      }
      if(fP && strcmp(secFileStr, "-"))
      {
	(void )fclose(fP);
      }
    }
    (void )WlzFree3DViewStruct(view);
    (void )WlzFreeObj(secObj);
  }
  (void )WlzFreeObj(tlObj);
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<output object>] [-h] [-o <file>] [-s] [-S <file>] [-t]\n"
    "                  [<input object>]\n"
    "Copied the input object to an object with tiled values.\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output tiled object.\n"
    "  -s  Cut section from tiled object.\n"
    "  -S  Output file for section object.\n"
    "  -t  Output timing information.\n");
  }
  return(!ok);
}
