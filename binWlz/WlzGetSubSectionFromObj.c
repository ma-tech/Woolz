#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGetSubSectionFromObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzGetSubSectionFromObj.c
* \author       Bill Hill
* \date         March 2009
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
* \brief	Cuts a section from an object with the domain of a given
* 		object.
* \ingroup	BinWlz
*/

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

static int      		WlzGSFOReadTriple(
				  char *str,
				  double *d0,
				  double *d1,
				  double *d2);

int		main(int argc, char *argv[])
{
  int		tI,
		ok,
  		option,
  		usage = 0;
  WlzObject	*inObj = NULL,
  		*refObj = NULL,
		*outObj = NULL;
  WlzThreeDViewStruct *view = NULL;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;
  char		*inFileStr,
  		*outFileStr,
  	        *inRefFileStr;
  const char	*errMsg;
  FILE          *fP = NULL;
  static char   optList[] = "a:f:d:h:m:u:o:r:",
  		fileStrDef[] = "-";


  ok = 1;
  opterr = 0;
  errMsg = "";
  inFileStr = inRefFileStr = fileStrDef;
  if((view = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum)) == NULL)
  {
    ok = 0;
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
                   "%s: Failed to create view structure (%s)\n",
		   *argv, errMsg);
  }
  if(ok)
  {
    view->dist = 0.0;
    view->scale = 1.0;
    view->ref_obj = NULL;
    view->view_mode = WLZ_UP_IS_UP_MODE;
    view->phi = view->theta = view->zeta = 0.0;
    view->up.vtX = view->up.vtY = 0.0; view->up.vtZ = 1.0;
    view->fixed.vtX = view->fixed.vtY = view->fixed.vtZ = 0.0;
    while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
    {
      switch(option)
      {
	case 'a':
	  usage = WlzGSFOReadTriple(optarg,
				    &(view->phi),
				    &(view->theta),
				    &(view->zeta)) < 1;
	  view->phi *= ALG_M_PI / 180;
	  view->theta *= ALG_M_PI / 180;
	  view->zeta *= ALG_M_PI / 180;
	  break;
	case 'f':
	  break;
	case 'd':
	  usage = sscanf(optarg, "%lg", &(view->dist)) != 1;
	  break;
	case 'm':
	  if(WlzStringMatchValue(&tI, optarg,
				 "up-is-up", WLZ_UP_IS_UP_MODE,
				 "statue", WLZ_STATUE_MODE,
				 "absolute", WLZ_ZETA_MODE,
				 NULL))
	  {
	    view->view_mode = (WlzThreeDViewMode )tI;
	  }
	  else
	  {
	    usage = 1;
	  }
	  break;
	case 'u':
	  usage = WlzGSFOReadTriple(optarg,
				    &(view->up.vtX),
				    &(view->up.vtY),
				    &(view->up.vtZ)) < 1;
	  break;
	case 'o':
	  outFileStr = optarg;
	  break;
	case 'r':
	  inRefFileStr = optarg;
	  break;
	case 'h': /* FALLTHROUGH */
	default:
	  usage = 1;
	  break;
      }
    }
    if(usage == 0)
    {
      if(optind < argc)
      {
	if(optind + 1 == argc)
	{
	  inFileStr = *(argv + optind);
	}
	else
	{
	  usage = 1;
	}
      }
    }
    ok = usage == 0;
  }
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
                     *argv, inFileStr, errMsg);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    if((inRefFileStr == NULL) ||
       (*inRefFileStr == '\0') ||
       ((fP = (strcmp(inRefFileStr, "-")?
              fopen(inRefFileStr, "r"): stdin)) == NULL) ||
       ((refObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
                     *argv, inRefFileStr, errMsg);
    }
    if(fP && strcmp(inRefFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    errNum = WlzInit3DViewStruct(view, inObj);
    if(errNum == WLZ_ERR_NONE)
    {
      outObj = WlzAssignObject(
	       WlzGetSubSectionFromObject(inObj, refObj, view, interp,
					  NULL, &errNum), NULL);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to cut section object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object to file %s (%s).\n",
                     *argv, outFileStr, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFree3DViewStruct(view);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(refObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-a<pitch,yaw,roll>] [-f <fx,fy,fz>]\n"
    "                [-d <dist>] [-h] [-m <mode>] [-u<ux,uy,uz>]\n"
    "                [-o<out obj>] [-r<ref obj>] <in obj>\n"
    "Creates a new 2D object which is a section cut from the input object\n"
    "within the domain of the reference object.\n"
    "Options are:\n"
    "  -a  Viewing angles: pitch (phi), yaw (theta) and roll (beta),\n"
    "      default 0.0,0.0,0.0 (all angles in degrees).\n"
    "  -f  Fixed point position, default 0.0,0.0,0.0.\n"
    "  -d  Distance parameter, default 0.0.\n"
    "  -h  Help, prints this usage message.\n"
    "  -m  Viewing mode, one of: up-is-up, statue or absolute, default\n"
    "      is up-is-up.\n"
    "  -u  Up vector, default 0.0,0.0,1.0.\n"
    "  -o  Output object file name.\n"
    "  -r  Reference object, within which cut section lies.\n",
    *argv,
    " -o out.wlz -r ref.wlz in.wlz\n"
    "Creates output object out.wlz by cutting a section through in.wlz\n"
    "from the plane z=0 of the input object in.wlz, such that it is\n"
    "contained within the domain of the object ref.wlz\n");
     

  }
  return(!ok);
}

/*!
* \return       Number of fields found, which will be -ve if there is a
*               parsing error.
* \ingroup      BinWlzApp
* \brief        Parses the input string for three coma seperated double
*               values, any of which may be missing, eg "4, ,2.367 " and
*               ",1.0 , " are both valid (white space is ignored).
* \param        str                     String to parse for double values.
* \param        d0                      Destination pointer for first value,
*                                       must not be NULL.
* \param        d1                      Destination pointer for second value,
*                                       must not be NULL.
* \param        d2                      Destination pointer for third value,
*                                       must not be NULL.
*/
static int      WlzGSFOReadTriple(char *str,
                                  double *d0, double *d1, double *d2)
{
  int           idx,
                cnt = 0;
  char          *subStr[3];
  double        *dP[3];

  if(str)
  {
    dP[0] = d0; dP[1] = d1; dP[2] = d2;
    while(*str && isspace(*str))
    {
      ++str;
    }
    if(*str == ',')
    {
      subStr[0] = NULL;
      subStr[1] = strtok(str, ",");
      subStr[2] = strtok(NULL, ",");
    }
    else
    {
      subStr[0] = strtok(str, ",");
      subStr[1] = strtok(NULL, ",");
      subStr[2] = strtok(NULL, ",");
    }
    if((subStr[0] == NULL) && (subStr[1] == NULL) && (subStr[2] == NULL))
    {
      cnt = -1;
    }
    else
    {
      for(idx = 0; idx < 3; ++idx)
      {
        if(subStr[idx] && (sscanf(subStr[idx], "%lg", dP[idx]) == 1))
        {
          ++cnt;
        }
        else
        {
          cnt = -1;
          break;
        }
      }
    }
  }
  return(cnt);
}

