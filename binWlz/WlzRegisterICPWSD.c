#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzRegisterICPWSD.c
* \author       Bill Hill
* \date         May 2004
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
* \brief	Computes a 3D grey valued object in which the values are
*		the weighted sums of distances as used for ICP registration.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzregistericpwsd "WlzRegisterICPWSD"
*/

/*!
\ingroup BinWlz
\defgroup wlzregistericpwsd WlzRegisterICPWSD
\par Name
WlzRegisterICPWSD - computes a 3D grey valued object using ICP weights.
\par Synopsis
\verbatim
WlzRegisterICPWSD - [-M #] [-R] [-i <init tr>] [-o<out obj>]
                    [-r #,#:#] [-x #,#:#] [-y #,#:#]
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
    <td>Output file name.</td>
  </tr>
  <tr> 
    <td><b>-M</b></td>
    <td>Minimum distance weight, range [0.0-1.0]: Useful values
        are 0.25 (default) for global matching and 0.0 for local
	matching.</td>
  </tr>
  <tr> 
    <td><b>-R</b></td>
    <td>Rotations are in radians not degrees.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Initial affine transform object.</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Range and step size for rotation in degrees.</td>
  </tr>
  <tr> 
    <td><b>-x</b></td>
    <td>Range and step size for translation through columns.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Range and step size for translation through lines.</td>
  </tr>
</table>
\par Description
Computes a 3D double precission valued object in which the values
are the weighted sums of distances as used for ICP registration.
In specifying the ranges the defaults are 0.0 -> 20.0, step 1.0.
The input objects are read from stdin and values are written to stdout
unless the filenames are given.
\par Examples
\verbatim
lzRegisterICPWSD -o out.wlz -r 0,90:0.5 -x 100,200: -y ,200: in0.wlz in1.wlz
\endverbatim
A weighted sums of distances object is computed for rotations of
0 to 90 degrees with a step size of 0.5 degrees, the x and y
translations are from 100 to 200 and 0 to 200 in unit steps.
The weighted sums of distances object is then written to out.wlz.
\par File
\ref WlzRegisterICPWSD.c "WlzRegisterICPWSD.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzregistericp "WlzRegisterICP(1)"
\ref WlzRegICPObjWSD2D "WlzRegICPObjWSD2D(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		idN,
		option,
		useRadians = 0,
		ok = 1,
		usage = 0;
  double	minDistWgt = 0.25,
  		rMin = 0.0,
  		rMax = 20.0,
		rStep = 1.0,
		xMin = 0.0,
		xMax = 20.0,
		xStep = 1.0,
		yMin = 0.0,
		yMax = 20.0,
		yStep = 1.0;
  char		*comma,
  		*colon;
  double	*range[3];
  char		*buf[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inTrObj = NULL,
  		*outObj = NULL;
  WlzObject	*inObj[2];
  FILE		*fP = NULL;
  char 		*inTrObjFileStr = NULL,
  		*outObjFileStr;
  char  	*inObjFileStr[2];
  const char	*errMsg;
  static char	optList[] = "M:i:o:r:x:y:Rh",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObj[0] = NULL;
  inObj[1] = NULL;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr[0] = inObjFileStrDef;
  inObjFileStr[1] = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'R':
        useRadians = 1;
	break;
      case 'M':
        if((sscanf(optarg, "%lg", &minDistWgt) != 1) ||
	   (minDistWgt < 0.0) || (minDistWgt > 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'i':
        inTrObjFileStr = optarg;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'r': /* FALLTHROUGH */
      case 'x': /* FALLTHROUGH */
      case 'y':
	if(optarg)
	{
	  comma = strchr(optarg, ',');
	  colon = strchr(optarg, ':');
	  if(comma && colon)
	  {
	    *comma = *colon = '\0';
	    switch(option)
	    {
	      case 'r':
		range[0] = &rMin; 
		range[1] = &rMax; 
		range[2] = &rStep; 
		break;
	      case 'x':
		range[0] = &xMin; 
		range[1] = &xMax; 
		range[2] = &xStep; 
		break;
	      case 'y':
		range[0] = &yMin; 
		range[1] = &yMax; 
		range[2] = &yStep; 
		break;
	    }
	    buf[0] = optarg;
	    buf[1] = comma + 1;
	    buf[2] = colon + 1;
	    for(idN = 0; (idN < 3) && (usage == 0); ++idN)
	    {
	      while(*buf[idN] && isspace(*buf[idN]))
	      {
	        ++buf[idN];
	      }
	      if(*buf[idN])
	      {
	        if(sscanf(buf[idN], "%lg", range[idN]) != 1)
		{
		  usage = 1;
		}
	      }
	    }
	  }
	  else
	  {
	    usage = 1;
	  }
	}
        break;
      case 'h':
      default:
        usage = 1;
	break;
    }
  }
  if(useRadians == 0)
  {
    rMin *= ALG_M_PI / 180.0;
    rMax *= ALG_M_PI / 180.0;
    rStep *= ALG_M_PI / 180.0;
  }
  if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
     (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
  {
    usage = 1;
  }
  ok = !usage;
  if(ok && (optind < argc))
  {
    idN = 0;
    while((idN < 2) && (optind < argc))
    {
      inObjFileStr[idN] = *(argv + optind);
      ++optind;
      ++idN;
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
	     "%s: failed to read initial affine transform from file %s (%s)\n",
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
      		     "%s: initial affine transform object invalid type\n",
		     *argv);
    }
  }
  if(ok)
  {
    /* Read objects. */
    idN = 0;
    while((errNum == WLZ_ERR_NONE) && (idN < 2))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((inObjFileStr[idN] == NULL) ||
	  (*inObjFileStr[idN] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idN], "-")?
		  fopen(inObjFileStr[idN], "r"): stdin)) == NULL) ||
	  ((inObj[idN] = WlzAssignObject(WlzReadObj(fP,
	  					    &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to read object %d from file %s (%s)\n",
		       *argv, idN, inObjFileStr[idN], errMsg);
      }
      if(fP && strcmp(inObjFileStr[idN], "-"))
      {
	fclose(fP);
      }
      ++idN;
    }
  }
  if(ok)
  {
    /* Check object types. */
    if(inObj[0]->type != inObj[1]->type)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if((inObj[0]->domain.core == NULL) ||
            (inObj[1]->domain.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      switch(inObj[0]->type)
      {
        case WLZ_2D_DOMAINOBJ:
          break;
        case WLZ_CONTOUR:
          if((inObj[0]->domain.core == NULL) ||
             (inObj[0]->domain.ctr->model == NULL))
          {
            errNum = WLZ_ERR_DOMAIN_NULL;
          }
          else
          {
            switch(inObj[0]->domain.ctr->model->type)
            {
              case WLZ_GMMOD_2I: /* FALLTHROUGH */
              case WLZ_GMMOD_2D:
              case WLZ_GMMOD_2N:
                break;
              default:
                errNum = WLZ_ERR_DOMAIN_TYPE;
                ok = 0;
                break;
            }
          }
          break;
        default:
          errNum = WLZ_ERR_OBJECT_TYPE;
          ok = 0;
          break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: input object(s) not appropriate\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    outObj = WlzRegICPObjWSD2D(inObj[0], inObj[1],
    			       inTrObj? inTrObj->domain.t: NULL,
			       xMin, xMax, xStep,
			       yMin, yMax, yStep,
			       rMin, rMax, rStep,
			       minDistWgt, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		   "%s: failed to compute weighted sums of distances (%s).\n",
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
		     "%s: failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
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
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-M #] [-R] [-i <init tr>] [-o<out obj>]\n"
    "                      [-r #,#:#] [-x #,#:#] [-y #,#:#]\n"
    "                      [<in obj 0>] [<in obj 1>]\n"
    "Options:\n"
    "  -M  Minimum distance weight, range [0.0-1.0]: Useful values are\n"
    "      0.25 (default) for global matching and 0.0 for local matching.\n"
    "  -R  Rotations are in radians not degrees.\n"
    "  -i  Initial affine transform object.\n"
    "  -o  Output file name.\n"
    "  -r  Range and step size for rotation in degrees.\n"
    "  -x  Range and step size for translation through columns.\n"
    "  -y  Range and step size for translation through lines.\n"
    "  -h  Help, prints this usage message.\n"
    "Computes a 3D double precission valued object in which the values\n"
    "are the weighted sums of distances as used for ICP registration.\n"
    "In specifying the ranges the defaults are 0.0 -> 20.0, step 1.0.\n"
    "The input objects are read from stdin and values are written to stdout\n"
    "unless the filenames are given.\n",
    *argv,
    " -o out.wlz -r 0,90:0.5 -x 100,200: -y ,200: in0.wlz in1.wlz\n"
    "A weighted sums of distances object is computed for rotations of\n"
    "0 to 90 degrees with a step size of 0.5 degrees, the x and y\n"
    "translations are from 100 to 200 and 0 to 200 in unit steps.\n"
    "The weighted sums of distances object is then written to out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
