#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDefGrdTensorFeatures_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzDefGrdTensorFeatures.c
* \author       Bill Hill
* \date         January 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Functions to extract features from deformation gradient
* 		tensor fields.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzdefgrdtensorfeatures "WlzDefGrdTensorFeatures"
*/

/*!
\ingroup BinWlz
\defgroup wlzdefgrdtensorfeatures WlzDefGrdTensorFeatures
\par Name
WlzDefGrdTensorFeatures - computes features of deformation gradient
			 tensor fields.
\par Synopsis
\verbatim
WlzDefGrdTensorFeatures [-d <x>,<y>,<z>] [-f <feature list>]
			[-g <x>,<y>,<z>] [-h] [-o<out file>]
			[-p] [-s <dst>] [-u] [-x] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-d</b></td>
    <td>Dither range for point positions (default no dithering,
        dithering implies a points output object).</td>
  </tr> <tr> 
    <td><b>-g</b></td>
    <td>Smooth values using a Gaussian with the given sigma values
        (default no smoothing).</td>
  </tr> <tr> 
    <td><b>-f</b></td>
    <td>Features for extraction identified by a comma separated list of
        single characters.</td>
  </tr> <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr> <tr> 
    <td><b>-p</b></td>
    <td>Output points objects rather than vector fields.</td>
  </tr> <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr> <tr> 
    <td><b>-s</b></td>
    <td>Minimum point seperation. If generating points objects then\n"
       this specifies the minimum distance between points.</td>
  </tr> <tr> 
    <td><b>-u</b></td>
    <td>If a single feature is requested do not use a compound object.</td>
  </tr> <tr>    
    <td><b>-x</b></td>
    <td>Use voxel size scaling for points objects.</td>
  </tr>
</table>
\par Description
Computes features of the input object which must be a deformation
gradient tensor field. The output is a compound object with a
component object per feature, with the features identifiable through
their name properties or ordering. The features for extraction may
be defined by a comma separated list of single characters. Available
features are:
<table border="0">
  <tr>
    <td>identifying character</td><td>property name</td><td>feature</td>
  </tr> <tr>
    <td>d</td><td>direction vectors</td><td>principle direction vectors</td>
  </tr> <tr>
    <td>j</td><td>jacobian</td><td>determinant of Jacobian tensor</td>
  </tr> <tr>
    <td>s</td><td>stretch values</td><td>principle stretch values</td>
  </tr>
Default feature selection is d,j,s.
</table>

By default all files are read from the standard input and are written
to the standard output.
\par Examples
\verbatim
WlzDefGrdTensorFeatures -f d -p -s 4 -d 20,20,20
                             -g 2,2,2 -o features.wlz dgf.wlz
\endverbatim
Creates a compound object with a single component points object in
which the point domain locations are scattered over the given field's
domain with a separation of 4 voxels but then dithered by an additional
range of +/-20 voxels. The principle directions are computed and smoothed
using a Gaussian with sigma values of 2,2,2 prior to sampling at the
point locations.  The associated point values will be vector values with
a vector for each principle direction.
The points object will also have it's name property set to "direction vectors".
The compound object will be written to the file features.wlz.
\par File
\ref WlzDefGrdTensorFeatures.c "WlzDefGrdTensorFeatures.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzDGTensorFeatures "WlzDGTensorFeatures(3)"
\ref WlzCMeshStrainTensor "WlzCMeshStrainTensor(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>


/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		option,
  		ok = 1,
		points = 0,
		singleObj = 0,
  		usage = 0,
		voxScaling = 0;
  int		nFeat = 3;
  unsigned int	features = 
                    WLZ_DGTENSOR_FEATURE_MASK(WLZ_DGTENSOR_FEATURE_DETJAC) |
  		    WLZ_DGTENSOR_FEATURE_MASK(WLZ_DGTENSOR_FEATURE_DIRVEC) |
		    WLZ_DGTENSOR_FEATURE_MASK(WLZ_DGTENSOR_FEATURE_STRVAL);
  double	dMin = 0.0;
  FILE		*fP = NULL;
  char		*inFileStr = NULL,
  		*outFileStr = NULL;
  WlzDVertex3   dither = {0.0},
		smooth = {0.0};
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char   optList[] = "hpuxd:f:g:o:s:";
  const char    fileStrDef[] = "-";

  opterr = 0;
  inFileStr = (char *)fileStrDef;
  outFileStr = (char *)fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
	points = 1;
        if(sscanf(optarg,
                  "%lg,%lg,%lg",
                  &(dither.vtX),
                  &(dither.vtY),
                  &(dither.vtZ)) < 3)
        {
          usage = 1;
        }
	break;
      case 'f':
	nFeat = 0;
	features = WLZ_DGTENSOR_FEATURE_NONE;
	{
	  char	*fStr;
	  
	  fStr = strtok(optarg, ",");
	  while(!usage && fStr)
	  {
	    if(strlen(fStr) != 1)
	    {
	      usage = 1;
	    }
	    else
	    {
	      ++nFeat;
	      switch(*fStr)
	      {
		case 'd':
		  features |= WLZ_DGTENSOR_FEATURE_MASK(
		  	      WLZ_DGTENSOR_FEATURE_DIRVEC);
		  break;
		case 'j':
		  features |= WLZ_DGTENSOR_FEATURE_MASK(
		  	      WLZ_DGTENSOR_FEATURE_DETJAC);
		  break;
		case 's':
		  features |= WLZ_DGTENSOR_FEATURE_MASK(
		  	      WLZ_DGTENSOR_FEATURE_STRVAL);
		  break;
		default:
		  usage = 1;
		  break;
	      }
	      fStr = strtok(NULL, ",");
	    }
	  }
	}
	break;
      case 'g':
        if(sscanf(optarg,
                  "%lg,%lg,%lg",
                  &(smooth.vtX),
                  &(smooth.vtY),
                  &(smooth.vtZ)) < 3)
        {
          usage = 1;
        }
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'p':
        points = 1;
        break;
      case 's':
        points = 1;
        if(sscanf(optarg, "%lg", &(dMin)) != 1)
        {
          usage = 1;
        }
        break;
      case 'u':
        singleObj = 1;
	break;
      case 'x':
        voxScaling = 1;
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
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inFileStr = *(argv + optind);
      }
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if(((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    outObj = WlzDGTensorFeatures(inObj, features, points, dMin, dither,
                                 smooth, voxScaling, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		"%s: Failed to compute tensor features object (%s)\n",
		argv[0],
		errMsgStr);
    }
  }
  WlzFreeObj(inObj);
  if(ok && (nFeat == 1) && singleObj)
  {
    WlzObject	*tObj;
    WlzCompoundArray *cObj;

    cObj = (WlzCompoundArray *)outObj;
    if((tObj = cObj->o[0]) != NULL)
    {
      cObj->o[0] = NULL;
      tObj->linkcount = 0;
    }
    (void )WlzFreeObj(outObj);
    outObj = tObj;
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to write object to file %s\n",
                     *argv, outFileStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
      "Usage: %s [-d <x>,<y>,<z>] [-f <feature list>]\n"
      "\t\t[-g <x>,<y>,<z] \t\t[-h] [-o<out file>] [-p]\n"
      "\t\t[-s <dst>] [-u] [-x] [<input file>]\n"
      "Computes features of the input object which must be a deformation\n"
      "gradient tensor field. The output is a compound object with a\n"
      "component object per feature, with the features identifiable through\n"
      "their name properties or ordering. The features for extraction may\n"
      "be defined by a comma separated list of single characters. Available\n"
      "features are:\n"
      "  *id char* *property name*     *feature*\n"
      "   d         direction vectors   principle direction vectors\n"
      "   j         jacobian            determinant of Jacobian tensor\n"
      "   s         stretch values      principle principle stretch values\n"
      "Default feature selection is d,j,s.\n"
      "By default all files are read from the standard input and are written\n"
      "to the standard output.\n"
      "Version: %s\n"
      "Options are:\n"
      "  -d  Dither range for point positions (default no dithering,\n"
      "      dithering implies a points output object).\n"
      "  -f  Features for extraction identified by a comma separated list of\n"
      "      single characters.\n"
      "  -g  Smooth values using a Gaussian with the given sigma values\n"
      "      (default no smoothing).\n"
      "  -h  Help, prints this usage message.\n"
      "  -o  Output object.\n"
      "  -p  Make the output objects points objects rather than vector\n"
      "      field objects.\n"
      "  -s  Minimum point seperation. If generating points objects then\n"
      "      this specifies the minimum distance between points.\n"
      "  -u  If a single feature is requested do not use a compound object.\n"
      "  -x  Use voxel size scaling for points objects.\n"
      "Example:\n"
      "%s -f d -p -s 4 -d 20,20,20\n"
      "\t\t-g 2,2,2 -o features.wlz dgf.wlz\n"
      "Creates a compound object with a single component points object\n"
      "in which the point domain locations are scattered over the given\n"
      "field's domain with a separation of 4 voxels but then dithered by\n"
      "an additional range of +/-20 voxels. The principle directions are\n"
      "computed and smoothed using a Gaussian with sigma values of 2,2,2\n"
      "prior to sampling at the point locations. The associated point\n"
      "values will be vector values with a vector for each principle\n"
      "direction. The points object will also have it's name property set\n"
      "to 'direction vectors'. The compound object will be written to the\n"
      "file features.wlz.\n",
      argv[0],
      WlzVersion(),
      argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
