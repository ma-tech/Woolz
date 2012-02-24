#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRandomAffineTransform_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzRandomAffineTransform.c
* \author       Bill Hill
* \date         August 2003
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
* \brief	Creates an affine transform with random transfrom primatives.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzrandomaffinetransform "WlzRandomAffineTransform"
*/

/*!
\ingroup BinWlz
\defgroup wlzrandomaffinetransform WlzRandomAffineTransform
\par Name
WlzRandomAffineTransform  -  creates an affine transform with random
                             transfrom primatives.
\par Synopsis
\verbatim
WlzRandomAffineTransform [-h] [-o<output object>] [-2] [-3]
                         [-R] [-g] [-r] [-t]
			 [-x #] [-X #] [-y #] [-Y #] [-z #] [-Z #]
			 [-a #] [-A #] [-b #] [-B #]
			 [-u #] [-U #] [-v #] [-V #] [-w #] [-W #]
			 [-s #] [-S #] [-e #]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output affine transform object.</td>
  </tr>
  <tr> 
    <td><b>-2</b></td>
    <td>2D transform, default.</td>
  </tr>
  <tr> 
    <td><b>-3</b></td>
    <td>3D transform.</td>
  </tr>
  <tr> 
    <td><b>-R</b></td>
    <td>Use radians for angles instead of degrees.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>General affine transform.</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Rigid body transform.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Translation only transform.</td>
  </tr>
  <tr> 
    <td><b>-x</b></td>
    <td>Mean column (x) translation.</td>
  </tr>
  <tr> 
    <td><b>-X</b></td>
    <td>Half range of column (x) translation.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Mean row (y) translation.</td>
  </tr>
  <tr> 
    <td><b>-Y</b></td>
    <td>Half range of row (y) translation.</td>
  </tr>
  <tr> 
    <td><b>-z</b></td>
    <td>Mean plane (z) translation.</td>
  </tr>
  <tr> 
    <td><b>-Z</b></td>
    <td>Half range of plane (z) translation.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Mean rotation about the z-axis.</td>
  </tr>
  <tr> 
    <td><b>-A</b></td>
    <td>Half range of rotation about the z-axis.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Mean rotation about the y-axis.</td>
  </tr>
  <tr> 
    <td><b>-B</b></td>
    <td>Half range of rotation about the y-axis.</td>
  </tr>
  <tr> 
    <td><b>-u</b></td>
    <td>Mean shear strength.</td>
  </tr>
  <tr> 
    <td><b>-U</b></td>
    <td>Half range of shear strength.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Mean shear angle in x-y plane.</td>
  </tr>
  <tr> 
    <td><b>-V</b></td>
    <td>Half range of shear angle in x-y plane.</td>
  </tr>
  <tr> 
    <td><b>-w</b></td>
    <td>Mean 3D shear angle.</td>
  </tr>
  <tr> 
    <td><b>-W</b></td>
    <td>Half range of 3D shear angle.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Mean scale factor.</td>
  </tr>
  <tr> 
    <td><b>-S</b></td>
    <td>Half range of scale factor.</td>
  </tr>
  <tr> 
    <td><b>-e</b></td>
    <td>Integer seed for random number generator.</td>
  </tr>
</table>
\par Description
Creates a random affine transform from random transfrom primatives.
The default mean and half range of the primatives is 0.0 except for
the scale which is 1.0 \f$\pm\f$ 0.0.
\par Examples
\verbatim
WlzRandomAffineTransform -t -2 -x 1 -X 4.2 -y 1 -Y 8
\endverbatim
Creates a 2D translation only transform with random translations of
1.0 \f$\pm\f$ 4.2 (x) and 1.0 \f$\pm\f$ 8.0 (y)
which is written to the standard output.
\par File
\ref WlzRandomAffineTransform.c "WlzRandomAffineTransform.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzaffinetransformobj "WlzAffineTransformObj(1)"
\ref WlzAffineTransformFromPrimVal "WlzAffineTransformFromPrimVal(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/time.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		usage = 0,
		dim = 2,
		radiansF = 0,
		seedF = 0;
  long		seed;
  WlzTransformType trType = WLZ_TRANSFORM_2D_AFFINE;
  WlzAffineTransform *trans = NULL;
  WlzValues     trVal;
  WlzDomain     trDom;
  WlzObject	*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outObjFileStr;
  WlzAffineTransformPrim prim0,
			 primR;
  struct timeval tp;
  static char	optList[] = "ho:23Rgrtx:X:y:Y:z:Z:a:A:b:B:u:U:v:V:w:W:s:S:e:",
		outObjFileStrDef[] = "-";
  
  opterr = 0;
  prim0.tx = 0.0;
  prim0.ty = 0.0;
  prim0.tz = 0.0;
  prim0.scale = 1.0;
  prim0.theta = 0.0;
  prim0.phi = 0.0;
  prim0.alpha = 0.0;
  prim0.psi = 0.0;
  prim0.xsi = 0.0;
  prim0.invert = 0;
  primR = prim0;
  primR.scale = 0.0;
  outObjFileStr = outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case '2':
	dim = 2;
        break;
      case '3':
	dim = 3;
        break;
      case 'g':
        trType = WLZ_TRANSFORM_2D_AFFINE;
	break;
      case 'r':
        trType = WLZ_TRANSFORM_2D_REG;
	break;
      case 't':
        trType = WLZ_TRANSFORM_2D_TRANS;
	break;
      case 'R':
        radiansF = 1;
	break;
      case 'x':
	if(sscanf(optarg, "%lg", &(prim0.tx)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'X':
	if(sscanf(optarg, "%lg", &(primR.tx)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'y':
	if(sscanf(optarg, "%lg", &(prim0.ty)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'Y':
	if(sscanf(optarg, "%lg", &(primR.ty)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'z':
	if(sscanf(optarg, "%lg", &(prim0.tz)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'Z':
	if(sscanf(optarg, "%lg", &(primR.tz)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'a':
	if(sscanf(optarg, "%lg", &(prim0.theta)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'A':
	if(sscanf(optarg, "%lg", &(primR.theta)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'b':
	if(sscanf(optarg, "%lg", &(prim0.phi)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'B':
	if(sscanf(optarg, "%lg", &(primR.phi)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'u':
	if(sscanf(optarg, "%lg", &(prim0.alpha)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'U':
	if(sscanf(optarg, "%lg", &(primR.alpha)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'v':
	if(sscanf(optarg, "%lg", &(prim0.psi)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'V':
	if(sscanf(optarg, "%lg", &(primR.psi)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'w':
	if(sscanf(optarg, "%lg", &(prim0.xsi)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'W':
	if(sscanf(optarg, "%lg", &(primR.xsi)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 's':
	if(sscanf(optarg, "%lg", &(prim0.scale)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'S':
	if(sscanf(optarg, "%lg", &(primR.scale)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'e':
	seedF = 1;
	if(sscanf(optarg, "%ld", &seed) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok && (radiansF == 0))
  {
    prim0.theta *= WLZ_M_PI / 180.0;
    prim0.phi *= WLZ_M_PI / 180.0;
    prim0.psi *= WLZ_M_PI / 180.0;
    prim0.xsi *= WLZ_M_PI / 180.0;
    primR.theta *= WLZ_M_PI / 180.0;
    primR.phi *= WLZ_M_PI / 180.0;
    primR.psi *= WLZ_M_PI / 180.0;
    primR.xsi *= WLZ_M_PI / 180.0;
  }
  if(ok && (dim == 3))
  {
    switch(trType)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
        trType = WLZ_TRANSFORM_3D_AFFINE;
	break;
      case WLZ_TRANSFORM_2D_REG:
        trType = WLZ_TRANSFORM_3D_REG;
	break;
      case WLZ_TRANSFORM_2D_TRANS:
        trType = WLZ_TRANSFORM_3D_TRANS;
	break;
      default:
        break;
    }
  }
  /* Compute the random primatives. */
  if(ok)
  {
    if(seedF == 0)
    {
#if defined(_WIN32)
      (void )time(&tp);
#else
      (void )gettimeofday(&tp, NULL);
#endif
      seed = (tp.tv_sec * 1000000) + tp.tv_usec;
    }
    AlgRandSeed(seed);
    if(fabs(primR.tx) > 0.0)
    {
      prim0.tx += ((2.0 * AlgRandUniform()) - 1.0) * primR.tx;
    }
    if(fabs(primR.ty) > 0.0)
    {
      prim0.ty += ((2.0 * AlgRandUniform()) - 1.0) * primR.ty;
    }
    if(fabs(primR.tz) > 0.0)
    {
      prim0.tz += ((2.0 * AlgRandUniform()) - 1.0) * primR.tz;
    }
    if(fabs(primR.scale) > 0.0)
    {
      prim0.scale += ((2.0 * AlgRandUniform()) - 1.0) * primR.scale;
    }
    if(fabs(primR.theta) > 0.0)
    {
      prim0.theta += ((2.0 * AlgRandUniform()) - 1.0) * primR.theta;
    }
    if(fabs(primR.phi) > 0.0)
    {
      prim0.phi += ((2.0 * AlgRandUniform()) - 1.0) * primR.phi;
    }
    if(fabs(primR.alpha) > 0.0)
    {
      prim0.alpha += ((2.0 * AlgRandUniform()) - 1.0) * primR.alpha;
    }
    if(fabs(primR.psi) > 0.0)
    {
      prim0.psi += ((2.0 * AlgRandUniform()) - 1.0) * primR.psi;
    }
    if(fabs(primR.xsi) > 0.0)
    {
      prim0.xsi += ((2.0 * AlgRandUniform()) - 1.0) * primR.xsi;
    }
    if((trans = WlzAffineTransformFromPrimVal(trType,
    				prim0.tx, prim0.ty, prim0.tz,
				prim0.scale, prim0.theta, prim0.phi,
				prim0.alpha, prim0.psi, prim0.xsi,
				prim0.invert, &errNum)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		  "%s: failed to make an affine transform from primitives\n",
		     *argv);
    }
  }
  if(ok)
  {
    trDom.t = trans;
    trVal.core = NULL;
    outObj = WlzMakeMain(WLZ_AFFINE_TRANS, trDom, trVal, NULL, NULL,
			 &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to make transform object.\n",
		     *argv);
    }
    if(ok)
    {
      if(((fP = (strcmp(outObjFileStr, "-")?
		fopen(outObjFileStr, "w"):
		stdout)) == NULL) ||
	  (WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to write output transform object\n",
		       *argv);
      }
      if(fP && strcmp(outObjFileStr, "-"))
      {
	fclose(fP);
      }
    }
  }
  if(outObj)
  {
    WlzFreeObj(outObj);
  }
  else if(trans)
  {
    WlzFreeAffineTransform(trans);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-o<output object>] [-2] [-3]\n"
    "        [-R] [-g] [-r] [-t]\n"
    "        [-x #] [-X #] [-y #] [-Y #] [-z #] [-Z #]\n"
    "        [-a #] [-A #] [-b #] [-B #]\n"
    "        [-u #] [-U #] [-v #] [-V #] [-w #] [-W #]\n"
    "        [-s #] [-S #]\n"
    "        [-e #]\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output affine transform object.\n" 
    "  -2  2D transform (default).\n"
    "  -3  3D transform instead of 2D.\n"
    "  -R  Use radians for angles instead of degrees.\n"
    "  -g  General affine transform.\n"
    "  -r  Rigid body transform.\n"
    "  -t  Translation only transform.\n"
    "  -x  Mean column (x) translation.\n"
    "  -X  Half range of column (x) translation.\n"
    "  -y  Mean row (y) translation.\n"
    "  -Y  Half range of row (y) translation.\n"
    "  -z  Mean plane (z) translation.\n"
    "  -Z  Half range of plane (z) translation.\n"
    "  -a  Mean rotation about the z-axis.\n"
    "  -A  Half range of rotation about the z-axis.\n"
    "  -b  Mean rotation about the y-axis.\n"
    "  -B  Half range of rotation about the y-axis.\n"
    "  -u  Mean shear strength.\n"
    "  -U  Half range of shear strength.\n"
    "  -v  Mean shear angle in x-y plane.\n"
    "  -V  Half range of shear angle in x-y plane.\n"
    "  -w  Mean 3D shear angle.\n"
    "  -W  Half range of 3D shear angle.\n"
    "  -s  Mean scale factor.\n" 
    "  -S  Half range of scale factor.\n" 
    "  -e  Integer seed for random number generator.\n" 
    "Creates a random affine transform from random transfrom primatives.\n"
    "The default mean and half range of the primatives is 0.0 except for\n"
    "the scale which is 1.0 +/- 0.0.\n",
    *argv,
    " -t -2 -x 1 -X 4.2 -y 1 -Y 8\n"
    "Creates a 2D translation only transform with random translations of\n"
    " 1.0 +/- 4.2 (x) and 1.0 +/- 8.0 (y) which is written to the standard\n"
    "output.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
