#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRmDirBleed_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzRmDirBleed.c
* \author       Bill Hill
* \date         November 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Removes directional bleed from a Woolz image object.
* \ingroup	wlzrmdirbleed "WlzRmDirBleed"
*/

/*!
\ingroup BinWlz
\defgroup wlzrmdirbleed WlzRmDirBleed
\par Name
WlzRmDirBleed - removes directional bleed from a Woolz image object
\par Synopsis
\verbatim
WlzRmDirBleed [-h] [-l] [-m] [-R] [-T] [-a#] [-b#] [-g#]
              [-o<out object file>] [<in object file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-a</b></td>
    <td>Alpha, fraction of buffer to remove from current plane.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Beta, fraction of previous plane add to buffer.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Gamma parameter, spread of bleed per section.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr>
    <td><b>-l</b></td>
    <td>Preserve plane image value limits.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Histogram match planes so that the input and output histograms are
        similar.</td>
  </tr>
  <tr>
    <td><b>-R</b></td>
    <td>Reverse direction (last plane to first).</td>
  </tr>
  <tr>
    <td><b>-T</b></td>
    <td>Make the output object a tiled object.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Attempts to remove directional bleed from a Woolz image object.
By default all files are read from the standard input and written
to the standard output.
 A single plane buffer (\f$f\f$) is updated at each plane using the contents
of the current and previous source planes (\f$s_p\f$ and \f$s_{p - 1}\f$):
\f[
  f = (1/(1 + \beta))  gauss(s_p + \beta  s_{p - 1}, \gamma)
\f]
At each destination plane (d_p) a fraction of the buffer is then
subtracted:
\f[
  d_p = s_p - \alpha f
\f]
By default all files are read from the standard input and written to
the standard output.
\par Examples
\verbatim
WlzRmDirBleed -oout.wlz -a 0.5 -b 0.6 -g 2.0 -m in.wlz
\endverbatim
a new image object with the same domain and grey type as the
input object (in.wlz) in which directional bleed has been removed
is written to out.wlz.
\par File
\ref WlzRmDirBleed.c "WlzRmDirBleed.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <Wlz.h>

static WlzErrorNum		WlzRmDirBleed(
				  WlzObject *dObj,
				  WlzObject *sObj,
				  double alpha,
				  double gamma,
				  int nrm,
				  WlzRasterDir dir);
static WlzObject		*WlzGetXYPlane(
				  WlzObject *obj3,
				  int p,
				  int ngv,
				  int *dstCpy,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzSetXYPlane(
				  WlzObject *obj3,
				  WlzObject *obj2,
				  int p);

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		nrm = 0,
  		ok = 1,
  		option,
		usage = 0,
		tiledOut = 0;
  char		*inFileStr,
		*outFileStr;
  double	alpha = 0.5,
  		gamma = 2.0;
  FILE		*fP = NULL;
  WlzRasterDir dir = WLZ_RASTERDIR_DPILIC;
  WlzGreyType 	gType = WLZ_GREY_ERROR;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hlmRTa:g:o:";
  const char    fileStrDef[] = "-";

  /* Parse the argument list and check for input files. */
  opterr = 0;
  outFileStr = (char *)fileStrDef;
  inFileStr = (char *)fileStrDef;
  while(ok && !usage && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'a':
	if((sscanf(optarg, "%lg", &alpha) != 1) || (alpha < 0.0))
	{
	  usage = 1;
	}
        break;
      case 'l':
        nrm = 2;
	break;
      case 'm':
        nrm = 1;
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'g':
	if((sscanf(optarg, "%lg", &gamma) != 1) || (gamma < 0.0))
	{
	  usage = 1;
	}
        break;
      case 'R':
        dir = WLZ_RASTERDIR_IPILIC;
	break;
      case 'T':
        tiledOut = 1;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if(ok)
  {
    ok = !usage;
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
      inFileStr = argv[optind];
    }
  }
  /* Read the input object. */
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read input object from file %s.\n",
                     argv[0], inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Check object type. */
  if(ok)
  {
    switch(inObj->type)
    {
      case WLZ_3D_DOMAINOBJ:
	if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(inObj->values.core == NULL)
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else
	{
          gType = WlzGreyTypeFromObj(inObj, &errNum);
	}
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
	  "%s: Input object must be a 3D spatial domain object with\n"
	  "values.\n",
	  argv[0]);
    }
  }
  if(ok)
  {
    if(tiledOut && (strcmp(outFileStr, "-") == 0))
    {
      ok = 0;
      (void )fprintf(stderr,
	  "%s: Can not write tiled value objects to the standard\n"
	  "output.\n",
	  argv[0]);
    }
  }
  /* Create output object. */
  if(ok)
  {
    WlzDomain dom;                     
    WlzValues inVal;

    dom = inObj->domain;
    inVal = inObj->values;
    if(tiledOut)
    {
      inVal = inObj->values;
      outObj = WlzAssignObject(
               WlzMakeTiledValuesFromObj(inObj, inVal.t->tileSz, 0,
	                                 gType, 0, NULL, inVal.t->bckgrnd,
					 &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	if((fP = fopen(outFileStr, "w")) == NULL)
	{
	  errNum = WLZ_ERR_WRITE_EOF;
	}
	else
	{
	  errNum = WlzWriteObj(fP, outObj);
	  (void )fclose(fP); fP = NULL;
	}
      }
      (void )WlzFreeObj(outObj);
      outObj = NULL;
      if(errNum == WLZ_ERR_NONE)
      {
        if((fP = fopen(outFileStr, "r+")) == NULL)
	{
	  errNum = WLZ_ERR_WRITE_EOF;
	}
	else
	{
	  outObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL);
	}
      }
    }
    else
    {
      WlzValues outVal;
      WlzPixelV	bgd;

      bgd = WlzGetBackground(inObj, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	outVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					 dom.p->plane1,
					 dom.p->lastpl,
					 bgd, NULL,
					 &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	int	p;
        WlzObjectType	vtt;

	vtt = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, gType, NULL);
        for(p = dom.p->plane1; p <= dom.p->lastpl; ++p)
	{
	  int	q;

	  q = p - dom.p->plane1;
	  if(dom.p->domains[q].core)
	  {
	    WlzObject	*inObj2 = NULL;
            WlzValues	nullValues;

	    nullValues.core = NULL;
	    inObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	                         dom.p->domains[q], nullValues,
				 NULL, NULL, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      WlzValues	outVal2;

	      outVal2.v = WlzNewValueTb(inObj2, vtt, bgd, &errNum);
	      outVal.vox->values[q] = WlzAssignValues(outVal2, NULL);
	    }
	    (void )WlzFreeObj(inObj2);
	    if(errNum != WLZ_ERR_NONE)
	    {
	      break;
	    }
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        outObj = WlzAssignObject(WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, outVal,
			         NULL, NULL, &errNum), NULL);
      }
    }
  }
  /* Attempt to remove bleed. */
  if(ok)
  {
    errNum = WlzRmDirBleed(outObj, inObj, alpha, gamma, nrm, dir);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to remove bleed from object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  /* Output the processed object. */
  if(ok)
  {
    if(!tiledOut)
    {
      if((fP = (strcmp(outFileStr, "-")?
	       fopen(outFileStr, "w"): stdout)) == NULL)
      {
	ok = 0;
	(void )fprintf(stderr,
	    "%s: Failed to open output file %s.\n",
	    argv[0], outFileStr);
      }
      if(ok)
      {
	errNum = WlzWriteObj(fP, outObj);
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	  (void )fprintf(stderr,
	      "%s: Failed to write output object, %s.\n",
	      argv[0], errMsgStr);
	}
      }
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-l] [-m] [-R] [-T] [-a#] [-b#] [-g#]\n"
    "\t\t[-o<output object>] [<input object>]\n"
    "Attempts to remove directional bleed from a Woolz image object.\n"
    "By default all files are read from the standard input and written\n"
    "to the standard output.\n"
    " A single plane buffer (f) is updated at each plane using the contents\n"
    "of the current and previous source planes (s_p and s_{p - 1}):\n"
    "  f = (1/(1 + beta))  gauss(s_p + beta  s_{p - 1}, gamma)\n"
    "At each destination plane (d_p) a fraction of the buffer is then\n"
    "subtracted:\n"
    "  d_p = s_p - \alpha f\n"
    "Version: %s\n"
    "Options:\n"
    "  -a  Alpha, fraction of buffer to remove from current plane.\n"
    "  -b  Beta, fraction of previous plane add to buffer.\n"
    "  -g  Gamma parameter, spread of bleed per section.\n"
    "  -o  Output object file.\n"
    "  -l  Preserve plane image value limits.\n"
    "  -m  Histogram match planes so that the input and output histograms\n"
    "      are similar.\n"
    "  -R  Reverse direction (last plane to first).\n"
    "  -T  Make the output object a tiled object.\n"
    "  -h  Help - prints this usage message\n"
    "By default all files are read from the standard input and written to\n"
    "the standard output.\n"
    "Example:\n"
    "  %s -oout.wlz -a 0.5 -b 0.6 -g 2.0 -m in.wlz\n"
    "a new image object with the same domain and grey type as the\n"
    "input object (in.wlz) in which directional bleed has been removed\n"
    "is written to out.wlz.\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}


/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Removes directional bleed from the source image object
* 		writting the result into the destination image object.
* 		The two objects must share the same domain and may also
* 		share the same values.
*
* 		At each plane the destination values are set using the
* 		values of the previous section:
* 		\f[
 		  d_p = (1/(1 - \alpha))(s_p  gauss(s_{p - 1}, \gamma))
 		\f]
* 		where \f$ d_p \f$ and \f$s_p\f$ are the destination and
* 		source images at plane \f$p\f$. Gauss is a Gaussian
* 		smoothing function \f$\alpha\f$ and \f$\gamma\f$  are
* 		given parameters.
* \param	dObj			Destination object.
* \param	sObj			Source object.
* \param	alpha			Fraction of the Gaussain blured
* 					previous plane to subtract from the
* 					current plane, range [0.0-1.0].
* \param	gamma			Gaussian smoothing parameter.
* \param	nrm			Normalisation:
* 					  0 - none,
* 					  1 - histogram match the destination
* 					      planes to the source planes.
*					  2 - use plane min/max values.
* \param	dir			Direction: WLZ_RASTERDIR_IPILIC or
* 					WLZ_RASTERDIR_DPILIC.
*/
static WlzErrorNum		WlzRmDirBleed(
				  WlzObject *dObj,
				  WlzObject *sObj,
				  double alpha,
				  double gamma,
				  int nrm,
				  WlzRasterDir dir)
{
  int		dIsTiled;
  WlzIBox3 	bb;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzObject	*bufObj = NULL;
  WlzGreyWSpace gWSp[4];
  WlzIntervalWSpace iWSp[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check given objects. */
  if((dObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((dObj->type != WLZ_3D_DOMAINOBJ) || (dObj->type != sObj->type))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((dObj->domain.core == NULL) || (sObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dObj->values.core == NULL) || (sObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((alpha < 0.0) || (alpha > 1.0) ||
          ((dir != WLZ_RASTERDIR_IPILIC) && (dir != WLZ_RASTERDIR_DPILIC)))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    gType = WlzGreyTypeFromObj(sObj, &errNum);
    dIsTiled = WlzGreyTableIsTiled(dObj->values.core->type);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzGreyType gTypeD;

      gTypeD = WlzGreyTypeFromObj(dObj, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        if(gType != gTypeD)
        {
	  errNum = WLZ_ERR_GREY_TYPE;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(gType)
      {
	case WLZ_GREY_INT:   /* FALLTHROUGH */
	case WLZ_GREY_SHORT: /* FALLTHROUGH */
	case WLZ_GREY_UBYTE:
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  /* Create a 2D buffer object with integer values and big enough for any
   * x-y plane. */
  if(errNum == WLZ_ERR_NONE)
  {
    bb = WlzBoundingBox3I(sObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		*inp = NULL;

    if(errNum == WLZ_ERR_NONE)
    {
      if((inp = (int *)AlcCalloc((bb.xMax - bb.xMin + 1) *
			         (bb.yMax - bb.yMin + 1),
				 sizeof(int))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        WlzPixelV	bgdV;

	bgdV.v.inv = 0;
	bgdV.type = WLZ_GREY_INT;
	bufObj = WlzAssignObject(
		 WlzMakeRect(bb.yMin, bb.yMax, bb.xMin, bb.xMax,
			     WLZ_GREY_INT, inp, bgdV, NULL, NULL,
			     &errNum), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	bufObj->values.v->freeptr = AlcFreeStackPush(
				        bufObj->values.v->freeptr,
					(void *)inp, NULL);
      }
      else
      {
	AlcFree(inp);
      }
    }
  }
  /* Work down through the planes x-y planes removing bleed. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		p,
		plane1,
		lastpl,
		planei;
    double	alpha1;
    WlzObject	*hObj[2] = {NULL};
    WlzPlaneDomain *pDom;

    alpha1 = 1.0 / (1.0 - alpha);
    pDom = sObj->domain.p;
    if(dir == WLZ_RASTERDIR_IPILIC)
    {
      planei = 1;
      plane1 = pDom->plane1;
      lastpl = pDom->lastpl;
    }
    else
    {
      planei = -1;
      plane1 = pDom->lastpl;
      lastpl = pDom->plane1;
    }
    for(p = plane1;
        ((planei > 0) && (p <= lastpl)) || ((planei < 0) && (p >= lastpl));
	p += planei)
    {
      int	skip = 0;
      WlzObject	*bObj2 = NULL,
      		*sObj2 = NULL,
      		*dObj2 = NULL,
		*gObj2 = NULL;
      WlzPixelV	dMin,
      		dMax,
		sMin,
      		sMax;

      /* Create 2D source and destination objects also intersection of
       * the buffer with the source and destination objects. */
      sObj2 = WlzGetXYPlane(sObj, p, 0, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	dObj2 =  WlzGetXYPlane(dObj, p, 1, NULL, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	bObj2 = WlzAssignObject(
		WlzMakeMain(WLZ_2D_DOMAINOBJ, sObj2->domain, bufObj->values,
			    NULL, NULL, &errNum), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	hObj[1] = WlzHistogramObj(sObj2, 0, 0.0, 1.0, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	skip |= (hObj[0] == NULL)? 2: 0; /* Something wrong with previous plane
					  * so all processing in this plane. */
	skip |= (hObj[1] == NULL)? 1: 0; /* Something wrong with current plane
					  * so no processing of this plane,
					  * just copy out. */
	if(skip == 0)
	{
	  double d;

	  d = WlzHistogramDistance(hObj[0], hObj[1], &errNum);
	  if(d < 0.7)
	  {
	    skip = 2;
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (skip == 2))
	{
	  errNum = WlzCopyObjectGreyValues(bObj2, sObj2);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  switch(nrm)
	  {
	    case 0:
	      break;
	    case 1:
	      break;
	    case 2:
	      errNum = WlzGreyRange(sObj2, &sMin, &sMax);
	      break;
	    default:
	      break;
	  }
	}
      }
      /* Compute gObj2 = gauss(bObj2, gamma) */
      if((errNum == WLZ_ERR_NONE) && (skip != 1))
      {
	gObj2 = WlzAssignObject(          
		WlzGauss2(bObj2, gamma, gamma, 0, 0, &errNum), NULL);
      }
      /* dObj2 = (1/(1 - alpha))gObj2
       * bObj2 = sObj2  ready for the next plane
       */
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzInitGreyScan(bObj2, &(iWSp[0]), &(gWSp[0]));
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzInitGreyScan(sObj2, &(iWSp[1]), &(gWSp[1]));
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzInitGreyScan(dObj2, &(iWSp[2]), &(gWSp[2]));
	}
	if((errNum == WLZ_ERR_NONE) && (skip != 1))
	{
	  errNum = WlzInitGreyScan(gObj2, &(iWSp[3]), &(gWSp[3]));
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  while((errNum == WLZ_ERR_NONE) &&
		((errNum = WlzNextGreyInterval(&(iWSp[0]))) == WLZ_ERR_NONE) &&
		((errNum = WlzNextGreyInterval(&(iWSp[1]))) == WLZ_ERR_NONE) &&
		((errNum = WlzNextGreyInterval(&(iWSp[2]))) == WLZ_ERR_NONE) &&
		((skip == 1) ||
		 ((errNum = WlzNextGreyInterval(&(iWSp[3]))) == WLZ_ERR_NONE)))
	  {
	    int	     i;
	    double   t;
	    int      *bp,
		     *gp;
	    WlzGreyP dp,
		     sp;

	    bp = gWSp[0].u_grintptr.inp;
	    sp = gWSp[1].u_grintptr;
	    dp = gWSp[2].u_grintptr;
	    if(skip == 1)
	    {
	      switch(gType)
	      {
		case WLZ_GREY_INT:
		  for(i = 0; i < iWSp[0].colrmn; ++i)
		  {
		    bp[i] = sp.inp[i];
		  }
		  break;
		case WLZ_GREY_SHORT:
		  for(i = 0; i < iWSp[0].colrmn; ++i)
		  {
		    bp[i] = sp.shp[i];
		  }
		  break;
		case WLZ_GREY_UBYTE:
		  for(i = 0; i < iWSp[0].colrmn; ++i)
		  {
		    bp[i] = sp.ubp[i];
		  }
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	    }
	    else
	    {
	      gp = gWSp[3].u_grintptr.inp;
	      switch(gType)
	      {
		case WLZ_GREY_INT:
		  for(i = 0; i < iWSp[0].colrmn; ++i)
		  {
		    t = alpha1 * (sp.inp[i] - (alpha * gp[i]));
		    dp.inp[i] = WLZ_CLAMP(t, INT_MIN, INT_MAX);
		    bp[i] = sp.inp[i];
		  }
		  break;
		case WLZ_GREY_SHORT:
		  for(i = 0; i < iWSp[0].colrmn; ++i)
		  {
		    t = alpha1 * (sp.shp[i] - (alpha * gp[i]));
		    dp.shp[i] = WLZ_CLAMP(t, SHRT_MIN, SHRT_MAX);
		    bp[i] = sp.shp[i];
		  }
		  break;
		case WLZ_GREY_UBYTE:
		  for(i = 0; i < iWSp[0].colrmn; ++i)
		  {
		    t = alpha1 * (sp.ubp[i] - (alpha * gp[i]));
		    dp.ubp[i] = WLZ_CLAMP(t, 0, 255);
		    bp[i] = sp.ubp[i];
		  }
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_EOO)
	  {
	    errNum = WLZ_ERR_NONE;
	  }
	}
      }
      if((errNum == WLZ_ERR_NONE) && (skip != 1))
      {
	switch(nrm)
	{
	  case 0:
	    break;
	  case 1:
	    errNum = WlzHistogramMatchObj(dObj2, hObj[1], 0, 0, -0.1, 1.1, 1);
	    break;
	  case 2:
	    errNum = WlzGreyRange(dObj2, &dMin, &dMax);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzGreySetRange(dObj2, dMin, dMax, sMin, sMax, 1);
	    }
	    break;
	  default:
	    break;

	}
      }
      if((errNum == WLZ_ERR_NONE) && dIsTiled)
      {
	errNum = WlzSetXYPlane(dObj, dObj2, p);
      }
      (void )WlzFreeObj(hObj[0]);
      hObj[0] = hObj[1];
      (void )WlzFreeObj(bObj2);
      (void )WlzFreeObj(dObj2);
      (void )WlzFreeObj(gObj2);
      (void )WlzFreeObj(sObj2);
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
    WlzFreeObj(hObj[0]);
  }
  (void )WlzFreeObj(bufObj);
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzDomainOps
* \brief	Creates a new Woolz object using the domain and (if they
* 		exist) the values with the given plane coordinate. The
* 		domain of the new object will be shared with the returned
* 		object. If the given object has values then these may also
* 		be shared. If the values are not shared (because the value
* 		table is incompatible) then a copy of the values may be made.
* \param	obj3			Given 3D domain object.
* \param	p			Plane coordinate. It is an error is
* 					the plane coordinate is outside of
* 					the given object's plane limits.
* \param	ngv			No grey values are copied to the
* 					grey table.
* \param	dstCpy			Destination pointer, set to non-zero
* 					value on return if the grey values
* 					are copied rather than shared.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzGetXYPlane(
				  WlzObject *obj3,
				  int p,
				  int ngv,
				  int *dstCpy,
				  WlzErrorNum *dstErr)
{
  int		cpy = 0;
  WlzObject	*obj2 = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(obj3 == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj3->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj3->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((p < obj3->domain.p->plane1) || (p > obj3->domain.p->lastpl))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    int	      pi;
    WlzValues val2,
              val3;

    val2.core = NULL;
    val3 = obj3->values;
    pi = p - obj3->domain.p->plane1;
    obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
		       obj3->domain.p->domains[pi], val2, NULL, NULL, &errNum);
    if(val3.core)
    {
      cpy = WlzGreyTableIsTiled(val3.core->type);
      if(cpy)
      {
	WlzGreyType	gt;
        WlzObjectType	vtt;
	WlzGreyValueWSpace *gVWSp = NULL;
	WlzGreyWSpace	gWSp;
	WlzIntervalWSpace iWSp;

	gt = WlzGreyTableTypeToGreyType(val3.core->type, NULL);
	vtt = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, gt, NULL);
        val2.v = WlzNewValueTb(obj2, vtt, val3.t->bckgrnd, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  obj2->values = WlzAssignValues(val2, NULL);
	}
	if((errNum == WLZ_ERR_NONE) && !ngv)
	{
          gVWSp = WlzGreyValueMakeWSp(obj3, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzInitGreyScan(obj2, &iWSp, &gWSp);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzIVertex3	q;

	    q.vtZ = p;
	    while((errNum == WLZ_ERR_NONE) &&
		  ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
	    {
	      WlzGreyP 	gp;

	      q.vtY = iWSp.linpos;
	      gp = gWSp.u_grintptr;
	      switch(gt)
	      {
		case WLZ_GREY_LONG:
		  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
		  {
		    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
		    *(gp.lnp)++ = gVWSp->gPtr[0].lnp[0];
		  }
		  break;
		case WLZ_GREY_INT:
		  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
		  {
		    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
		    *(gp.inp)++ = gVWSp->gPtr[0].inp[0];
		  }
		  break;
		case WLZ_GREY_SHORT:
		  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
		  {
		    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
		    *(gp.shp)++ = gVWSp->gPtr[0].shp[0];
		  }
		  break;
		case WLZ_GREY_UBYTE:
		  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
		  {
		    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
		    *(gp.ubp)++ = gVWSp->gPtr[0].ubp[0];
		  }
		  break;
		case WLZ_GREY_FLOAT:
		  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
		  {
		    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
		    *(gp.flp)++ = gVWSp->gPtr[0].flp[0];
		  }
		  break;
		case WLZ_GREY_DOUBLE:
		  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
		  {
		    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
		    *(gp.dbp)++ = gVWSp->gPtr[0].dbp[0];
		  }
		  break;
		case WLZ_GREY_RGBA:
		  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
		  {
		    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
		    *(gp.rgbp)++ = gVWSp->gPtr[0].rgbp[0];
		  }
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	    }
	    if(errNum == WLZ_ERR_EOO)
	    {
	      errNum = WLZ_ERR_NONE;
	    }
	  }
	}
      }
      else
      {
        obj2->values = WlzAssignValues(val3.vox->values[pi], NULL);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstCpy)
    {
      *dstCpy = cpy;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(obj2);
    obj2 = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj2);
}

static WlzErrorNum		WlzSetXYPlane(
				  WlzObject *obj3,
				  WlzObject *obj2,
				  int p)
{
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzGreyWSpace	gWSp;
  WlzIntervalWSpace iWSp;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  if((obj3 == NULL) || (obj2 == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((obj3->domain.core == NULL) || (obj2->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((obj3->type != WLZ_3D_DOMAINOBJ) || (obj2->type != WLZ_2D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((p < obj3->domain.p->plane1) || (p > obj3->domain.p->lastpl))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    gVWSp = WlzGreyValueMakeWSp(obj3, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(obj2, &iWSp, &gWSp );
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzIVertex3	q;

    q.vtZ = p;
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
    {
      WlzGreyP	gP;

      q.vtY = iWSp.linpos;
      gP = gWSp.u_grintptr;
      switch(gWSp.pixeltype)
      {
        case WLZ_GREY_INT:
	  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
	  {
	    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
	    gVWSp->gPtr[0].inp[0] = *(gP.inp)++;
	  }
	  break;
        case WLZ_GREY_SHORT:
	  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
	  {
	    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
	    gVWSp->gPtr[0].shp[0] = *(gP.shp)++;
	  }
	  break;
        case WLZ_GREY_UBYTE:
	  for(q.vtX = iWSp.lftpos; q.vtX <= iWSp.rgtpos; ++(q.vtX))
	  {
	    WlzGreyValueGet(gVWSp, q.vtZ, q.vtY, q.vtX);
	    gVWSp->gPtr[0].ubp[0] = *(gP.ubp)++;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  WlzGreyValueFreeWSp(gVWSp);
  return(errNum);
}

