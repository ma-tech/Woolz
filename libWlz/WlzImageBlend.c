#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzImageBlend_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzImageBlend.c
* \author       Bill Hill
* \date         December 2016
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
* \brief	Functions for blending images.
* \ingroup	WlzValueFilters
*/

#include <sys/time.h>
#include <stdio.h>
#include <Wlz.h>

#define WLZ_FAST_CODE

#ifndef restrict
#define restrict
#endif

static void			WlzImageBlendItvToBufRGB(
				  WlzUByte *buffer,
				  WlzUInt iLin,
				  WlzUInt iRgt,
				  WlzUInt iLft,
				  WlzIVertex2  bufOrg,
				  WlzUInt  bufWidth,
				  WlzUInt color);
static void			WlzImageBlendItvToBufRGBA(
				  WlzUByte *buffer,
				  WlzUInt iLin,
				  WlzUInt iRgt,
				  WlzUInt iLft,
				  WlzIVertex2  bufOrg,
				  WlzUInt  bufWidth,
				  WlzUInt color);


/*!
* \return	Woolz error code.
* \ingroup	WlzValueFilters
* \brief	Blends the give 2D object's domain into the given buffer
* 		using Porter Duff source over destination. Object values
* 		are ignored.
* \param	buffer		Given buffer with either 3 or (when
* 				alpha is used) 4 bytes per pixel.
* \param	obj		Given object to blend into the buffer.
* \param	bufOrg		Position of the buffer origin.
* \param	bufSz		Size of the buffer.
* \param	useAlpha	Use alpha channel if non-zero.
* \param	color		Colour to blend within the object's domain.
*/
WlzErrorNum			WlzImageBlendObjToBufRGBA(
				  WlzUByte *buffer,
				  WlzObject *obj,
				  WlzIVertex2 bufOrg,
				  WlzIVertex2 bufSz,
				  WlzUInt useAlpha,
				  WlzUInt color)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(buffer == NULL)
  {
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type == WLZ_2D_DOMAINOBJ)
  {
    WlzIntervalWSpace iWSp;

    if((errNum = WlzInitRasterScan(obj, &iWSp,
                                   WLZ_RASTERDIR_ILIC)) == WLZ_ERR_NONE)
    {
      while((errNum = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE)
      {
	int	bs1,
		lin,
		lft,
		rgt;

	/* Make all coordinates relative to the buffer. */
	bs1 = bufSz.vtX - 1;
	lin = iWSp.linpos - bufOrg.vtY;
	lft = iWSp.lftpos - bufOrg.vtX;
	rgt = iWSp.rgtpos - bufOrg.vtX;
#ifdef WLZ_FAST_CODE
	if(((unsigned int )lin < bufSz.vtY) && (lft < bufSz.vtX) && (rgt > 0))
#else
	if((lin > 0) && (lin < bufSz.vtY) && (lft < bufSz.vtX) && (rgt > 0))
#endif
	{
	  if(lft < 0)
	  {
	    lft = 0;
	  }
	  if(rgt > bs1)
	  {
	    rgt = bs1;
	  }
	  if(useAlpha)
	  {
	    WlzImageBlendItvToBufRGBA(buffer,
			 lin, rgt, lft, bufOrg, bufSz.vtX, color);
	  }
	  else
	  {
	    WlzImageBlendItvToBufRGB(buffer,
			 lin, rgt, lft, bufOrg, bufSz.vtX, color);
	  }
	}
      }
      if(errNum == WLZ_ERR_EOO)
      {
        errNum = WLZ_ERR_NONE;
      }
    }
  }
  return(errNum);
}

static void			WlzImageBlendItvToBufRGB(
				  WlzUByte *buffer,
				  WlzUInt iLin,
				  WlzUInt iRgt,
				  WlzUInt iLft,
				  WlzIVertex2  bufOrg,
				  WlzUInt  bufWidth,
				  WlzUInt color)
{

  int     i,
	  lnOff,
	  iWidth;
  WlzUInt r,
	  g,
	  b,
	  a,
	  a1;
  WlzUByte *restrict cLnBuf;

  a = WLZ_RGBA_ALPHA_GET(color);
  r = WLZ_RGBA_RED_GET(color) * a;
  g = WLZ_RGBA_GREEN_GET(color) * a;
  b = WLZ_RGBA_BLUE_GET(color) * a;
  a1 = 255 - a;
  iWidth = iRgt - iLft + 1;
  lnOff = (bufWidth * iLin) + iLft;
  cLnBuf = buffer + (lnOff * 3);
  for(i = 0; i < iWidth; ++i)
  {
    int         j;
    WlzUByte	*cb;

    j = i * 3;
    cb = cLnBuf + j;
    cb[0] = (r + (cb[0] * a1)) >> 8;
    cb[1] = (g + (cb[1] * a1)) >> 8;
    cb[2] = (b + (cb[2] * a1)) >> 8;
  }
}

static void			WlzImageBlendItvToBufRGBA(
				  WlzUByte *buffer,
				  WlzUInt iLin,
				  WlzUInt iRgt,
				  WlzUInt iLft,
				  WlzIVertex2  bufOrg,
				  WlzUInt  bufWidth,
				  WlzUInt color)
{
  int     i,
	  lnOff,
	  iWidth;
  WlzUInt r,
	  g,
	  b,
	  a,
	  a1,
	  a2;
  WlzUByte *restrict cLnBuf;

  a = WLZ_RGBA_ALPHA_GET(color);
  r = WLZ_RGBA_RED_GET(color) * a;
  g = WLZ_RGBA_GREEN_GET(color) * a;
  b = WLZ_RGBA_BLUE_GET(color) * a;
  a1 = 255 - a;
  a2 = a << 8;
  iWidth = iRgt - iLft + 1;
  lnOff = (bufWidth * iLin) + iLft;
  cLnBuf = buffer + (lnOff * 4);
  for(i = 0; i < iWidth; ++i)
  {
    int         j;
    WlzUByte    *cb;

    j = i * 4;
    cb = cLnBuf + j;
    cb[0] = (r  + (cb[0] * a1)) >> 8;
    cb[1] = (g  + (cb[1] * a1)) >> 8;
    cb[2] = (b  + (cb[2] * a1)) >> 8;
    cb[3] = (a2 + (cb[3] * a1)) >> 8;
  }
}

#ifdef WLZ_IMAGEBLEND_TEST

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;


static int			WlzImageBlendReadObj(
				  WlzObject **dstObj,
				  char *fileStr,
				  char *typeStr,
				  char *prog)
{
  int		ok = 1;
  FILE        	*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((fP = (strcmp(fileStr, "-"))?
	    fopen(fileStr, "r"): stdin) == NULL) ||
     ((*dstObj = WlzAssignObject(
		   WlzReadObj(fP, &errNum), NULL)) == NULL) ||
     (errNum != WLZ_ERR_NONE))
  {
    ok = 0;
    (void )fprintf(stderr,
		   "%s: failed to read %s object from file %s\n",
		   prog, typeStr, fileStr);
  }
  if(fP && strcmp(fileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  return(ok);
}

int				main(
				  int argc,
				  char *argv[])
{
  int           option,
                ok = 1,
                usage = 0,
		useTimer = 0;
  char          *inObjFileStr,
                *outObjFileStr,
		*refObjFileStr;
  int		refWidth;
  WlzUInt	color;
  WlzUInt	**ary = NULL;
  WlzIVertex2	bufOrg = {0},
  		bufSz = {0};
  const char    *errMsgStr;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzObject     *inObj = NULL,
  		*refObj = NULL,
		*outObj = NULL;
  struct timeval times[3];
  const char    *errMsg;
  static char   optList[] = "hTc:o:r:";
  const char    fileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = (char *)fileStrDef;
  outObjFileStr = (char *)fileStrDef;
  refObjFileStr = (char *)fileStrDef;
  WLZ_RGBA_RGBA_SET(color, 255, 255, 255, 255);
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'T':
        useTimer = 1;
	break;
      case 'c':
        {
	  WlzUInt c[4];

	  if((sscanf(optarg, "%d,%d,%d,%d", c + 0, c + 1, c + 2, c + 3)< 1) ||
	     ((c[0] > 255) || (c[1] > 255) || (c[2] > 255) || (c[3] > 255)))
	  {
	    usage = 1;
	  }
	  else
	  {
	    WLZ_RGBA_RGBA_SET(color, c[0], c[1], c[2], c[3]);
	  }
	}
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'r':
        refObjFileStr = optarg;
	break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  ok = (usage == 0);
  if(ok)
  {
    if(optind < argc)
    {
      if((optind + 1) != argc)
      {
	ok = 0;
        usage = 1;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  ok = (ok && WlzImageBlendReadObj(&refObj, refObjFileStr, "reference", *argv));
  ok = (ok && WlzImageBlendReadObj(&inObj, inObjFileStr, "input", *argv));
  if(ok)
  {
    WlzIBox2 iBox,
    	     rBox;

    iBox = WlzBoundingBox2I(inObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      rBox = WlzBoundingBox2I(refObj, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      iBox = WlzBoundingBoxUnion2I(rBox, iBox);
      bufSz.vtX = iBox.xMax - iBox.xMin + 1;
      bufSz.vtY = iBox.yMax - iBox.yMin + 1;
      bufOrg.vtX = iBox.xMin;
      bufOrg.vtY = iBox.yMin;
      errNum = WlzToArray2D((void ***)&ary, refObj, bufSz, bufOrg, 0,
      			    WLZ_GREY_RGBA);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	       "%s: failed to cut array from file reference object (%s)\n",
	       *argv, errMsg);
    }
  }
  if(ok)
  {
    gettimeofday(times + 0, NULL);
    errNum = WlzImageBlendObjToBufRGBA((WlzUByte *)*ary, inObj, bufOrg, bufSz,
    				       1, color);
    gettimeofday(times + 1, NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	      "%s: failed to blend input object over reference object (%s)\n",
	      *argv, errMsg);
    }
  }
  if(ok)
  {
    if(useTimer)
    {
      timersub(times + 1, times + 0, times + 2);    
      (void )printf("%s: Elapsed time = %g\n",
		    *argv, times[2].tv_sec + (0.000001 * times[2].tv_usec));
    }
    outObj = WlzFromArray2D((void **)ary, bufSz, bufOrg,
                            WLZ_GREY_RGBA, WLZ_GREY_RGBA,
                            0.0, 1.0, 0, 0, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	      "%s: failed make output object from array (%s)\n",
	      *argv, errMsg);
    }
  }
  if(ok)
  {
    FILE	*fP = NULL;

    if(((fP = (strcmp(outObjFileStr, "-")?
        fopen(outObjFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      "%s: Failed to write output object to file %s (%s).\n",
      *argv, outObjFileStr, errMsg);
    }
    if(fP != NULL)
    {
      (void )fclose(fP);
    }
  }
  WlzFreeObj(inObj);
  WlzFreeObj(refObj);
  WlzFreeObj(outObj);
  Alc2Free((void **)ary);
  if(usage)
  {
  static char   optList[] = "ahTc:o:r:";
    (void )fprintf(stderr,
    "Usage: %s%s%s%s",
    *argv,
    " [-o<output object>] [-h] [-T] [-c #,#,#,#]\n"
    "\t\t[-o <output file>] [-r <reference file>] [<input file>]\n"
    "Blends input object over reference object.\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -T  Report execution time to stderr.\n"
    "  -c  Colour value.\n"
    "  -o  Output object.\n"
    "  -r  Reference object.\n");
  }
  return(!ok);
}
#endif /* WLZ_IMAGEBLEND_TEST */
