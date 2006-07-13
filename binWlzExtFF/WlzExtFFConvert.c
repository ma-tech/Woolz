#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFConvert_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzExtFF/WlzExtFFConvert.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Woolz binary which uses the external file formats
* 		extension to convert object data file formats.
* \ingroup	BinWlzExtFF
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzextffconvert "WlzExtFFConvert"
*/

/*!
\ingroup BinWlzExtFF
\defgroup wlzextffconvert WlzExtFFConvert
\par Name
WlzExtFFConvert - converts image objects between file formats.
\par Synopsis
\verbatim
WlzExtFFConvert [-h] [-s] [-b<background>]
                [-f<input format>] [-F<output format>]
		[-x<x size>] [-y<y size>] [-z<z size>]
		[-o<output file>] [<input file>)]

\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Split labeled volumes into domains.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Set background to value.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Set size of minimum dimension, i.e. min of width or height.</td>
  </tr>
  <tr> 
    <td><b>-D</b></td>
    <td>Set size of maximum dimension, i.e. max of width or height.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Input file format.</td>
  </tr>
  <tr> 
    <td><b>-F</b></td>
    <td>Ouput file format.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
  <tr> 
    <td><b>-x</b></td>
    <td>Column (x) voxel/pixel size.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Line (y) voxel/pixel size.</td>
  </tr>
  <tr> 
    <td><b>-z</b></td>
    <td>Plane (z) voxel/pixel size.</td>
  </tr>
</table>
\par Description
Converts objects between one file format and another,
neither of which need be the Woolz data file format.
By default objects are read from the standard input and written to
the standard output.

The valid file formats are:
<table width="500" border="0">
  <tr>
  <td><b>Description</b></td> <td><b>Format</b></td> <td><b>Extension</b></td>
  </tr>
  <tr>
  <td>Amira Lattice</td> <td>am</td> <td>.am</td>
  </tr>
  <tr>
  <td>Analyze</td> <td>anl</td> <td>.hdr/.img</td>
  </tr>
  <tr>
  <td>Microsoft Bitmap</td> <td>bmp</td> <td>.bmp</td>
  </tr>
  <tr>
  <td>Stanford Density</td> <td>den</td> <td>.den</td>
  </tr>
  <tr>
  <td>ICS</td> <td>ics</td> <td>.ics/.ids</td>
  </tr>
  <tr>
  <td>PNM</td> <td>pnm</td> <td>.pgm</td>
  </tr>
  <tr>
  <td>BioRad Confocal</td> <td>pic</td> <td>.pic</td>
  </tr>
  <tr>
  <td>SLC</td> <td>slc</td> <td>.slc</td>
  </tr>
  <tr>
  <td>Sunvision VFF</td> <td>vff</td> <td>.vff</td>
  </tr>
  <tr>
  <td>Visualization Toolkit</td> <td>vtk</td> <td>.vtk</td>
  </tr>
  <tr>
  <td>IPLab</td> <td>ipl</td> <td>.ipl</td>
  </tr>
  <tr>
  <td>TIFF</td> <td>tif</td> <td>.tif</td>
  </tr>
  <tr>
  <td>JPEG</td> <td>jpg</td> <td>.jpg</td>
  </tr>
  <tr>
  <td>GIF</td> <td>gif</td> <td>.gif</td>
  </tr>
  <tr>
  <td>MRC HGU Woolz</td> <td>wlz</td> <td>.wlz</td>
  </tr>
</table>

The TIFF file format must be read/written from/to a file
and not from/to the standard input/output.

Not all formats can retain position information i.e. they can
only keep the size of the bounding box. In these formats the
size of the bounding box is maintained but the position is set to
(0,0). This implies that conversion back to woolz will, in general,
result in a shifted image, i.e. registration is lost. Most 3D
formats encode this data, of the 2D formats only woolz can retain
all offsets, TIFF can only encode positive offsets.

\par Examples

\verbatim
WlzExtFFConvert -f wlz -F slc <in.wlz >out.slc
\endverbatim
Converts the Woolz object in.wlz to an SLC data file out.slc.

\verbatim
WlzExtFFConvert -f den -F pnm -o out.pgm in.den
\endverbatim
Converts the Stanford density file in.den to a series of PGM files
each with a name of the form out000001.pgm where the number
encodes the image plane and a control file which specifies the
volume origin, size and voxel dimensions.
\par File
\ref WlzFacts.c "WlzFacts.c"
\par See Also
\ref BinWlzExtFF "WlzIntro(1)"
\ref WlzEffReadObj "WlzEffReadObj(3)"
\ref WlzEffWriteObj "WlzEffWriteObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Wlz.h>
#include <WlzExtFF.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

static WlzEffFormat WlzEffStringExtFromFileName(const char *fNameStr)
{
  int		len;
  char		*dot,
  		*ext,
		*sep;
  WlzEffFormat	fmt = WLZEFF_FORMAT_NONE;

  if(((len = strlen(fNameStr)) >= 3) &&  /* Minimum of 3 chars. */
     (*(fNameStr + len - 1) != '/'))     /* Directory not plain file. */
  {
    sep = strrchr(fNameStr, '/');
    if(((dot = strrchr(fNameStr, '.')) != NULL) &&
       ((sep == NULL) || ((dot - sep) > 1)))
    {
      ext = ++dot;
      if(ext != NULL)
      {
        fmt = WlzEffStringExtToFormat(ext);
      }
    }
  }
  return(fmt);
}

int             main(int argc, char **argv)
{
  int		bgdFlag = 0,
		split = 0,
  		option,
		ok = 1,
		usage = 0;
  int		minDimension, maxDimension;
  int		linearFlg;
  WlzEffFormat  inFmt = WLZEFF_FORMAT_NONE,
  		outFmt = WLZEFF_FORMAT_NONE;
  WlzFVertex3	voxelSz;
  WlzIVertex3	voxelSzFlags;
  FILE		*fP;
  WlzObject     *inObj = NULL;
  char		*fStr,
  		*outObjFileStr,
  		*inObjFileStr;
  WlzPixelV	tV0,
  		bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char    *errMsg;
  static char	optList[] = "b:d:D:f:F:Lo:x:y:z:hs",
		outObjFileStrDef[] = "-",
		inObjFileStrDef[] = "-";
 
  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  voxelSzFlags.vtX = voxelSzFlags.vtY = voxelSzFlags.vtZ = 0;
  minDimension = 0;
  maxDimension = 0;
  linearFlg = 0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        if(sscanf(optarg, "%lg", &(bgdV.v.dbv)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	else
	{
	  bgdV.type = WLZ_GREY_DOUBLE;
	  bgdFlag = 1;
	}
        break;
      case 'd':
	if(sscanf(optarg, "%d", &minDimension) != 1){
	  usage = 1;
	  ok = 0;
	}
	if( minDimension <= 0 ){
	  fprintf(stderr, "%s: minDimension must be > zero\n");
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'D':
	if(sscanf(optarg, "%d", &maxDimension) != 1){
	  usage = 1;
	  ok = 0;
	}
	if( maxDimension <= 0 ){
	  fprintf(stderr, "%s: maxDimension must be > zero\n");
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'f':
        if((inFmt = WlzEffStringExtToFormat(optarg)) == 0)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'F':
        if((outFmt = WlzEffStringExtToFormat(optarg)) == 0)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'L':
	linearFlg = 1;
	break;
      case 'x':
	voxelSzFlags.vtX = 1;
	if(sscanf(optarg, "%g", &(voxelSz.vtX)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'y':
	voxelSzFlags.vtY = 1;
	if(sscanf(optarg, "%g", &(voxelSz.vtY)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'z':
	voxelSzFlags.vtZ = 1;
	if(sscanf(optarg, "%g", &(voxelSz.vtZ)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 's':
        split = 1;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
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
  if(ok)
  {
    /* Try and determine file formats from extensions if not known already. */
    if(inFmt == WLZEFF_FORMAT_NONE)
    {
      if((inFmt = WlzEffStringExtFromFileName(inObjFileStr)) == 0)
      {
        usage = 1;
	ok = 0;
      }
    }
    if(outFmt == WLZEFF_FORMAT_NONE)
    {
      if((outFmt = WlzEffStringExtFromFileName(outObjFileStr)) == 0)
      {
        usage = 1;
	ok = 0;
      }
    }
  }
  if(ok)
  {
    /* if the output format is TIFF then a file name must be
       given */
    if((outFmt == WLZEFF_FORMAT_TIFF) && !strcmp(outObjFileStr, "-"))
    {
      usage = 1;
      ok = 0;
    }
  }
  if(ok)
  {
    if(strcmp(inObjFileStr, "-") == 0)
    {
      fP = stdin;
      fStr = NULL;
    }
    else
    {
      fP = NULL;
      fStr = inObjFileStr;
    }
    errNum = WLZ_ERR_READ_EOF;
    if((inObj= WlzAssignObject(WlzEffReadObj(fP, fStr, inFmt, split,
    					     &errNum), NULL)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: failed to read object from file(s) %s (%s)\n",
                     *argv, inObjFileStr, errMsg);
    }
  }
  if(ok && bgdFlag)
  {
    /* Set background if requested on command line */
    tV0 = WlzGetBackground(inObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzValueConvertPixel(&bgdV, bgdV, tV0.type);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum =  WlzSetBackground(inObj, bgdV);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
               "%s: failed to set background of object from file(s) %s (%s)\n",
		     *argv, inObjFileStr, errMsg);
    }
  }
  if(ok && (inObj->type == WLZ_3D_DOMAINOBJ) && inObj->domain.core &&
     (inObj->domain.core->type == WLZ_PLANEDOMAIN_DOMAIN))
  {
    if(voxelSzFlags.vtX)
    {
      inObj->domain.p->voxel_size[0] = voxelSz.vtX;
    }
    if(voxelSzFlags.vtY)
    {
      inObj->domain.p->voxel_size[1] = voxelSz.vtY;
    }
    if(voxelSzFlags.vtZ)
    {
      inObj->domain.p->voxel_size[2] = voxelSz.vtZ;
    }
  }
  /* option to change output dimension, only 2D for now */
  if(ok && (maxDimension || minDimension) && (inObj->type == WLZ_2D_DOMAINOBJ)){
    double 	scale;
    int		xSize, ySize;
    WlzAffineTransform	*trans;
    WlzObject	*tmpObj;

    xSize = inObj->domain.i->lastkl - inObj->domain.i->kol1 + 1;
    ySize = inObj->domain.i->lastln - inObj->domain.i->line1 + 1;
    if( maxDimension ){
      scale = ((double) maxDimension) / WLZ_MAX(xSize, ySize);
    }
    else {
      scale = ((double) minDimension) / WLZ_MIN(xSize, ySize);
    }
    trans = WlzAffineTransformFromScale(WLZ_TRANSFORM_2D_AFFINE,
					scale, scale, scale, &errNum);
    if( tmpObj = WlzAffineTransformObj(inObj, trans,
				       linearFlg?WLZ_INTERPOLATION_LINEAR:
				       WLZ_INTERPOLATION_NEAREST,
				       &errNum) ){
      tmpObj = WlzAssignObject(tmpObj, &errNum);
      WlzFreeObj(inObj);
      inObj = tmpObj;
    }
    else {
      ok = 0;
    }
    WlzFreeAffineTransform( trans );
  }
      
  if(ok)
  {
    if(strcmp(outObjFileStr, "-") == 0)
    {
      fP = stdout;
      fStr = NULL;
    }
    else
    {
      fP = NULL;
      fStr = outObjFileStr;
    }
    if((errNum = WlzEffWriteObj(fP, fStr,
    				inObj, outFmt)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: failed to write object to file(s) %s (%s).\n",
                     *argv, outObjFileStr, errMsg);
    }
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(!ok)
  {
    if(usage)
    {
      (void )fprintf(
	stderr,
	"Usage: %s%s%s%s%s%s%s%s\n",
	*argv,
	" [-h] [-s] [-b<background>]\n"
	"        [-d<min-dimension>] [-D<max-dimension>]\n"
	"        [-f<input format>] [-F<output format>]\n"
	"        [-x<x size>] [-y<y size>] [-z<z size>]\n"
	"        [-o<output file>] [<input file>)]\n"
	"Converts objects between one file format and another,\n"
	"neither of which need be the Woolz data file format.\n"
	"Options are:\n"
	"  -h    Help, prints this usage information.\n"
	"  -s    Split labeled volumes into domains.\n"
	"  -b#   Set background to value,\n"
	"  -d#   Set size of minimum dimension, i.e. min of width or height\n"
	"  -D#   Set size of maximum dimension, i.e. max of width or height\n"
	"  -f#   Input file format.\n"
	"  -F#   Ouput file format.\n"
	"  -o#   Output file name.\n"
	"  -x#   X voxel/pixel size.\n"
	"  -y#   Y voxel/pixel size.\n"
	"  -z#   Z voxel/pixel size.\n"
	"Valid formats are:\n"
	"  Description                 Fmt     Extension\n"
	"  ***********                 ***     *********\n"
	"  Amira Lattice               am      .am\n"
	"  Analyze                     anl     .hdr/.img\n"
	"  Microsoft Bitmap            bmp     .bmp\n"
	"  Stanford Density            den     .den\n"
	"  ICS                         ics     .ics/.ids\n"
	"  PNM                         pnm     .pgm\n"
	"  BioRad Confocal             pic     .pic\n"
	"  SLC                         slc     .slc\n"
	"  Sunvision VFF               vff     .vff\n"
	"  Visualization Toolkit VTK   vtk     .vtk\n"
	"  IPLab                       ipl     .ipl\n"
	"  TIFF                        tif     .tif\n"
	"  JPEG                        jpg     .jpg\n"
	"  GIF                         gif     .gif\n"
	"  MRC HGU Woolz               wlz     .wlz\n",
	"Simple example:\n  ",
	*argv,
	" -f wlz -F slc <in.wlz >out.slc\n"
	"  Converts the Woolz object in.wlz to an SLC data file out.slc\n",
	"More complex example:\n  ",
	*argv,
	" -f den -F pnm -o out.pgm in.den\n"
	"  Converts the Stanford density file in.den to a series of PGM files\n"
	"  each with a name of the form out000001.pgm where the number\n"
	"  encodes the image plane and a control file which specifies the\n"
	"  volume origin, size and voxel dimensions.\n"
	"By default objects are read from the standard input and written to\n"
	"the standard output.\n"
	"File formats which use more than one file can not be read or written\n"
	"using the standard input or standard output.\n"
	"The TIFF file format must be read/written from/to a file i.e. not\n"
	"from/to stdin or stdout\n"
	"Not all formats can retain position information i.e. they can\n"
	"only keep the size of the bounding box. In these formats the\n"
	"size of the bounding box is maintained but the position is set to\n"
	"(0,0). This implies that conversion back to woolz will, in general,\n"
	"result in a shifted image, i.e. registration is lost. Most 3D\n"
	"formats encode this data, of the 2D formats only woolz can retain\n"
	"all offsets, TIFF can only encode positive offsets.\n"
	);
    }
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
