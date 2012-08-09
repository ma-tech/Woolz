#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCutObjToBox_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCutObjToBox.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Cuts  out  a  filled  box from the given object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcutobjtobox "WlzCutObjToBox"
*/

/*!
\ingroup BinWlz
\defgroup wlzcutobjtobox WlzCutObjToBox
\par Name
WlzCutObjToBox - Cuts  out  a  filled  box (rectangle or
cuboid) from the given object. Background values filled
with single a value or random distribution.
\par Synopsis
\verbatim
WlzCutObjToBox [-dfhisuN] [-o<output object file>]
           [-M<noise mean>] [-S<noise std dev>]
           [-x<x min>,<x max>] [-y<y min>,<y max>]
           [-z<z min>,<z max>] [<input object file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file</td>
  </tr>
  <tr> 
    <td><b>-N</b></td>
    <td>Fill Background with Gaussian noise.</td>
  </tr>
  <tr> 
    <td><b>-M</b></td>
    <td>Mean of Gaussian noise.</td>
  </tr>
  <tr> 
    <td><b>-S</b></td>
    <td>Standard deviation of Gaussian noise.</td>
  </tr>
  <tr> 
    <td><b>-x</b></td>
    <td>Column box limits.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Line box limits.</td>
  </tr>
  <tr> 
    <td><b>-z</b></td>
    <td>Plane box limits</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Output object has <b>double</b> value table
	    elements.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Output object has <b>float</b> value table
	    elements.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Output object has <b>integer</b> value table
	    elements.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Output object has <b>short</b> value table
	    elements.</td>
  </tr>
  <tr> 
    <td><b>-u</b></td>
    <td>Output object has <b>unsigned byte</b> value table
	    elements.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation (not enabled).</td>
  </tr>
</table>
By  default  the  input  object  is read from the standard
       input and the output object is  written  to  the  standard
       output.   Box limit values are comma separated and default
       to the bounding box of  the  input  objects  domain.   The
       default  grey  value  type is the same as the input object
       grey value type.
\par Description
Cuts a WLZ_2D_DOMAINOBJ or WLZ_3D_DOMAINOBJ rectangular or
cuboid object so that it fills the given box.
\par Examples
A simple example of using WlzCutObjToBox to cut a rectangular or
cuboid object from the object that is read from the standard input.
The output object, which is written to the standard output,
fills the bounding box of the input object.
\verbatim
WlzCutObjToBox <infile.wlz >outfile.wlz
\endverbatim
An example in which the input object is used to create a
0-256 rectangular or cuboid object.
\verbatim
WlzCutObjToBox -x0,256, -y0,256 -z0,256 -o outfile.wlz infile.wlz
\endverbatim  

\par File
\ref WlzCutObjToBox.c "WlzCutObjToBox.c"
\par See Also
\ref WlzCutObjToBox2D "WlzCutObjToBox2D(3)"
\ref WlzCutObjToBox3D "WlzCutObjToBox3D(3)"
\ref WlzClipObjToBox2D "WlzClipObjToBox2D(3)"
\ref WlzClipObjToBox3D "WlzClipObjToBox3D(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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

int             main(int argc, char **argv)
{
  int		idx,
  		option,
		noiseFlg = 0,
		ok = 1,
		usage = 0;
  double	noiseMu = 0.0,
  		noiseSigma = 0.0;
  WlzIBox3	cutBox,
		cutSet;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzGreyType	dstGreyType = WLZ_GREY_ERROR;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  int		cutVal[2];
  char		*cutStr[2];
  const char	*errMsg;
  static char	optList[] = "iNsufdM:S:o:x:y:z:h",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  cutSet.xMin = 0;
  cutSet.xMax = 0;
  cutSet.yMin = 0;
  cutSet.yMax = 0;
  cutSet.zMin = 0;
  cutSet.zMax = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'i':
        dstGreyType = WLZ_GREY_INT;
	break;
      case 's':
        dstGreyType = WLZ_GREY_SHORT;
	break;
      case 'u':
        dstGreyType = WLZ_GREY_UBYTE;
	break;
      case 'f':
        dstGreyType = WLZ_GREY_FLOAT;
	break;
      case 'd':
        dstGreyType = WLZ_GREY_DOUBLE;
	break;
      case 'N':
        noiseFlg = 1;
        break;
      case 'M':
        if(sscanf(optarg, "%lg", &noiseMu) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'S':
        if(sscanf(optarg, "%lg", &noiseSigma) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'x':
      case 'y':
      case 'z':
	if(optarg)
	{
	  while(*optarg && isspace(*optarg))
	  {
	    ++optarg;
	  }
	  if(*optarg == ',')
	  {
	    cutStr[0] = NULL;
	    cutStr[1] = strtok(optarg, ",");
	  }
	  else
	  {
	    cutStr[0] = strtok(optarg, ",");
	    cutStr[1] = strtok(NULL, ",");
	  }
	  if((cutStr[0] == NULL) && (cutStr[1] == NULL))
	  {
	    usage = 1;
	    ok = 0;
	  }
	  else
	  {
	    idx = 0;
	    while(ok && (idx < 2))
	    {
	      if(cutStr[idx] && (sscanf(cutStr[idx], "%d",
					 cutVal + idx) != 1))
	      {
		usage = 1;
		ok = 0;
	      }
	      ++idx;
	    }
	  }
	}
	if(ok)
	{
	  switch(option)
	  {
	    case 'x':
	      if(cutStr[0])
	      {
		cutSet.xMin = 1;
		cutBox.xMin = cutVal[0];
	      }
	      if(cutStr[1])
	      {
		cutSet.xMax = 1;
		cutBox.xMax = cutVal[1];
	      }
	      break;
	    case 'y':
	      if(cutStr[0])
	      {
		cutSet.yMin = 1;
		cutBox.yMin = cutVal[0];
	      }
	      if(cutStr[1])
	      {
		cutSet.yMax = 1;
		cutBox.yMax = cutVal[1];
	      }
	      break;
	    case 'z':
	      if(cutStr[0])
	      {
		cutSet.zMin = 1;
		cutBox.zMin = cutVal[0];
	      }
	      if(cutStr[1])
	      {
		cutSet.zMax = 1;
		cutBox.zMax = cutVal[1];
	      }
	      break;
	  }
	}
	break;
      case 'h':
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
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzReadObj(fP, &errNum)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    switch(inObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	if(cutSet.xMin == 0)
	{
          cutBox.xMin = inObj->domain.i->kol1;
	}
	if(cutSet.xMax == 0)
	{
          cutBox.xMax = inObj->domain.i->lastkl;
	}
	if(cutSet.yMin == 0)
	{
          cutBox.yMin = inObj->domain.i->line1;
	}
	if(cutSet.yMax == 0)
	{
          cutBox.yMax = inObj->domain.i->lastln;
	}
	cutBox.zMin = 0;
	cutBox.zMax = 0;
        break;
      case WLZ_3D_DOMAINOBJ:
	if(cutSet.xMin == 0)
	{
          cutBox.xMin = inObj->domain.p->kol1;
	}
	if(cutSet.xMax == 0)
	{
          cutBox.xMax = inObj->domain.p->lastkl;
	}
	if(cutSet.yMin == 0)
	{
          cutBox.yMin = inObj->domain.p->line1;
	}
	if(cutSet.yMax == 0)
	{
          cutBox.yMax = inObj->domain.p->lastln;
	}
	if(cutSet.zMin == 0)
	{
          cutBox.zMin = inObj->domain.p->plane1;
	}
	if(cutSet.zMax == 0)
	{
          cutBox.zMax = inObj->domain.p->lastpl;
	}
        break;
      default:
        ok = 0;
	errMsg = WlzStringFromObjType(inObj, NULL);
	(void )fprintf(stderr,
		       "%s: invalid object type %s.\n", 
		       *argv,
		       errMsg? errMsg: "unknown");
    }
  }
  if(ok)
  {
    if(inObj->domain.core == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: invalid object, null domain\n",
		     *argv);
    }
  }
  if(ok)
  {
    if(dstGreyType == WLZ_GREY_ERROR)
    {
      if(inObj->values.core == NULL)
      {
        dstGreyType = WLZ_GREY_UBYTE;
      }
      else
      {
	dstGreyType = WlzGreyTypeFromObj(inObj, NULL);
      }
    }
    if(((outObj = WlzCutObjToBox3D(inObj, cutBox, dstGreyType,
    			           noiseFlg, noiseMu, noiseSigma,
				   &errNum)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to cut object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write cut object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(outObj)
  {
    WlzFreeObj(outObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-dfhiNsu] [-o<out object>]\n"
    "       [-M <mean>] [-S <std dev>]\n"
    "       [-x<x min>,<x max>] [-y<y min>,<y max>] [-z<z min>,<z max>]\n"
    "       [<in object>]\n"
    "Cuts a  Woolz object using specified cut box.\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -d  Double output values.\n"
    "  -f  Float output values.\n"
    "  -h  Help, prints this usage message.\n"
    "  -i  Int output values.\n"
    "  -N  Use gausian noise for background.\n"
    "  -s  Short output values.\n"
    "  -u  Unsigned byte output values.\n"
    "  -M  Mean of gaussian noise.\n"
    "  -o  Output object file name.\n"
    "  -S  Standard deviation of gaussian noise.\n"
    "  -x  Column cut box limits.\n"
    "  -y  Line cut box limits.\n"
    "  -z  Plane cut box limits.\n"
    "All cut limits default to the input objects domain limits.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -i -o cut.wlz -x 100,800 -y ,800 myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, cut and written\n"
    "to cut.wlz, with int values.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
