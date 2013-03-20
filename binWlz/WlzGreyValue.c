#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyValue_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzGreyValue.c
* \author       Bill Hill
* \date         May 2004
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
* \brief	Gets a single grey value in an object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzgreyvalue "WlzGreyValue"
*/

/*!
\ingroup BinWlz
\defgroup wlzgreyvalue WlzGreyValue
\par Name
WlzGreyValue - gets a single grey value in an object.
\par Synopsis
\verbatim
WlzGreyValue [-b] [-v#] [-x#] [-y#] [-z#] [-h] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Indicate background or foreground with b or f following value
        written, default false.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Grey value to be set, default false.</td>
  </tr>
  <tr> 
    <td><b>-</b></td>
    <td>Column position, default 0.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Line position.</td>
  </tr>
  <tr> 
    <td><b>-z</b></td>
    <td>Plane position.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
</table>
\par Description
Either sets or gets a grey value within the given objects values.
If the -v flag is given the objects value is set and the modified
object is written to the output file, but if the flag is not set
then the grey value is written out in ascii form with a new line
character.
\par Examples
\verbatim
WlzGreyValue -x 100 -y 100 toucan.wlz
Outputs as ascii the value of the pixel at coordinate (100,100).
\endverbatim
\par File
\ref WlzGreyValue.c "WlzGreyValue.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzGreyValueGet "WlzGreyValueGet(3)"
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
  int		option,
		bgdFlag = 0,
		setVal = 0,
  		ok = 1,
		usage = 0;
  WlzDVertex3	pos;
  WlzPixelV	pix0,
  		pix1;
  WlzGreyValueWSpace *gVWSp = NULL;
  char		*outFileStr,
  		*inObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL;
  const char	*errMsg;
  static char	optList[] = "bo:v:x:y:z:h",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  pos.vtX = 0;
  pos.vtY = 0;
  pos.vtZ = 0;
  opterr = 0;
  pix0.type = WLZ_GREY_DOUBLE;
  pix0.v.dbv = 0.0;
  inObjFileStr = inObjFileStrDef;
  outFileStr = outFileStrDef;
  while(ok && (usage == 0) &&
  	((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        bgdFlag = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'v':
	setVal = 1;
        if(sscanf(optarg, "%lg", &(pix0.v.dbv)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'x':
        if(sscanf(optarg, "%lg", &(pos.vtX)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'y':
        if(sscanf(optarg, "%lg", &(pos.vtY)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'z':
        if(sscanf(optarg, "%lg", &(pos.vtZ)) != 1)
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
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (outFileStr == NULL) || (*outFileStr == '\0'))
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
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
                     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    gVWSp = WlzGreyValueMakeWSp(inObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to make workspace (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    WlzGreyValueGet(gVWSp, pos.vtZ, pos.vtY, pos.vtX);
    if(setVal)
    {
      errNum = WlzValueConvertPixel(&pix1, pix0, gVWSp->gType);
      if(errNum == WLZ_ERR_NONE)
      {
	switch(gVWSp->gType)
	{
	  case WLZ_GREY_LONG:
	    *(gVWSp->gPtr[0].lnp) = pix1.v.lnv;
	    break;
	  case WLZ_GREY_INT:
	    *(gVWSp->gPtr[0].inp) = pix1.v.inv;
	    break;
	  case WLZ_GREY_SHORT:
	    *(gVWSp->gPtr[0].shp) = pix1.v.shv;
	    break;
	  case WLZ_GREY_UBYTE:
	    *(gVWSp->gPtr[0].ubp) = pix1.v.ubv;
	    break;
	  case WLZ_GREY_FLOAT:
	    *(gVWSp->gPtr[0].flp) = pix1.v.flv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    *(gVWSp->gPtr[0].dbp) = pix1.v.dbv;
	    break;
	  case WLZ_GREY_RGBA:
	    *(gVWSp->gPtr[0].rgbp) = pix1.v.rgbv;
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
    }
    else
    {
      pix1.type = gVWSp->gType;
      pix1.v = gVWSp->gVal[0];
      errNum = WlzValueConvertPixel(&pix0, pix1, WLZ_GREY_DOUBLE);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(setVal == 0)
    {
      if((fP = (strcmp(outFileStr, "-")? fopen(outFileStr, "w"):
		                         stdout)) != NULL)
      {
	const char *borf;

	borf = (bgdFlag)? (((gVWSp->bkdFlag & 1) != 0)? " b": " f"): "";
	errNum = WLZ_ERR_NONE;
        switch(pix1.type)
	{
	  case WLZ_GREY_LONG:
	    if(fprintf(fP, "%lld%s\n", pix1.v.lnv, borf) < 2)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    break;
	  case WLZ_GREY_INT:
	    if(fprintf(fP, "%d%s\n", pix1.v.inv, borf) < 2)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    if(fprintf(fP, "%d%s\n", pix1.v.shv, borf) < 2)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    if(fprintf(fP, "%d%s\n", pix1.v.ubv, borf) < 2)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    if(fprintf(fP, "%g%s\n", pix1.v.flv, borf) < 2)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    if(fprintf(fP, "%lg%s\n", pix1.v.dbv, borf) < 2)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    if(fprintf(fP, "%d,%d,%d,%d%s\n",
	               WLZ_RGBA_RED_GET(pix1.v.rgbv),
	               WLZ_RGBA_GREEN_GET(pix1.v.rgbv),
	               WLZ_RGBA_BLUE_GET(pix1.v.rgbv),
	               WLZ_RGBA_ALPHA_GET(pix1.v.rgbv),
		       borf) < 2)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to write output file (%s).\n",
		       *argv, errMsg);
      }
    }
    else
    {
      if(((fP = (strcmp(outFileStr, "-")?
		fopen(outFileStr, "w"):
		stdout)) == NULL) ||
	 ((errNum = WlzWriteObj(fP, inObj)) != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to write output object (%s).\n",
		       *argv, errMsg);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(gVWSp)
  {
    WlzGreyValueFreeWSp(gVWSp);
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: "
    "%s [-b] [-v#] [-x#] [-y#] [-z#] [-h] [<in object>]\n"
    "Either sets or gets a grey value within th given objects values.\n"
    "If the -v flag is given the objects value is set and the modified\n"
    "object is written to the output file, but if the flag is not set\n"
    "then the grey value is written out in ascii form with a new line\n"
    "character.\n"
    "Version: %s\n"
    "Options:\n"
    "  -b    Indicate background or foreground with b or f following value\n"
    "        written (set to %s).\n"
    "  -v#   Grey value to be set (set to %g).\n"
    "  -x#   Column position (set to %g).\n"
    "  -y#   Line position (set to %g).\n"
    "  -z#   Plane position (set to %g).\n"
    "  -o#   Output file name.\n"
    "  -h    Display this usage information.\n",
    *argv,
    WlzVersion(),
    (bgdFlag)? "true": "false",
    pix0.v.dbv,
    pos.vtX, pos.vtY, pos.vtZ);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
