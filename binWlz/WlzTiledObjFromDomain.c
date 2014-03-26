#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTiledObjFromDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzTiledObjFromDomain.c
* \author       Bill Hill
* \date         April 2010
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
* \brief	Creates an object with a tiled value table from an
* 		object with a valid spatial domain.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlztiledobjfromdomain "WlzTiledObjFromDomain"
*/

/*!
\ingroup BinWlz
\defgroup wlztiledobjfromdomain WlzTiledObjFromDomain
\par Name
WlzTiledObjFromDomain  - creates an object with a tiled value table.
\par Synopsis
\verbatim
WlzTiledObjFromDomain  [-b #] [-c] [-g ] [-h] [-o<output file>]
                       [-s #,#,#] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-b</b></td>
    <td>Background value.</td>
  </tr>
  <tr>
    <td><b>-c</b></td>
    <td>Copy grey values from the given object.</td>
  </tr>
  <tr>
    <td><b>-g</b></td>
    <td>Grey type specified using one of the characters:
        l, i, s, u, f, d, r for long, intm shortm unsigned byte,
	float, double or red-green-blue-alpha.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
WlzTiledObjFromDomain creates an object with a tiled value table from
an object with a valid spatial domain.
\par Examples
\verbatim
WlzTiledObjFromDomain -o tiled.wlz in.wlz
WlzCopyToTiledObj -t tiled.wlz in.wlz
\endverbatim
Creates a new object (tiled.wlz) with the same domain (and voxel dimensions)
as the input object (in.wlz) but a with tiled value table . The values in
the tiles are set in this example using WlzCopyToTiledObj.
\par File
\ref WlzTiledObjFromDomain.c "WlzTiledObjFromDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCopyToTiledObj "WlzCopyToTiledObj(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
		copy = 0,
  		ok = 1,
		usage = 0,
		voxSzSet = 0;
  WlzGreyType	gType = WLZ_GREY_UBYTE;
  WlzFVertex3	voxSz;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  char		*inFileStr,
		*outFileStr;
  int		iBuf[4];
  const char	*errMsg;
  const size_t	tlSz = 4096;
  static char	optList[] = "chb:g:o:s:",
  		inFileStrDef[] = "-",
		outFileStrDef[] = "-";

  opterr = 0;
  inFileStr = inFileStrDef;
  outFileStr = outFileStrDef;
  bgdV.v.dbv = 0.0;
  bgdV.type = WLZ_GREY_DOUBLE;
  voxSz.vtX = voxSz.vtY = voxSz.vtZ = 1.0f;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'c':
        copy = 1;
	break;
      case 'b':
	bgdV.type = gType;
	switch(gType)
	{

	  case WLZ_GREY_LONG:
	    if(sscanf(optarg, "%lld", &(bgdV.v.lnv)) != 1)
	    {
	      usage = 1;
	    }
	    break;
	  case WLZ_GREY_INT:
	    if(sscanf(optarg, "%d", &(bgdV.v.inv)) != 1)
	    {
	      usage = 1;
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    if((sscanf(optarg, "%d", iBuf) != 1) ||
	       (iBuf[0] < SHRT_MIN) || (iBuf[0] > SHRT_MAX))
	    {
	      usage = 1;
	    }
	    else
	    {
	      bgdV.v.shv = iBuf[0];
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    if((sscanf(optarg, "%d", iBuf) != 1) ||
	       (iBuf[0] < 0) || (iBuf[0] > 255))
	    {
	      usage = 1;
	    }
	    else
	    {
	      bgdV.v.ubv = iBuf[0];
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    if(sscanf(optarg, "%g", &(bgdV.v.flv)) != 1)
	    {
	      usage = 1;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    if(sscanf(optarg, "%lg", &(bgdV.v.dbv)) != 1)
	    {
	      usage = 1;
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    iBuf[0] = iBuf[1] = iBuf[2] = iBuf[3] = 0;
	    if((sscanf(optarg, "%d,%d,%d,%d",
	               iBuf + 0, iBuf + 1, iBuf + 2, iBuf + 3) < 3) ||
	       (iBuf[0] < 0) || (iBuf[0] > 255) ||
	       (iBuf[1] < 0) || (iBuf[1] > 255) ||
	       (iBuf[2] < 0) || (iBuf[2] > 255) ||
	       (iBuf[3] < 0) || (iBuf[3] > 255))
	    {
	      usage = 1;
	    }
	    else
	    {
	      WLZ_RGBA_RGBA_SET(bgdV.v.rgbv,
	                        iBuf[0], iBuf[1], iBuf[2], iBuf[3]);
	    }
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'g':
        switch(*optarg)
	{
	  case 'l':
	    gType = WLZ_GREY_LONG;
	    break;
	  case 'i':
	    gType = WLZ_GREY_INT;
	    break;
	  case 's':
	    gType = WLZ_GREY_SHORT;
	    break;
	  case 'u':
	    gType = WLZ_GREY_UBYTE;
	    break;
	  case 'f':
	    gType = WLZ_GREY_FLOAT;
	    break;
	  case 'd':
	    gType = WLZ_GREY_DOUBLE;
	    break;
	  case 'r':
	    gType = WLZ_GREY_RGBA;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 's':
	voxSzSet = 1;
	if(sscanf(optarg, "%g,%g,%g",
	          &(voxSz.vtX), &(voxSz.vtY), &(voxSz.vtZ)) != 3)
	{
	  usage = 1;
	}
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
        case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
        case WLZ_3D_DOMAINOBJ:
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
    WlzObject *domObj = NULL;
    WlzValues nullVal;

    nullVal.core = NULL;
    domObj = WlzMakeMain(inObj->type, inObj->domain,
                         (copy == 0)? nullVal: inObj->values,
    			 NULL, NULL, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if((voxSzSet != 0) && (domObj->type == WLZ_3D_DOMAINOBJ))
      {
	domObj->domain.p->voxel_size[0] = voxSz.vtX;
	domObj->domain.p->voxel_size[1] = voxSz.vtY;
	domObj->domain.p->voxel_size[2] = voxSz.vtZ;
      }
      outObj = WlzMakeTiledValuesFromObj(domObj, tlSz, copy, gType, bgdV,
      					 &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to create object with tiled values (%s).\n",
		     *argv, errMsg);
    }
    (void )WlzFreeObj(domObj);
  }
  if(ok && (outFileStr != NULL))
  {
    if(((fP = (strcmp(outFileStr, "-")? fopen(outFileStr, "w"):
                                        stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write tiled object to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    if(fP!= NULL)
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(outObj);
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%s",
    *argv,
    " [-o<output object>] [-h] [-b #] [-g #] [-s #,#,#]\n"
    "                             [<input object>]\n"
    "Creates an object with a tiled value table from an object with a\n"
    "valid spatial domain.\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -b  Background value. If the grey type is RGBA then the channel\n"
    "      background values should be comma seperated and in the order\n"
    "      RGBA.\n"
    "  -c  Copy the values from the given object to the tiled object.\n"
    "  -g  Grey type specified using one of the characters:\n"
    "      l, i, s, u, f, d, r for long, intm shortm unsigned byte,\n"
    "      float, double or RGBA.\n"
    "  -h  Prints this usage information.\n"
    "  -s  Voxel size (x,y,z).\n"
    "  -o  Output tiled object.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
