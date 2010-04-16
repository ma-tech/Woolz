#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTiledObjFromDomain_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzTiledObjFromDomain.c
* \author       Bill Hill
* \date         April 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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
*/

/*!
\ingroup BinWlz
\defgroup wlztiledobjfromdomain WlzTiledObjFromDomain
\par Name
WlzTiledObjFromDomain  - creates an object with a tiled value table.
\par Synopsis
\verbatim
WlzTiledObjFromDomain  [-h] [-o<output file>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-b</b></td>
    <td>Background value.</td>
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
  		ok = 1,
		usage = 0;
  WlzGreyType	gType = WLZ_GREY_UBYTE;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*domObj = NULL,
  		*inObj = NULL,
  		*outObj = NULL;
  char		*inFileStr,
		*outFileStr;
  const char	*errMsg;
  const size_t	tlSz = 4096;
  static char	optList[] = "hb:g:o:",
  		inFileStrDef[] = "-",
		outFileStrDef[] = "-";

  opterr = 0;
  inFileStr = inFileStrDef;
  outFileStr = outFileStrDef;
  bgdV.v.dbv = 0.0;
  bgdV.type = WLZ_GREY_DOUBLE;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        if(sscanf(optarg, "%lg", &(bgdV.v.dbv)) != 1)
	{
	  usage = 1;
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
    WlzValues nullVal;

    nullVal.core = NULL;
    domObj = WlzMakeMain(inObj->type, inObj->domain, nullVal,
    			 NULL, NULL, &errNum);
    outObj = WlzMakeTiledValuesFromObj(inObj, tlSz, 0, gType, bgdV, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to create object with tiled values (%s).\n",
		     *argv, errMsg);
    }
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
  (void )WlzFreeObj(domObj);
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<output object>] [-h] [-o <file>] [<input object>]\n"
    "Creates an object with a tiled value table from an	object with a\n"
    "valid spatial domain.\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -b  Background value.\n"
    "  -g  Grey type specified using one of the characters:\n"
    "      l, i, s, u, f, d, r for long, intm shortm unsigned byte,\n"
    "      float, double or red-green-blue-alpha.\n"
    "  -s  Voxel size (x,y,z).\n"
    "  -o  Output tiled object.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
