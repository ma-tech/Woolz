#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCopyToTiledObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCopyToTiledObj.c
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
* \brief	Sets values in a tiled object using the values from
* 		another object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcopytotiledobj "WlzCopyToTiledObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzcopytotiledobj WlzCopyToTiledObj
\par Name
WlzCopyToTiledObj  - sets values in a tiled object from another object.
\par Synopsis
\verbatim
WlzCopyToTiledObj  [-h] -t <tiled fobject> [-x #] [-y #] [-z #] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>tiled object.</td>
  </tr>
  <tr>
    <td><b>-x</b></td>
    <td>Column offset.</td>
  </tr>
  <tr>
    <td><b>-y</b></td>
    <td>Line offset.</td>
  </tr>
  <tr>
    <td><b>-z</b></td>
    <td>Plane offset.</td>
  </tr>
</table>
\par Description
WlzCopyToTiledObj sets values in a tiled object from another object.
\par Examples
\verbatim
WlzTiledObjFromDomain -o tiled.wlz in.wlz
WlzCopyToTiledObj -t tiled.wlz in.wlz
\endverbatim
Creates a new object (tiled.wlz) with the same domain as the input object
(in.wlz) but a tiled value table. The values in the tiles are set using
WlzCopyToTiledObj.
\par File
\ref WlzCopyToTiledObj.c "WlzCopyToTiledObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzTiledObjFromDomain "WlzTiledObjFromDomain(1)"
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
  		nFiles = 0,
		usage = 0;
  WlzIVertex3	offset;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*tlObj = NULL;
  char		*tlFileStr;
  const char	*errMsg;
  static char	optList[] = "ht:x:y:z:",
		tlFileStrDef[] = "-";

  opterr = 0;
  tlFileStr = tlFileStrDef;
  WLZ_VTX_3_SET(offset, 0, 0, 0);
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 't':
        tlFileStr = optarg;
	break;
      case 'x':
        if(sscanf(optarg, "%d", &(offset.vtX)) != 1)
	{
	  usage = 1;
	}
	break;
      case 'y':
        if(sscanf(optarg, "%d", &(offset.vtY)) != 1)
	{
	  usage = 1;
	}
	break;
      case 'z':
        if(sscanf(optarg, "%d", &(offset.vtZ)) != 1)
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
  if(usage == 0)
  {
    if((nFiles = argc - optind) <= 0)
    {
      usage = 1;
    }
  }
  ok = !usage;
  if(ok)
  {
    FILE	*fP = NULL;

    errNum = WLZ_ERR_READ_EOF;
    if(((fP = fopen(tlFileStr, "r+")) == NULL) ||
       ((tlObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to open tiled object in file %s (%s).\n",
		     *argv, tlFileStr, errMsg);
    }
    if(fP)
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    if(tlObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if(tlObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else if(tlObj->values.core == NULL)
    {
      errNum = WLZ_ERR_VALUES_NULL;
    }
    else if(WlzGreyTableIsTiled(tlObj->values.core->type) == 0)
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: invalid tiled object in file %s (%s).\n",
		     *argv, tlFileStr, errMsg);
    }
  }
  if(ok)
  {
    int	idx;

    for(idx = 0; idx < nFiles; ++idx)
    {
      char	*inFileStr;

      if((inFileStr = *(argv + optind + idx)) != NULL)
      {
        FILE	*fP = NULL;

        if((fP = (strcmp(inFileStr, "-")?
	         fopen(inFileStr, "r"): stdin)) == NULL)
        {
	  ok = 0;
	  (void )fprintf(stderr,
	                 "%s: Failed to open input file %s.\n",
			 argv[0], inFileStr);
	}
	if(ok)
	{
	  WlzObject	*obj1 = NULL,
	  		*obj2 = NULL,
			*inObj = NULL;

	  inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(inObj == NULL)
	    {
	      errNum = WLZ_ERR_OBJECT_NULL;
	    }
	    else
	    {
	      switch(inObj->type)
	      {
	        case WLZ_2D_DOMAINOBJ:
		  obj1 = WlzConstruct3DObjFromObj(1, &inObj, 0,
		  				  1.0f, 1.0f, 1.0f,
						  &errNum);
		  break;
	        case WLZ_3D_DOMAINOBJ:
		  obj1 = WlzAssignObject(inObj, NULL);
		  break;
	        default:
		  errNum = WLZ_ERR_OBJECT_TYPE;
		  break;
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    obj2 = WlzShiftObject(obj1, offset.vtX, offset.vtY, offset.vtZ,
	                          &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzCopyObjectGreyValues(tlObj, obj2);
	  }
          (void )WlzFreeObj(obj1);
          (void )WlzFreeObj(obj2);
          (void )WlzFreeObj(inObj);
	}
      }
    }
  }
  (void )WlzFreeObj(tlObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-h] [-t <tiled object>] [-x #] [-y #] [-z #]\n"
    "                         [<input objects>]\n"
    "Sets the values in the tiled object using the values in the given\n"
    "input object(s).\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -t  The tiled object.\n"
    "  -x  Column offset.\n"
    "  -y  Line offset.\n"
    "  -z  Plane offset.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
