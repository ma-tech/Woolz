#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstThreshold_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstThreshold.c
* \author       Bill Hill
* \date         May 2009
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
* \brief	Simple test for WlzThreshold, may be useful to verify
* 		memory access is OK and that there are no leaks.
* \ingroup	BinWlzTst
*/


#include <stdio.h>
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
  		usage = 0;
  FILE		*fP = NULL;
  char		*iFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*iObj = NULL;
  static char   optList[] = "h";
  const char    defFile[] = "-";

  opterr = 0;
  iFileStr = (char *)defFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if((iFileStr == NULL) || (*iFileStr == '\0'))
    {
      usage = 1;
    }
    if((usage == 0) && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
      }
      else
      {
        iFileStr = *(argv + optind);
      }
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if((iFileStr == NULL) ||
       (*iFileStr == '\0') ||
       ((fP = (strcmp(iFileStr, "-")? fopen(iFileStr, "r"): stdin)) == NULL) ||
       ((iObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read object from file (%s)\n",
                     *argv, iFileStr);
    }
    if(fP && strcmp(iFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    switch(iObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	break;
      default:
        ok = 0;
	errNum = WLZ_ERR_OBJECT_TYPE;
        (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
               "%s: Invalid object type, must be WLZ_[23]D_DOMAINOBJ (%s),\n",
		       argv[0],
		       errMsgStr);
        break;
    }
  }
  if(ok)
  {
    WlzPixelV	tV;

    tV.type = WLZ_GREY_INT;
    tV.v.inv = 0;
    do
    {
      WlzObject	*obj0;

      obj0 = WlzAssignObject(
             WlzThreshold(iObj, tV, WLZ_THRESH_HIGH, &errNum), NULL);
      if(WlzIsEmpty(obj0, NULL))
      {
	errNum = WLZ_ERR_EOO;
      }
      (void )WlzFreeObj(obj0);
      ++(tV.v.inv);
    }
    while(errNum == WLZ_ERR_NONE);
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	             "%s: WlzThreshold failed (%s),\n",
		     argv[0],
		     errMsgStr);
    }
  }
  (void )WlzFreeObj(iObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [<input object>]\n"
    "Thresholds the input object to test WlzThreshold for memory leaks etc.\n"
    "Options are:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file.\n",
    argv[0]);
  }
  return(!ok);
}
