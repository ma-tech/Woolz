#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzContourCut_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzContourCut.c
* \author       Bill Hill
* \date         February 2010
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
* \brief	Creates a new contour from the two given contours in which
* 		no simplices of the given model contour which intersect
* 		the given knife contour are included.
* \ingroup	BinWlz	
* \ref wlzcontourcut "WlzContourCut"
*/

/*!
\ingroup BinWlz
\defgroup wlzcontourcut WlzContourCut
\par Name
WlzContourCut - cuts one contour uning a second.
\par Synopsis
\verbatim
WlzContourCut [-k <knife>] [-o <out>] [-h] [<model>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr> 
    <td><b>-k</b></td>
    <td>The knife contour object.</td>
  </tr>
</table>
\par Description
Given model and knife contours; creates a new contour which consists
of all simplices of the model that do not intersect simplices of the
knife. This may be used to cut into a contour removing sections from
it.
\par Examples
\verbatim
WlzContourCut -k plate.wlz iso.wlz > out.wlz
\endverbatim
Reads model contour from iso.wlz and knife contour from plate.wlz,
creates a new contour object in which all simplices of the model
are present unless they are intersected by simplices of the knife.
The cut contour object is then written to the output file out.wlz.
\par File
\ref WlzContourCut.c "WlzContourCut.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzGMModelCut "WlzGMModelCut(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok,
  		option,
		usage = 0;
  FILE		*fP = NULL;
  WlzDomain	oDom;
  WlzValues	oVal;
  char		*mFileStr,
  		*kFileStr,
		*oFileStr;
  WlzObject	*mObj = NULL,
		*kObj = NULL,
		*oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hk:o:";
  const char    fileStrDef[] = "-";

  oDom.core = NULL;
  oVal.core = NULL;
  /* Parse the argument list and check for input files. */
  opterr = 0;
  mFileStr = (char *)fileStrDef;
  kFileStr = (char *)fileStrDef;
  oFileStr = (char *)fileStrDef;
  while((option = getopt(argc, argv, optList)) != EOF)
  {
    switch(option)
    {
      case 'k':
	kFileStr = optarg;
	break;
      case 'o':
	oFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      mFileStr = argv[optind];
    }
  }
  /* Read the model contour object. */
  if(ok)
  {
    if((mFileStr == NULL) ||
       (*mFileStr == '\0') ||
       ((fP = (strcmp(mFileStr, "-")?
              fopen(mFileStr, "r"): stdin)) == NULL) ||
       ((mObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read model contour object from file %s.\n",
                     argv[0], mFileStr);
    }
    if(fP && strcmp(mFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read the knife object. */
  if(ok)
  {
    if((kFileStr == NULL) ||
       (*kFileStr == '\0') ||
       ((fP = (strcmp(kFileStr, "-")?
              fopen(kFileStr, "r"): stdin)) == NULL) ||
       ((kObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read knife contour object from file %s.\n",
                     argv[0], kFileStr);
    }
    if(fP && strcmp(kFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    oDom.ctr = WlzMakeContour(&errNum);
    switch(kObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	if(errNum == WLZ_ERR_NONE)
	{
	  oDom.ctr->model = WlzGMModelCutDom(mObj->domain.ctr->model,
					     kObj, &errNum);
	}
	break;
      case WLZ_CONTOUR:
	if(errNum == WLZ_ERR_NONE)
	{
	  oDom.ctr->model = WlzGMModelCut(mObj->domain.ctr->model,
					  kObj->domain.ctr->model, &errNum);
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      oObj = WlzMakeMain(WLZ_CONTOUR, oDom, oVal, NULL, NULL, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(oDom.ctr != NULL)
      {
	(void )WlzFreeContour(oDom.ctr);
      }
      else if(oDom.ctr->model != NULL)
      {
	(void )WlzGMModelFree(oDom.ctr->model);
      }
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to cut contour (%s).\n",
		     argv[0], errMsgStr);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(oFileStr, "-")?
	     fopen(oFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], oFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, oObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
                     argv[0], errMsgStr);
    }
  }
  if(fP && strcmp(oFileStr, "-"))
  {
      (void )fclose(fP);
  }
  (void )WlzFreeObj(mObj);
  (void )WlzFreeObj(kObj);
  (void )WlzFreeObj(oObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-k <knife>] [-o <out>] [-h] [<model>]\n"
    "Given model and knife contours; creates a new contour which consists\n"
    "of all simplices of the model that do not intersect simplices of the\n"
    "knife. This may be used to cut into a contour removing sections from\n"
    "it.\n"
    "Options are:\n"
    "  -k  The knife contour object.\n"
    "  -o  Output object file.\n"
    "  -h  Help - prints this usage message\n",
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
