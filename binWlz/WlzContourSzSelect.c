#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzContourSzSelect_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzContourSzSelect.c
* \author       Bill Hill
* \date         March 2001
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
* \brief	Filter to remove small shells from contours.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcontourszselect "WlzContourSzSelect"
*/

/*!
\ingroup BinWlz
\defgroup wlzcontourszselect WlzContourSzSelect
\par Name
WlzContourSzSelect  -  Filter remove small shells from contours.
\par Synopsis
\verbatim
WlzContourSzSelect  [-h] -s<threshold> [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Threshold shell size (edges in 2D, loops in 3D).</td>
  </tr>
</table>
\par Description
Removes all shells in the contour which have less than the threshold
number of elements.
\par Examples
\verbatim
ContourSzSelect -s 1000 lobster.wlz >out.wlz
\endverbatim
Removes all fragments with less than 1000 faces from  the 3D contour model
read from lobster.wlz and writes the output to out.wlz.
\par File
\ref WlzContourSzSelect.c "WlzContourSzSelect.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzcontourobj "WlzContourObj(1)"
\ref WlzGMFilterRmSmShells "WlzGMFilterRmSmShells(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <Wlz.h>
#include <string.h>


extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           option,
  		ok = 1,
		usage = 0,
		smShellSz = -1;
  FILE		*fP = NULL;
  char		*inObjFileStr,
		*outObjFileStr;
  WlzObject	*obj;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char    *errMsgStr;
  static char	optList[] = "ho:s:";
  const char	outObjFileStrDef[] = "-",
		inObjFileStrDef[] = "-";

  obj = NULL;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
	outObjFileStr = optarg;
	break;
      case 's':
        if(sscanf(optarg, "%d", &smShellSz) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(smShellSz <= 1)
  {
    ok = 0;
    usage = 1;
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0'))
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
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
	((fP = (strcmp(inObjFileStr, "-")?
		fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP);
      fP = NULL;
    }
  }
  if(ok)
  {
    if(obj->type != WLZ_CONTOUR)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Input object isn't a contour.\n",
		     argv[0]);
    }
  }
  if(ok)
  {
    errNum = WlzGMFilterRmSmShells(obj->domain.ctr->model, smShellSz);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		     "%s: failed to filter small shells from contour (%s).\n",
		     *argv, errMsgStr);
    }
  }
  if(ok)
  {
  if((fP = (strcmp(outObjFileStr, "-")?
           fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outObjFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object (%s).\n",
                     argv[0], errMsgStr);
    }
  }
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
  }
  if(obj)
  {
    (void )WlzFreeObj(obj);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s%s %s",
      *argv,
      " [-h] -s<threshold> [<input object>]\n"
      "Version: ",
      WlzVersion(),
      "\n"
      "Options:\n"
      "  -s  Threshold shell size (edges in 2D, loops in 3D).\n"
      "  -h  Prints this usage information.\n"
      "Removes all shells in the contour which have less than the threshold\n"
      "number of elements.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
