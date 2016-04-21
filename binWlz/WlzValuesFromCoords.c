#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzValuesFromCoords_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzValuesFromCoords.c
* \author       Bill Hill
* \date         April 2016
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
* \brief	Adds integer coordinate values to the domain of the
* 		input object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzvaluesfromcoords "WlzValuesFromCoords"
*/


/*!
\ingroup BinWlz
\defgroup wlzvaluesfromcoords  WlzValuesFromCoords
\par Name
WlzValuesFromCoords - creates an object in which the pixel/voxel values are
		      those of their coordinates.
\par Synopsis
\verbatim
WlzValuesFromCoords [-h] [-S] [-o<output file>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name, default standard output.</td>
  </tr>
  <tr> 
    <td><b>-S</b></td>
    <td>Split the comound object into seperate files, appending
        X, Y or Z to the body of the output file name.</td>
  </tr>
</table>
\par Description
Reads a 2D or 3D spatial domain object and computes a compound
object in which each component object has the domain of the input
object and the values of it's pixel/voxel coordinates. When a
compound objec is output, the components are the column (x)
coordinates followed by line (Y) and (for 3D objects), plane
(Z). By default the input object is read from the standard input
and the output is writen to the standard output.
\par File
\ref binWlz/WlzValuesFromCoords.c "WlzValuesFromCoords.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzValuesFromCoords "WlzValuesFromCoords(1)"
\ref WlzValuesFromCoords "WlzValuesFromCoords(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <Wlz.h>

#define WLZ_CFP_READLN_LEN	(1024)

static WlzVertexP 		WlzMTDReadVtxArray(
				  FILE *fP,
				  int dim,
				  int *dstNVtx, 
				  WlzVertexType *dstVType,
				  WlzErrorNum *dstErr);

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
		split = 0,
  		usage = 0;
  FILE		*fP = NULL;
  char		*inFile,
  		*outFile;
  const char	*errMsg;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "hSo:";
  const char    fileDef[] = "-";

  opterr = 0;
  inFile = (char *)fileDef;
  outFile = (char *)fileDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outFile = optarg;
	break;
      case 'S':
        split = 1;
	break;
      case 'h': /* FALLTROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(!usage)
  {
    if((outFile == NULL) || (*outFile == '\0') ||
       (inFile == NULL) || (*inFile == '\0'))
    {
      usage = 1;
    }
  }
  if(!usage && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      inFile = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    if((inFile == NULL) ||
       (*inFile == '\0') ||
       ((fP = (strcmp(inFile, "-")?
              fopen(inFile, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_READ_EOF;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to read input file %s (%s).\n",
		     *argv, inFile, errMsg);
    }
  }
  if(fP)
  {
    if(strcmp(inFile, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  if(ok)
  {
    outObj = WlzValuesFromCoords(inObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to create coordinate valued object (%s).\n",
      		     argv[0],
		     errMsg);
    }
  }
  if(ok)
  {
    int		idx;
    WlzCompoundArray *cpd;

    cpd = (WlzCompoundArray *)outObj;
    if(split && strcmp(outFile, "-"))
    {
      size_t	len;
      char 	*dot,
	      	*splitFile = NULL;

      len = strlen(outFile);
      if((splitFile = (char *)AlcCalloc(len + 6, sizeof(char))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	dot = strrchr(outFile, '.');
	if(dot)
	{
	  *dot = '\0';
	}
	for(idx = 0; idx < cpd->n; ++idx)
	{
	  int	last;
	  char  dim;

	  switch(idx)
	  {
	    case 0:
	      dim = 'X';
	      break;
	    case 1:
	      dim = 'Y';
	      break;
	    default:
	      dim = 'Z';
	      break;
	  }
	  last = sprintf(splitFile, "%s%c", outFile, dim);
	  if(dot && *(dot + 1))
	  {
	    (void )sprintf(splitFile + last, "%c%s", '.', dot + 1);
	  }
	  if((fP = fopen(splitFile, "w")) == NULL)
	  {
	    errNum = WLZ_ERR_FILE_OPEN;
	  }
	  else
	  {
	    errNum = WlzWriteObj(fP, cpd->o[idx]);
	    (void )fclose(fP);
	  }
	}
	AlcFree(splitFile);
      }
    }
    else
    {
      if((fP = (strcmp(outFile, "-")? fopen(outFile, "w"): stdout)) == NULL)
      {
        errNum = WLZ_ERR_FILE_OPEN;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(split)
	{
	  for(idx = 0; idx < cpd->n; ++idx)
	  {
	    if((errNum = WlzWriteObj(fP, cpd->o[idx])) != WLZ_ERR_NONE)
	    {
	      break;
	    }
	  }
	}
	else
	{
          errNum = WlzWriteObj(fP, outObj);
	}
      }
      if(fP && strcmp(outFile, "-"))
      {
	(void )fclose(fP);
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsg);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-S] [-o<output file>] [<input file>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -S  Split the comound object into seperate files, appending\n"
    "      X, Y or Z to the body of the output file name.\n"
    "Reads a 2D or 3D spatial domain object and computes a compound\n"
    "object in which each component object has the domain of the input\n"
    "object and the values of it's pixel/voxel coordinates. When a\n"
    "compound objec is output, the components are the column (x)\n"
    "coordinates followed by line (Y) and (for 3D objects), plane\n"
    "(Z). By default the input object is read from the standard input\n"
    "and the output is writen to the standard output.\n",
    argv[0],
    WlzVersion());
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
