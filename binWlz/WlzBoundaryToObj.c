#pragma ident "MRC HGU $Id$"
/*!
* \file		binWlz/WlzBoundaryToObj.c
* \author       Richard Baldock
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
* \brief 	Converts boundary lists to domain objects.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzboundarytoobj "WlzBoundaryToObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzboundarytoobj WlzBoundaryToObj
\par Name
WlzBoundaryToObj  -  converts boundary lists to domain objects.
       objects.
\par Synopsis
\verbatim
WlzBoundaryToObj [-f s|e|o] [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-f</b></td>
    <td>Sets the boundary fill method to be one of: simple (s),
        even/odd (e) or vertex (v) fill, default is even/odd..</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Converts a boundary list woolz object to a corresponding domain
using the given fill method.
\par Examples
\verbatim
WlzBoundaryToObj bound.wlz >dom.wlz
\endverbatim
Reads a boundary list from bound.wlz and computes a domain objects,
which is then written to dom.wlz.
\par File
\ref WlzBoundaryToObj.c "WlzBoundaryToObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzobjtoboundary "WlzObjToBoundary(1)"
\ref WlzBoundaryToObj "WlzBoundaryToObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void ShowUsage(char *arg)
{
  fprintf(stderr,
	  "Usage: %s [-f s|e|o] [-h] [<input file>]\n"
	  "Options:\n"
	  "  -f Sets the boundary fill method to be one of: simple (s),\n"
	  "     even/odd (e) or vertex (v) fill, default is even/odd.\n"
	  "  -h Help - prints this usage message\n"
	  "Converts a boundary list woolz object to a corresponding domain\n"
	  "using the given fill method.\n",
	  arg);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "f:h";
  int		option,
  		usage = 0;
  WlzPolyFillMode fill = WLZ_EVEN_ODD_FILL;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  if((option = getopt(argc, argv, optList)) != EOF)
  {
    switch(option)
    {
      case 'f':
	switch(*optarg)
	{
	  case 'e':
	    fill = WLZ_EVEN_ODD_FILL;
	    break;
	  case 's':
	    fill = WLZ_SIMPLE_FILL;
	    break;
	  case 'v':
	    fill = WLZ_VERTEX_FILL;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }

  if(usage)
  {
    ShowUsage(argv[0]);
    return 1;
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      ShowUsage(argv[0]);
      return 1;
    }
  }

  /* Read objects and convert. */
  while((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) 
  {
    switch( obj->type )
    {
    case WLZ_BOUNDLIST:
      if((nobj = WlzBoundaryToObj(obj, fill, NULL)) != NULL)
      {
	(void )WlzWriteObj(stdout, nobj);
	WlzFreeObj(nobj);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if(obj->domain.p->type == WLZ_PLANEDOMAIN_BOUNDLIST)
      {
	if((nobj = WlzBoundaryToObj(obj, fill, NULL)) != NULL)
	{
	  (void )WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
      }
      else {
	(void )WlzWriteObj(stdout, obj);
      }
      break;
	
    default:
      (void )WlzWriteObj(stdout, obj);
      break;
    }

    WlzFreeObj(obj);
  }
  return WLZ_ERR_NONE;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
