#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTstVxInSimplex_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTstVxInSimplex.c
* \author       Bill Hill
* \date         July 2007
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2007 Medical research Council, UK.
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
* \brief	Test program for WlzGeomVxInTriangle() and
* 		WlzGeomVxInTetrahedron().
* \ingroup	WlzTst
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <Wlz.h>

static int			WlzTstInGetPosition(
				  int nVx,
				  WlzDVertex3 *vx,
                                  int dim,
				  char *pos[]);
static char			*WlzTstInValueToStr(
				  int inV);

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		inside,
  		option,
  		dim = 2,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  WlzDVertex3	vP;
  WlzDVertex2	vP2;
  WlzDVertex3	v[4];
  WlzDVertex2	v2[3];
  static char   optList[] = "23hvp:";


  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case '2':
        dim = 2;
	break;
      case '3':
        dim = 3;
	break;
      case 'p':
	if(WlzTstInGetPosition(1, &vP, dim, &optarg) != 1)
	{
	  usage = 1;
	}
        break;
      case 'v':
        verbose = 1;
	break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if(((optind + dim + 1) != argc) ||
       (WlzTstInGetPosition(dim + 1, &(v[0]), dim, argv + optind) != dim + 1))
    {
      usage = 1;
      ok = 0;
    }
  }
  if(ok)
  {
    if(dim == 2)
    {
      v2[0].vtX = v[0].vtX; v2[0].vtY = v[0].vtY;
      v2[1].vtX = v[1].vtX; v2[1].vtY = v[1].vtY;
      v2[2].vtX = v[2].vtX; v2[2].vtY = v[2].vtY;
      vP2.vtX = vP.vtX; vP2.vtY = vP.vtY;
      inside = WlzGeomVxInTriangle(v2[0], v2[1], v2[2], vP2);
      if(verbose)
      {
	(void )printf("(%g,%g) is %s triangle (%g,%g), "
	              "(%g,%g), (%g,%g).\n",
		       vP2.vtX, vP2.vtY,
		       WlzTstInValueToStr(inside),
		       v2[0].vtX, v2[0].vtY,
		       v2[1].vtX, v2[1].vtY,
		       v2[2].vtX, v2[2].vtY);

      }
      else
      {
        (void )printf("%d", inside);
      }
    }
    else /* dim == 3 */
    {
      inside = WlzGeomVxInTetrahedron(v[0], v[1], v[2], v[3],
                                      vP);
      if(verbose)
      {
	(void )printf("(%g,%g,%g) is %s tetrahedron (%g,%g,%g), "
	               "(%g,%g,%g), (%g,%g,%g), (%g,%g,%g).\n",
		       vP.vtX, vP.vtY, vP.vtZ,
		       WlzTstInValueToStr(inside),
		       v[0].vtX, v[0].vtY, v[0].vtZ,
		       v[1].vtX, v[1].vtY, v[1].vtZ,
		       v[2].vtX, v[2].vtY, v[2].vtZ,
		       v[3].vtX, v[3].vtY, v[3].vtZ);

      }
      else
      {
        (void )printf("%d", inside);
      }
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-2|3] [-p <query vertex position>] [-v]\n"
    "                         <vertex position list>\n"
    "Options are:\n"
    " -2  Simplex is a triangle.\n"
    " -3  simplex is a tetrahedron.\n"
    " -p  Query position.\n"
    " -v  Verbose output.\n"
    "Each psoition must be given as a comma seperated list of floating point\n"
    "point values and the vertices of the vertex position list must be\n"
    "white space seperated. White spaces are not allowed within a psoition.\n"
    "Example:\n"
    "  %s -2 -p 3,4 0,0 10,0 10,10\n"
    "determines whether the point 3,4 (x = 3, y = 4) is within the triangle\n"
    "having vertices (0,0), (10,0), (10,10).\n",
    argv[0], argv[0]);
  }
  exit(!ok);
}

static int	WlzTstInGetPosition(int nVx, WlzDVertex3 *vx,
                                    int dim, char *pos[])
{
  int		ok,
  		idx;

  ok = 1;
  idx = 0;
  while(ok && (idx < nVx))
  {
    if(dim == 2)
    {
      ok = sscanf(pos[idx], "%lg,%lg",
                  &((vx + idx)->vtX), &((vx + idx)->vtY)) == 2;
    }
    else /* dim == 3 */
    {
      ok = sscanf(pos[idx], "%lg,%lg,%lg",
                  &((vx + idx)->vtX), &((vx + idx)->vtY),
		  &((vx + idx)->vtZ)) == 3;
    }
    if(ok)
    {
      ++idx;
    }
  }
  return(idx);
}

static char	*WlzTstInValueToStr(int inV)
{
  char		*str;
  static char	*inStr = "inside",
  		*onStr = "on",
		*outStr = "outside";
 
 if(inV < 0)
 {
   str = outStr;
 }
 else if(inV > 0)
 {
   str = inStr;
 }
 else
 {
   str = onStr;
 }
 return(str);
}
