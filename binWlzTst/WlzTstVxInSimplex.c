#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstVxInSimplex_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstVxInSimplex.c
* \author       Bill Hill
* \date         July 2007
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
* \brief	Test program for WlzGeomVxInTriangle() and
* 		WlzGeomVxInTetrahedron().
* \ingroup	BinWlzTst
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
		tri = 0,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  double	eps = 1.0e-6;
  WlzDVertex3	vP;
  WlzDVertex2	vP2;
  WlzDVertex3	v[4];
  WlzDVertex2	v2[3];
  static char   optList[] = "23htvp:";


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
      case 't':
        tri = 1;
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
    int	nv;

    if(dim == 2)
    {
      tri = 1;
    }
    nv = (tri != 0)? 3: 4;
    if(((optind + nv) != argc) ||
       (WlzTstInGetPosition(nv, &(v[0]), dim, argv + optind) != nv))
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
      inside = WlzGeomVxInTriangle2D(v2[0], v2[1], v2[2], vP2);
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
      if(tri == 0)
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
      else
      {
	inside = WlzGeomVxInTriangle3D(v[0], v[1], v[2], vP, eps);
	if(verbose)
	{
	  (void )printf("(%g,%g,%g) is %s triangle (%g,%g,%g), "
			 "(%g,%g,%g), (%g,%g,%g).\n",
			 vP.vtX, vP.vtY, vP.vtZ,
			 WlzTstInValueToStr(inside),
			 v[0].vtX, v[0].vtY, v[0].vtZ,
			 v[1].vtX, v[1].vtY, v[1].vtZ,
			 v[2].vtX, v[2].vtY, v[2].vtZ);

	}
	else
	{
	  (void )printf("%d", inside);
	}
      }
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-2|3] [-p <query vertex position>] [-v]\n"
    "                         <vertex position list>\n"
    "Options are:\n"
    " -2  Simplex is in 2D, this implies that the simplex is a triangle.\n"
    " -3  Simplex is in 3D and may either be a triangle or tetrahedron.\n"
    " -p  Query position.\n"
    " -t  Forces the simplex to be a triangle (rather than the default\n"
    "     tetrahedron) if dimension is 3.\n"
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
