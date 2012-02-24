#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstGeomLSqOPlane_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstGeomLSqOPlane.c
* \author       Bill Hill
* \date         June 2008
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
* \brief	Test for WlzGeometryLSqOPlane().
* \ingroup	BinWlzTst
*/

#include <string.h>
#include <float.h>
#include <stdio.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

static int 			WlzTstReadVtxList(
				  WlzDVertex3 **vtx,
				  FILE *fP);

int		main(int argc, char *argv[])
{
  int		option,
  		ok = 1,
		nVtx = 0,
  		usage = 0,
		verbose = 0;
  WlzDVertex3	c,
  		n;
  char		*inFileStr,
  		*outFileStr;
  FILE          *fP = NULL;
  WlzDVertex3	*vtx = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char   optList[] = "ho:v";
  const char	inFileStrDef[] = "-",
  		outFileStrDef[] = "-";

  inFileStr = (char *)inFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'v':
        verbose = 1;
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
  if(usage == 0)
  {
    if((inFileStr == NULL) || (*inFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
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
        inFileStr = *(argv + optind);
      }
    }
  }
  ok = usage == 0;
  if(ok)
  {
    fP = NULL;
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_READ_EOF;
      (void )fprintf(stderr,
                     "%s: Failed to open vertex file (%s)\n",
                     *argv, inFileStr);
    }
    if(ok)
    {
      nVtx = WlzTstReadVtxList(&vtx, fP);
      if(nVtx < 1)
      {
        ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )fprintf(stderr,
	               "%s: Failed to read vertices from file (%s)\n",
		       *argv, inFileStr);
      }
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    errNum = WlzGeometryLSqOPlane(&n, &c, nVtx, vtx);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to compute least squares plane (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    fP = NULL;
    if((outFileStr == NULL) ||
       (*outFileStr == '\0') ||
       ((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_WRITE_EOF;
      (void )fprintf(stderr,
                     "%s: Failed to open output file (%s)\n",
                     *argv, outFileStr);
    }
    if(ok)
    {
      if(verbose)
      {
	(void )fprintf(fP,
		       "normal = %lg %lg %lg\n"
		       "centroid = %lg %lg %lg\n",
		       n.vtX, n.vtY, n.vtZ,
		       c.vtX, c.vtY, c.vtZ);
      }
      else
      {
	(void )fprintf(fP,
		       "%lg %lg %lg\n"
		       "%lg %lg %lg\n",
		       n.vtX, n.vtY, n.vtZ,
		       c.vtX, c.vtY, c.vtZ);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
    
  }
  AlcFree(vtx);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<file>] [-v] [<file>]\n"
    "Options are:\n"
    "  -h  Outputs this usage message.\n"
    "  -o  Output file.\n"
    "  -v  Verbose output.\n"
    "Computes the plane with least square orthogonal distance to the given\n"
    "3D vertices.\n",
    argv[0]);
  }
  exit(!ok);
}

/*!
* \return	Number of vertices read.
* \ingroup	binWlzTst
* \brief	Reads input vertices from file in the format:
*                <x> <y> <z>.
* \param	vtx			Destination pointer for vertices.
* \param	fP			Input file pointer.
*/
static int 	WlzTstReadVtxList(WlzDVertex3 **vtx, FILE *fP)
{
  int		ok = 1,
  		inR = 0,
  		inC = 0,
  		nVtx = 0;
  double        **inData = NULL;

  if((AlcDouble2ReadAsci(fP, &inData,
                         (size_t *)&inR, (size_t *)&inC) != ALC_ER_NONE) ||
     (inC != 3) ||
     (inR < 1) ||
     ((*vtx = AlcMalloc(inR * sizeof(WlzDVertex3))) == NULL))
  {
    ok = 0;
  }
  if(ok)
  {
    int		idx;

    nVtx = inR;
    for(idx = 0; idx < nVtx; ++idx)
    {
      (*vtx + idx)->vtX = inData[idx][0];
      (*vtx + idx)->vtY = inData[idx][1];
      (*vtx + idx)->vtZ = inData[idx][2];
    }
  }
  AlcDouble2Free(inData);
  return(nVtx);
}
