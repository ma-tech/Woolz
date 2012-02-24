#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBasisFnTransformVertices_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzBasisFnTransformVertices.c
* \author       Jianguo Rao
* \date         January 2004
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
* \brief	Applies a basis function to vertices.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzbasisfntransformvertices "WlzBasisFnTransformVertices"
*/

/*!
\ingroup BinWlz
\defgroup wlzbasisfntransformvertices WlzBasisFnTransformVertices
\par Name
WlzBasisFnTransformVertices  -  applies a basis function to vertices.
\par Synopsis
\verbatim
WlzBasisFnTransformVertices [-o<out object>] [-p<tie points file>]
                            [-m<min mesh dist>] [-M<max mesh dist>]
			    [-t<basis fn transform>] [-Y<order of polynomial>]
			    [-g] [-h] [-q] [-Q] [-s] [-y] [-B] [-D] [-G]
			    [-L] [-T] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-o</b></td>
    <td>Output vertices file name.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Tie point file.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Vertex file.</td>
  </tr>
  <tr> 
    <td><b>-c</b></td>
    <td>Use conformal polynomial basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Use Gaussian basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-q</b></td>
    <td>Use multi-quadric basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-Q</b></td>
    <td>Use inverse-multi-quadric basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Use thin plate spline basis function (default) if tie points
        are given.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Use polynomianl basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-Y</b></td>
    <td>Polynomial order for oolynomianl basis function (default 3).</td>
  </tr>
</table>
\par Description
Computes and applies Woolz basis function transforms.
Tie points may be read from an ascii file with the format:
\verbatim
  <vertex x> <vertex y> <displacement x> <displacement y>
\endverbatim
\par Examples
\verbatim
WlzBasisFnTransformVertices -o outVertices.dat -p points.tie -v vertices.dat
\endverbatim
A thin plate spline basis function transform is computed from the
tie points read from the file points.tie. This transform is then
applied to the vertices read from vertices.dat.
The resulting vertices are then written to outVertices.dat.
\par File
\ref WlzBasisFnTransformVertices.c "WlzBasisFnTransformVertices.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzbasisfntransformobj "WlzBasisFnTransformObj(1)"
\ref WlzBasisFnTransformVertexI "WlzBasisFnTransformVertexI(3)"
\ref WlzBasisFnTransformVertexF "WlzBasisFnTransformVertexF(3)"
\ref WlzBasisFnTransformVertexD "WlzBasisFnTransformVertexD(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

#define	IN_RECORD_MAX   (1024)


extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		nTiePP,
  		nVertP,
		option,
		vxCount = 0,
		vxLimit = 0,
		basisFnPolyOrder = 3,
		ok = 1,
		ic = 0,
		usage = 0;
  WlzDVertex2	*vx0,
  		*vx1,
		*vxVec0  = NULL,
		*vxVec1  = NULL,
		*vertVec = NULL,
		*vertVecTr = NULL,
		srcV,
		transV;
  WlzMeshTransform *meshTr = NULL;
  WlzBasisFnTransform *basisTr = NULL;
  WlzFnType basisFnType = WLZ_FN_BASIS_2DMQ;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*rec,
  		*inObjFileStr,
		*tiePtFileStr  = NULL,
		*vertptFileStr = NULL,
  		*outVerticesFileStr;
  const char    *errMsg;
  static char	optList[] = "o:p:v:t:Y:cghqsy",
  		inObjFileStrDef[] = "-",
		outVerticesFileStrDef[] = "-",
  		inRecord[IN_RECORD_MAX];

  opterr = 0;
  outVerticesFileStr = outVerticesFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outVerticesFileStr = optarg;
	break;
      case 'c':
        basisFnType = WLZ_FN_BASIS_2DCONF_POLY;
	break;
      case 'g':
        basisFnType = WLZ_FN_BASIS_2DGAUSS;
	break;
      case 'p':
        tiePtFileStr  = optarg;
	break;
      case 'v':
        vertptFileStr = optarg;
	break;
      case 'q':
        basisFnType = WLZ_FN_BASIS_2DMQ;
	break;
      case 'Q':
        basisFnType = WLZ_FN_BASIS_2DIMQ;
	break;
      case 's':
        basisFnType = WLZ_FN_BASIS_2DTPS;
	break;
      case 'y':
        basisFnType = WLZ_FN_BASIS_2DPOLY;
	break;
      case 'Y':
        if(sscanf(optarg, "%d", &basisFnPolyOrder) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if( (inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     ((tiePtFileStr == NULL) || (*tiePtFileStr == '\0')) || 
     (outVerticesFileStr == NULL) || (*outVerticesFileStr == '\0'))
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
  if(ok)
  {
    if(tiePtFileStr)
    {
      if((fP = (strcmp(tiePtFileStr, "-")?
                fopen(tiePtFileStr, "r"): stdin)) == NULL)
      {
        ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open tie points file %s (%s).\n",
		       *argv, tiePtFileStr, errMsg);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	while((errNum == WLZ_ERR_NONE) &&
	      (fgets(inRecord, IN_RECORD_MAX - 1, fP) != NULL))
	{
	  inRecord[IN_RECORD_MAX - 1] = '\0';
	  rec = inRecord;
	  while(*rec && isspace(*rec))
	  {
	    ++rec;
	  }
	  if(*rec && (*rec != '#'))
	  {
	    if(vxCount >= vxLimit)
	    {
	      vxLimit = (vxLimit + 1024) * 2;
	      if(((vxVec0 = (WlzDVertex2 *)AlcRealloc(vxVec0,
				     vxLimit * sizeof(WlzDVertex2))) == NULL) ||
		 ((vxVec1 = (WlzDVertex2 *)AlcRealloc(vxVec1,
				     vxLimit * sizeof(WlzDVertex2))) == NULL))
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		vx0 = vxVec0 + vxCount;
		vx1 = vxVec1 + vxCount;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(sscanf(rec, "%lg %lg %lg %lg", &(vx0->vtX), &(vx0->vtY),
			&(vx1->vtX), &(vx1->vtY)) != 4)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		vx1->vtX += vx0->vtX;
		vx1->vtY += vx0->vtY;
		++vx0;
		++vx1;
		++vxCount;
	      }
	    }
	  }
	}
	nTiePP = vxCount;
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
			 "%s: failed to read tie points file %s (%s).\n",
			 *argv, tiePtFileStr, errMsg);
	}
      }
      if(fP && strcmp(tiePtFileStr, "-"))
      {
	fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    /* get the basis function transform */
      basisTr = WlzBasisFnTrFromCPts2D(basisFnType, basisFnPolyOrder,
					nTiePP, vxVec0,
					nTiePP, vxVec1, NULL, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		     "%s: failed to compute basis function transform (%s).\n",
		     *argv, errMsg);
      }
  }
  if(ok)
  {
    /* get the vertices to be transformed */
    vxCount = 0;
    vxLimit = 0;
    if(vertptFileStr)
    {
      if((fP = (strcmp(vertptFileStr, "-")?
                fopen(vertptFileStr, "r"): stdin)) == NULL)
      {
        ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open vertices file %s (%s).\n",
		       *argv, vertptFileStr, errMsg);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	while((errNum == WLZ_ERR_NONE) &&
	      (fgets(inRecord, IN_RECORD_MAX - 1, fP) != NULL))
	{
	  inRecord[IN_RECORD_MAX - 1] = '\0';
	  rec = inRecord;
	  while(*rec && isspace(*rec))
	  {
	    ++rec;
	  }
	  if(*rec && (*rec != '#'))
	  {
	    if(vxCount >= vxLimit)
	    {
	      vxLimit = (vxLimit + 1024) * 2;
	      if(((vertVec = (WlzDVertex2 *)AlcRealloc(vertVec,
				     vxLimit * sizeof(WlzDVertex2))) == NULL)  ||
		 ((vertVecTr = (WlzDVertex2 *)AlcRealloc(vertVecTr,
				     vxLimit * sizeof(WlzDVertex2))) == NULL))

	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		vx0 = vertVec + vxCount;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(sscanf(rec, "%lg %lg", &(vx0->vtX), &(vx0->vtY)
		       ) != 2)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		++vx0;
		++vxCount;
	      }
	    }
	  }
	}
	nVertP = vxCount;
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
			 "%s: failed to read vertices file %s (%s).\n",
			 *argv, tiePtFileStr, errMsg);
	}
      }
      if(fP && strcmp(tiePtFileStr, "-"))
      {
	fclose(fP);
      }
      fP = NULL;
    }

  }
  if(ok)
  {
    /* Transform the vertices and then write it out. */
    while(ic < nVertP)
    {
	if(errNum == WLZ_ERR_NONE)
	{
          srcV.vtX  = (vertVec + ic )->vtX;
	  srcV.vtY  = (vertVec + ic )->vtY;
	  
	  transV = WlzBasisFnTransformVertexD( basisTr,
				 	   srcV,
					   &errNum);
	  (vertVecTr + ic )->vtX = transV.vtX;
	  (vertVecTr + ic )->vtY = transV.vtY;

	  if(errNum != WLZ_ERR_NONE)
	  {
	    ok = 0;
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
			   "%s: failed to transfer the vertices  (%s).\n",
			   *argv, errMsg);
	  }
	  /*
	  vertVecTr++;
	  vertVec++;
	  */
	}
	ic++;
    }
    /* output the transformed vertices */
    if(( fP = fopen(outVerticesFileStr, "w")) == NULL )
    {
       printf("cannot open output file.\n");
       exit(1);
    }
    ic = 0; 
    while(ic < nVertP)
    {
	(void )fprintf(fP,
		       "%g %g\n",
		       ( vertVecTr + ic )->vtX, ( vertVecTr + ic )->vtY
		       );
		       ic++;
		       /*
		       vertVecTr++;
		       */
    }
    fclose(fP);
    fP = NULL;



  }
  (void )WlzMeshFreeTransform(meshTr);
  (void )WlzBasisFnFreeTransform(basisTr);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out object>] [-p<tie points file>]\n"
    "                  [-m<min mesh dist>] [-M<max mesh dist>]\n"
    "                  [-t<basis fn transform>] [-Y<order of polynomial>]\n"
    "                  [-g] [-h] [-q] [-Q] [-s] [-y] [-B] [-D] [-G] [-L]\n"
    "                  [-T] [<in object>]\n"
    "Options:\n"
    "  -o  Output vertices file name.\n"
    "  -p  Tie point file.\n"
    "  -v  vertices file.\n"
    "  -c  Use conformal polynomial basis function if tie points are given.\n"
    "  -g  Use Gaussian basis function if tie points are given.\n"
    "  -h  Help, prints this usage message.\n"
    "  -q  Use multi-quadric basis function if tie points are given.\n"
    "  -Q  Use inverse-multi-quadric basis function if tie points are given.\n"
    "  -s  Use thin plate spline basis function (default) if tie points\n"
    "      are given.\n"
    "  -y  Use polynomianl basis function if tie points are given.\n"
    "  -Y  Polynomial order for oolynomianl basis function (default 3).\n"
    "Computes and applies Woolz basis function transforms.\n"
    "Tie points may be read from an ascii file with the format:\n"
    "  <vertex x> <vertex y> <displacement x> <displacement y>\n"
    "vertices are read from the file\n",
    *argv,
    " -o outVertices.dat -p points.tie -v vertices.dat\n"
    "A thin plate spline basis function transform is computed from the\n"
    "tie points read from the file points.tie. This transform is then\n"
    "applied to the vertices read from vertices.dat. The resulting vertices\n"
    "is then written to outVertices.dat.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
