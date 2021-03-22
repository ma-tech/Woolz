#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTransformProduct_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzTransformProduct.c
* \author       Bill Hill
* \date         August 2003
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
* \brief	Computes the product of a pair of transforms.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlztransformproduct "WlzTransformProduct"
*/

/*!
\ingroup BinWlz
\defgroup wlztransformproduct WlzTransformProduct
\par Name
WlzTransformProduct - computes the product of a pair of transforms.
\par Synopsis
\verbatim
WlzTransformProduct [-h] [-o<output object>] [transform 0] [transform 1]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
</table>
\par Description
Computes the product of a pair of transforms.
\par Examples
\verbatim
WlzTransformProduct -o t2.wlz t0.wlz t1.wlz
\endverbatim
Computes the product of \f$\mathbf{T_0}\f$ and \f$\mathbf{T_1}\f$,
such that \f$\mathbf{T_2}\f$ = \f$\mathbf{T_1} \mathbf{T_0}\f$,
i.e. transform \f$\mathbf{T_0}\f$ is applied first then \f$\mathbf{T_1}\f$.
\par File
\ref WlzTransformProduct.c "WlzTransformProduct.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzaffinetransformobj "WlzAffineTransformObj(1)"
\ref WlzTransformProduct "WlzTransformProduct(3)"
\ref WlzAffineTransformProduct "WlzAffineTransformProduct(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		usage = 0;
  WlzTransform  t[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outObjFileStr = NULL;
  char		*inObjFileStr[2];
  const char	*errMsgStr;
  static char	optList[] = "ho:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outObjFileStr = outFileStrDef;
  t[0].core = t[1].core = t[2].core = NULL;
  inObjFileStr[0] = inObjFileStr[1] = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if(optind < argc)
    {
      if((optind + 1) == argc)
      {
	inObjFileStr[0] = *(argv + optind);
      }
      else if((optind + 2) == argc)
      {
	inObjFileStr[0] = *(argv + optind);
	inObjFileStr[1] = *(argv + optind + 1);
      }
      else
      {
	usage = 1;
	ok = 0;
      }
    }
  }
  if(ok)
  {
    int		i;
    WlzObject	*obj = NULL;

    for(i = 0; i < 2; ++i)
    {
      if((*inObjFileStr[i] == '\0') ||
	 ((fP = (strcmp(inObjFileStr[i], "-")?
		fopen(inObjFileStr[i], "r"): stdin)) == NULL) ||
	 ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	 (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
      }
      if(fP && strcmp(inObjFileStr[i], "-"))
      {
	(void )fclose(fP);
      }
      if(ok)
      {
        switch(obj->type)
	{
	  case WLZ_EMPTY_OBJ:
	    t[i].empty = WlzMakeEmptyTransform(&errNum);
	    (void )WlzAssignTransform(t[i], NULL);
	    break;
	  case WLZ_TRANS_OBJ:             /* FALLTHROUGH */
	  case WLZ_AFFINE_TRANS:
	    t[i].affine = obj->domain.t;
	    (void )WlzAssignTransform(t[i], NULL);
	    break;
	  case WLZ_MESH_TRANS:
	    t[i].mesh = obj->domain.mt;
	    (void )WlzAssignTransform(t[i], NULL);
	    break;
	  case WLZ_TRANSFORM_2D5_MESH:    /* FALLTHROUGH */
	  case WLZ_TRANSFORM_3D_MESH:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	  case WLZ_CMESH_2D:              /* FALLTHROUGH */
	  case WLZ_CMESH_2D5:             /* FALLTHROUGH */
	  case WLZ_CMESH_3D:
	    t[i].obj = WlzAssignObject(obj, NULL);
	    break;
	  default:
	    errNum = WLZ_ERR_TRANSFORM_TYPE;
	    break;

	}
      }
      (void )WlzFreeObj(obj);
      if(!ok)
      {
	(void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		       "%s: Failed to read transform from file %s (%s).\n",
		       *argv, inObjFileStr[i], errMsgStr);
        break;
      }
    }
  }
  if(ok)
  {
    t[3] = WlzTransformProduct(t[0], t[1], &errNum);
    if((t[3].core == NULL) || (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s: Failed to compute product (%s).\n",
                     *argv, errMsgStr);
    }
  }
  if(ok)
  {
    WlzObject	*obj = NULL;
    WlzDomain	dom;
    WlzValues	val;

    switch(t[3].core->type)
    {
      case WLZ_TRANSFORM_EMPTY:
        obj = WlzAssignObject(
	      WlzMakeEmpty(&errNum),
	      NULL);
        break;
      case WLZ_TRANSFORM_2D_REG:  	/* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_TRANS:      /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_NOSHEAR:    /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_AFFINE:     /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_REG:        /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_TRANS:      /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_NOSHEAR:    /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_AFFINE:
	dom.t = t[3].affine;
	val.core = NULL;
	obj = WlzAssignObject(
	      WlzMakeMain(WLZ_AFFINE_TRANS, dom, val, NULL, NULL, &errNum),
	      NULL);
        break;
      case WLZ_TRANSFORM_2D_MESH:
	dom.mt = t[3].mesh;
	val.core = NULL;
	obj = WlzAssignObject(
	      WlzMakeMain(WLZ_MESH_TRANS, dom, val, NULL, NULL, &errNum),
	      NULL);
        break;
      case WLZ_TRANSFORM_2D_CMESH:	/* FALLTHROUGH */
      case WLZ_TRANSFORM_2D5_CMESH:	/* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_CMESH:
	obj = WlzAssignObject(t[3].obj, NULL);
        break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		     "%s: Failed to make transform object (%s).\n",
		     *argv, errMsgStr);
    }
    else
    {
      if(((fP = (strcmp(outObjFileStr, "-")?
		fopen(outObjFileStr, "w"):
		stdout)) == NULL) ||
	  (WlzWriteObj(fP, obj) != WLZ_ERR_NONE))
      {
	ok = 0;
        (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		       "%s: Failed to write transform product object (%s).\n",
		       *argv, errMsgStr);
      }
      if(fP && strcmp(outObjFileStr, "-"))
      {
	fclose(fP);
      }
      (void )WlzFreeObj(obj);
    }
  }
  (void )WlzFreeTransform(t[0]);
  (void )WlzFreeTransform(t[1]);
  (void )WlzFreeTransform(t[2]);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-h] [-o<output object>] [transform 0] [transform 1]\n" 
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file name.\n"
    "Computes the product of a pair of transforms.\n",
    *argv,
    " -o t2.wlz t0.wlz t1.wlz\n"
    "Computes the product of t0 and t1 such that t2 = t1 t0 i.e.\n"
    "t0 is applied first then t1\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
