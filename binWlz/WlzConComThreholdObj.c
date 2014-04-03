#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConComThreholdObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzConComThreholdObj.c
* \author       Bill Hill
* \date         March 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Connected component thresholding.
* \ingroup	BinWlz
* \ref		wlzconcomthreholdobj "WlzConComThreholdObj"
*/

/*!
\ingroup 	BinWlz
\defgroup	wlzconcomthreholdobj "WlzConComThreholdObj"
\par Name
WlzConComThreholdObj - computes  connected component threshold domain give
		       seed points.
\par Synopsis
\verbatim
WlzConComThreholdObj [-d 2|3] [-h] [-r <radius>] [-R l|m|h]
                     [-s<seed point file>] [-S<seed point list>]
                     [-H|L|E] [-x <extra>] [-o<output file>] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-d</b></td>
    <td>Seed point dimension, 2 or 3 for 2 or 3D seeds.</td>
  </tr>
  <tr>
    <td><b>-E</b></td>
    <td>Threshold equal, keep pixels/voxels with the threshold value
        at the seed points.</td>
  </tr>
  <tr>
    <td><b>-H</b></td>
    <td>Threshold high, keep pixels/voxels at or above the threshold
        value at the seed points (default).</td>
  </tr>
  <tr>
    <td><b>-L</b></td>
    <td>Threshold low, keep pixels/voxels at or below the threshold
        value at the seed points (watch out that this is different
        to WlzThreshold).</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Show this usage message.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name (default -, stdout).</td>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Radius of sampling region around seed points (default = 0,
        ie single pixel/voxel).</td>
  </tr>
  <tr>
    <td><b>-R</b></td>
    <td>Use the lowest (l), mean (m) or highest (h) value in the
        sampling region as the threhold value (default mean).</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>File containing seed point coordinates, one per line.</td>
  </tr>
  <tr>
    <td><b>-S</b></td>
    <td>Seed point coordinate list (comma separated).</td>
  </tr>
  <tr>
    <td><b>-x</b></td>
    <td>Integer percentage p to add/subtract to/from the grey value at
        each seed,  with v(s) * ((100 +/- p) / 100). The sign is
        choosen to keep the seed within the thresholded region also
        p is ignored for threshold equal (default p = 0).</td>
  </tr>
</table>
\par Description
Computes a connected component threshold domain for the given 2 or
3D image object, in which the threshold domain is the union of all
connected regions, thresholded using a value at a seed point, which
include that seed point.
A command line seed point list must be give the seed points as comma
separated values, while seed points in a file may comma or space
separated but must each be on a seperate line. The dimension of the
seed points is only needed for those specified using a command line
list. Valid syntax for the seed points in a seed point list is:
\verbatim
  -d 2 -S <x1>,<y1>,<x2>,<y2>...,<xn>,<yn>
\endverbatim
or
\verbatim
  -d 3 -S<x1>,<y1>,<z1><x2>,<y2>,<z2>...,<xn>,<yn>,<zn>
\endverbatim
Objects are read from stdin and written to stdout unless filenames
are given.
\par Examples
\verbatim
WlzConComThreholdObj -d 2 -S 3,6,102,203 -r 2 -R l -s seedsfile -o out.wlz
                     -L in.wlz
\endverbatim
The 2D seed points (3,6) and (102,203) are appended to the seed points
read from the file \"seedsfile\". These are then used to threshold
the 2D image read from the file "in.wlz" using the lowest value
in a radius 2 circle around each seed point for the theshold value
and forming the connected components using image values below these
values.
The output domain is then written to the file "out.wlz".
\par File
\ref WlzConComThreholdObj.c "WlzConComThreholdObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzConComThreshold(3) "WlzConComThreshold(3)"
\ref WlzLabel(3) "WlzLabel(3)"
\ref WlzThreshold(3) "WlzThreshold(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

#define WLZ_CCT_READLN_LEN      (1024)

static WlzVertexP		WlzCCTParseSeedList(
				  char *seedListStr,
				  int dim,
				  int *nSeeds,
				  int *usage);
static WlzVertexP		WlzCCTParseSeedFile(
				  FILE *fP,
				  WlzVertexP seeds,
				  int *dim,
				  int *nSeeds,
				  int *err);

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		dim = 0,
		usage = 0,
  		nSeeds = 0,
		xtra = 0;
  double	rad = 0.0;
  WlzVertexP	seeds;
  WlzObject	*iObj = NULL,
		*oObj = NULL;
  WlzThresholdType rHiLo = WLZ_THRESH_EQUAL,
  		   tHiLo = WLZ_THRESH_HIGH;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*seedFileStr = NULL,
  		*seedListStr = NULL,
		*iObjFileStr,
  		*oObjFileStr;
  const char    *errMsg;
  static char	optList[] = "hHLE:d:o:r:R:s:S:x:",
  		fileStrDef[] = "-";

  opterr = 0;
  seeds.v = NULL;
  oObjFileStr = fileStrDef;
  iObjFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'd':
        if((sscanf(optarg, "%d", &dim) != 1) ||
	   ((dim != 2) && (dim != 3)))
	{
	  usage = 1;
	}
	break;
      case 'H':
        tHiLo = WLZ_THRESH_HIGH;
	break;
      case 'L':
        tHiLo = WLZ_THRESH_LOW;
	break;
      case 'E':
        tHiLo = WLZ_THRESH_EQUAL;
	break;
      case 'o':
        oObjFileStr = optarg;
	break;
      case 'r':
        if(sscanf(optarg, "%lg", &rad) != 1)
	{
	  usage = 1;
	}
	break;
      case 'R':
        if(strlen(optarg) != 1)
	{
	  usage = 1;
	}
	else
	{
	  switch(*optarg)
	  {
	    case 'l':
	      rHiLo = WLZ_THRESH_LOW;
	      break;
	    case 'm':
	      rHiLo = WLZ_THRESH_EQUAL;
	      break;
	    case 'h':
	      rHiLo = WLZ_THRESH_HIGH;
	      break;
	    default:
	      usage = 1;
	      break;
	  }
	}
	break;
      case 's':
        seedFileStr = optarg;
        break;
      case 'S':
        seedListStr = optarg;
	break;
      case 'x':
        if((sscanf(optarg, "%d", &xtra) != 1) || (xtra < 0) || (xtra > 100))
	{
	  usage = 1;
	}
	break;
      case 'h':
      default:
        usage = 1;
	break;
    }
  }
  if((usage == 0) && (seedListStr != NULL))
  {
    if(dim == 0)
    {
      usage = 1;
      errNum = WLZ_ERR_PARAM_DATA;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	  "%s: Dimension required for command line seed list (%s).\n",
	  *argv, errMsg);
    }
    else
    {
      seeds = WlzCCTParseSeedList(seedListStr, dim, &nSeeds, &usage);
    }
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      iObjFileStr = *(argv + optind);
      if((iObjFileStr == NULL) || (*iObjFileStr == '\0'))
      {
        usage = 1;
      }
    }
  }
  ok = !usage;
  if(ok && (seedFileStr != NULL))
  {
    if((*seedFileStr == '\0') || 
       ((fP = (strcmp(seedFileStr, "-")?
	      fopen(seedFileStr, "r"): stdin)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_FILE_OPEN;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to open seed file %s (%s).\n",
	  *argv, seedFileStr, errMsg);
    }
    else
    {
      int		err = 0;

      seeds = WlzCCTParseSeedFile(fP, seeds, &dim, &nSeeds, &err);
      if(fP && strcmp(seedFileStr, "-"))
      {
        (void )fclose(fP);
      }
      if(err > 0)
      {
        ok = 0;
	errNum = WLZ_ERR_READ_INCOMPLETE;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	    "%s: Failed to read seeds from file %s,\n"
	    "possible error at or near line %d (%s).\n",
	    *argv, seedFileStr, err, errMsg);
      }
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if(((fP = (strcmp(iObjFileStr, "-")?
              fopen(iObjFileStr, "r"): stdin)) == NULL) ||
       ((iObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to read object from file %s (%s).\n",
	  *argv, iObjFileStr, errMsg);
    }
    if(fP && strcmp(iObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    oObj = WlzConComThreshold(iObj, nSeeds,
                              (dim == 2)? WLZ_VERTEX_I2: WLZ_VERTEX_I3,
			      seeds, tHiLo, xtra, rad, rHiLo, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to compute connected component threshold (%s).\n",
	  *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(oObjFileStr, "-")? fopen(oObjFileStr, "w"):
                                         stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, oObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to write output object (%s).\n",
	  *argv, errMsg);
    }
    if(fP && strcmp(oObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  AlcFree(seeds.v);
  (void )WlzFreeObj(iObj);
  (void )WlzFreeObj(oObj);
  if(usage != 0)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-d 2|3] [-h] [-r <radius>] [-R l|m|h]\n"
    "       [-s<seed point file>] [-S<seed point list>]\n"
    "       [-H|L|E] [-x <extra>] [-o<output file>] [<input object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -d  Seed point dimension, 2 or 3 for 2 or 3D seeds.\n"
    "  -E  Threshold equal, keep pixels/voxels with the threshold value\n"
    "      at the seed points.\n"
    "  -H  Threshold high, keep pixels/voxels at or above the threshold\n"
    "      value at the seed points (default).\n"
    "  -L  Threshold low, keep pixels/voxels at or below the threshold\n"
    "      value at the seed points (watch out that this is different\n"
    "      to WlzThreshold).\n"
    "  -h  Show this usage message.\n"
    "  -o  Output object file name (default -, stdout).\n"
    "  -r  Radius of sampling region around seed points (default = 0,\n"
    "      ie single pixel/voxel).\n"
    "  -R  Use the lowest (l), mean (m) or highest (h) value in the\n"
    "      sampling region as the threhold value (default mean).\n"
    "  -s  File containing seed point coordinates, one per line.\n"
    "  -S  Seed point coordinate list (comma separated).\n"
    "  -x  Integer percentage p to add/subtract to/from the grey value at\n"
    "      each seed,  with v(s) * ((100 +/- p) / 100). The sign is\n"
    "      choosen to keep the seed within the thresholded region also\n"
    "      p is ignored for threshold equal (default p = 0).\n"
    "Computes a connected component threshold domain for the given 2 or\n"
    "3D image object, in which the threshold domain is the union of all\n"
    "connected regions, thresholded using a value at a seed point, which\n"
    "include that seed point.\n"
    "A command line seed point list must be give the seed points as comma\n"
    "separated values, while seed points in a file may comma or space\n"
    "separated but must each be on a seperate line. The dimension of the\n"
    "seed points is only needed for those specified using a command line\n"
    "list. Valid syntax for the seed points in a seed point list is:\n"
    "  -d 2 -S <x1>,<y1>,<x2>,<y2>...,<xn>,<yn>\n"
    "or\n"
    "  -d 3 -S<x1>,<y1>,<z1><x2>,<y2>,<z2>...,<xn>,<yn>,<zn>\n"
    "Objects are read from stdin and written to stdout unless filenames\n"
    "are given.\n",
     *argv,
    " -d 2 -S 3,6,102,203 -r 2 -R l -s seedsfile -o out.wlz -L in.wlz\n"
    "The 2D seed points (3,6) and (102,203) are appended to the seed points\n"
    "read from the file \"seedsfile\". These are then used to threshold\n"
    "the 2D image read from the file \"in.wlz\" using the lowest value\n"
    "in a radius 2 circle around each seed point for the theshold value\n"
    "and forming the connected components using image values below these\n"
    "values.\n"
    "The output domain is then written to the file \"out.wlz\".\n");
  }
  return(!ok);
}

/*!
* \return	Array of vertices, either WlzIVertex2 or WlzIVertex3
* 		depending on the dimension.
* \ingroup	BinWlz
* \brief	Parses a comma separated vertex string (vertices and their
* 		coordinates are all comma separated) for the vertices.
* \param	seedListStr		Given seed list string.
* \param	dim			Dimension of the vertices.
* \param	nSeeds			Return for the number of seed vertices.
* \param	dstErr			Return for parse error, will be
* 					non-zero on error.
*/
static WlzVertexP		WlzCCTParseSeedList(
				  char *seedListStr,
				  int dim,
				  int *nSeeds,
				  int *dstErr)
{
  int		maxSeed,
  		idx = 0,
  		err = 0,
		seedCnt = 0;
  int		*buf = NULL;
  WlzVertexP 	seeds;

  seeds.v = NULL;
  if((dim != 2) && (dim != 3))
  {
    err = 1;
  }
  else if((seedListStr == NULL) ||
          ((maxSeed = strlen(seedListStr) / 2) < 2))
  {
  }
  else if((buf = (int *)AlcMalloc(sizeof(int) * maxSeed * dim)) == NULL)
  {
    err = 1;
  }
  else
  {
    char 	*pStr[3];

    pStr[0] = strtok(seedListStr, ",");
    while((err == 0) && (pStr[0] != NULL) && (seedCnt < maxSeed)) 
    {
      if((pStr[1] = strtok(NULL, ",")) == NULL)
      {
        err = 1;
      }
      else if((dim == 3) && ((pStr[2] = strtok(NULL, ",")) == NULL))
      {
        err = 1;
      }
      else if((sscanf(pStr[0], "%d", buf + idx) != 1) ||
	      (sscanf(pStr[1], "%d", buf + idx + 1) != 1) ||
	      ((dim == 3) && (sscanf(pStr[2], "%d", buf + idx + 2) != 1)))
      {
        err = 1;
      }
      else
      {
	++seedCnt;
        idx += dim;
	pStr[0] = strtok(NULL, ",");
      }
    }
  }
  if(!err && (seedCnt > 0))
  {
    size_t	seedSz;

    seedSz = (dim == 3)? sizeof(WlzIVertex3): sizeof(WlzIVertex2);
    if((seeds.v = AlcMalloc(seedCnt *  seedSz)) == NULL)
    {
      err = 1;
    }
    else
    {
      int	bIdx = 0;

      for(idx = 0; idx < seedCnt; ++idx)
      {

	if(dim == 3)
	{
	  seeds.i3[idx].vtX = buf[bIdx];
	  seeds.i3[idx].vtY = buf[bIdx + 1];
	  seeds.i3[idx].vtZ = buf[bIdx + 2];
	  bIdx += 3;
	}
	else
	{
	  seeds.i2[idx].vtX = buf[bIdx];
	  seeds.i2[idx].vtY = buf[bIdx + 1];
	  bIdx += 2;
	}
      }
    }
  }
  AlcFree(buf);
  if(err)
  {
    AlcFree(seeds.v);
    seeds.v = NULL;
    *nSeeds = 0;
    *dstErr = 1;
  }
  else
  {
    *nSeeds = seedCnt;
    *dstErr = 0;
  }
  return(seeds);
}

/*!
* \return	Array of vertices, either WlzIVertex2 or WlzIVertex3
* 		depending on the dimension.
* \ingroup	BinWlz
* \brief	Parses a comma separated vertex string from each line of
* 		the given (opened) file. The vertices must have comma
* 		separated coordinates and there must be one vertex per line.
* 		Lines starting with a has character are ignored.
* 		The parsed vertices are added to any existing vertices.
* \param	fP			Given seed file.
* \param	seeds			Existing seeds.
* \param	dim			Dimension of the vertices may be
* 					zero if unknown.
* \param	nSeeds			Number of seeds and the return for 
* 					the total number of seeds.
* \param	dstErr			Return for parse error, will be
* 					non-zero on error.
*/
static WlzVertexP		WlzCCTParseSeedFile(
				  FILE *fP,
				  WlzVertexP seeds,
				  int *dim,
				  int *nSeeds,
				  int *dstErr)
{
  int		maxSeed,
  		err = 0,
  		lnCnt = 0,
		first = 1;
  WlzVertex	seed;
  char		*ln;
  char		lnBuf[WLZ_CCT_READLN_LEN];
  const int	seedInc = 1024;

  maxSeed = *nSeeds;
  while((err == 0) && 
        (fgets(lnBuf, WLZ_CCT_READLN_LEN, fP) != NULL))
  {
    ++lnCnt;
    lnBuf[WLZ_CCT_READLN_LEN - 1] = '\0';
    ln = WlzStringWhiteSpSkip(lnBuf);
    if(ln && *ln && (*ln != '#'))
    {
      switch(*dim)
      {
	case 0:
	  if(first)
	  {
	    int        tDim;
	    WlzIVertex3 tSeed;
	    
	    tDim = sscanf(ln, "%d,%d,%d",
	                  &(tSeed.vtX), &(tSeed.vtY), &(tSeed.vtZ));
	    switch(tDim)
	    {
	      case 2:
		*dim = tDim;
	        seed.i2.vtX = tSeed.vtX;
		seed.i2.vtY = tSeed.vtY;
		break;
	      case 3:
		*dim = tDim;
	        seed.i3 = tSeed;
		break;
	      default:
	        err = lnCnt;
	    }
	    first = 0;
	  }
	  else
	  {
	    err = lnCnt;
	  }
        case 2:
	  if(sscanf(ln, "%d,%d",
	            &(seed.i2.vtX), &(seed.i2.vtY)) != 2)
	  {
	    err = lnCnt;
	  }
	  break;
        case 3:
	  if(sscanf(ln, "%d,%d,%d",
	            &(seed.i3.vtX), &(seed.i3.vtY), &(seed.i3.vtZ)) != 3)
	  {
	    err = lnCnt;
	  }
	  break;
      }
      if(err == 0)
      {
	/* Is there room for more seeds? */
	if(*nSeeds >= maxSeed)
	{
	  size_t  sz;

	  sz = (*dim == 3)? sizeof(WlzIVertex3): sizeof(WlzIVertex2);
	  maxSeed += seedInc;
	  if((seeds.v = AlcRealloc(seeds.v, maxSeed * sz)) == NULL)
	  {
	    err = lnCnt;
	  }
	}
      }
      /* Add seed. */
      if(err == 0)
      {
	if(*dim == 3)
	{
	  seeds.i3[*nSeeds] = seed.i3;
	}
	else
	{
	  seeds.i2[*nSeeds] = seed.i2;
	}
	++*nSeeds;
      }
    }
  }
  *dstErr = err;
  return(seeds);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
