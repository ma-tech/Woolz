#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstBuildObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzTstBuildObj.c
* \author       Bill Hill
* \date         August 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Test program for WlzBuildObj3().
* \ingroup	BinWlzTst
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

static int			ParseDimLim(
				  int *dstMin,
				  int *dstMax,
				  char *str);

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;

int		main(int argc, char *argv[])
{
  int		tSz1,
  		option,
		nT = 255,
		debug = 0,
  		ok = 1,
  		usage = 0;
  char		*outFileS;
  WlzGreyP	bufP;
  WlzIVertex2	tSz;
  WlzIBox3	bBox;
  FILE		*fP = NULL;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgS = NULL;
  static char   optList[] = "dhn:o:s:x:y:z:";

  opterr = 0;
  tSz.vtX = tSz.vtY = 64;
  bBox.xMin = 0;
  bBox.yMin = 0;
  bBox.zMin = 0;
  bBox.xMax = 255;
  bBox.yMax = 255;
  bBox.zMax = 255;
  bufP.v = NULL;
  outFileS = "-";
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
	debug = 1;
	break;
      case 'o':
	outFileS = optarg;
	break;
      case 'n':
	if((sscanf(optarg, "%d", &nT) != 1) || (nT < 0))
	{
	  usage = 1;
	}
	break;
      case 's':
	if((sscanf(optarg, "%d,%d", &(tSz.vtX), &(tSz.vtY)) != 2) ||
	   (tSz.vtX < 0) || (tSz.vtY < 0))
	{
	  usage = 1;
	}
	break;
      case 'x':
        usage = ParseDimLim(&(bBox.xMin), &(bBox.xMax), optarg);
	break;
      case 'y':
        usage = ParseDimLim(&(bBox.yMin), &(bBox.yMax), optarg);
	break;
      case 'z':
        usage = ParseDimLim(&(bBox.zMin), &(bBox.zMax), optarg);
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  tSz1 = tSz.vtX * tSz.vtY;
  if(debug)
  {
    (void )fprintf(stderr, "nT   = %d\n",
		   nT);
    (void )fprintf(stderr, "tSz  = %d,%d\n",
		   tSz.vtX, tSz.vtY);
    (void )fprintf(stderr, "bBox = %d,%d,%d,%d,%d,%d\n",
		   bBox.xMin, bBox.xMax,
		   bBox.yMin, bBox.yMax,
		   bBox.zMin, bBox.zMax);
  }
  ok = !usage;
  if(ok)
  {
    if((bufP.v = AlcMalloc(sizeof(WlzUByte) * tSz1)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to allocate buffer.\n",
                     *argv);
    }
  }
  if(ok)
  {
    int		idN;

    AlgRandSeed(0L);
    for(idN = 0; idN < nT; ++idN)
    {
      WlzIVertex3 og;
      WlzObject *nObj = NULL;

      og.vtX = bBox.xMin +
               (int )floor(AlgRandUniform() *
	                   (bBox.xMax - bBox.xMin - tSz.vtX));
      og.vtY = bBox.yMin +
               (int )floor(AlgRandUniform() *
	                   (bBox.yMax - bBox.yMin - tSz.vtY));
      og.vtZ = bBox.zMin + 
               (int )floor(AlgRandUniform() *
	                   (bBox.zMax - bBox.zMin));
      (void )memset(bufP.v, idN % 255, tSz1);
      nObj = WlzBuildObj3(obj, og, tSz, WLZ_GREY_UBYTE, tSz1, bufP,
      			  &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsgS);
	(void )fprintf(stderr,
	               "%s: Failed to add tile #%d to object (%s)\n",
		       *argv, idN, errMsgS);
        break;
      }
      else
      {
        (void )WlzFreeObj(obj);
	obj = nObj;
      }
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outFileS, "-")?
        fopen(outFileS, "w"): stdout)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_WRITE_EOF;
      (void )WlzStringFromErrorNum(errNum, &errMsgS);
      (void )fprintf(stderr,
                     "%s: Failed to open file %s (%s).\n",
		     *argv, outFileS, errMsgS);
    }
  }
  if(ok)
  {
    if((errNum = WlzWriteObj(fP, obj)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgS);
      (void )fprintf(stderr,
                     "%s: Failed to write to file %s (%s).\n",
		     *argv, outFileS, errMsgS);
    }
  }
  if(fP && strcmp(outFileS, "-"))
  {
    (void )fclose(fP);
  }
  AlcFree(bufP.v);
  (void )WlzFreeObj(obj);
  if(usage)
  {
  static char   optList[] = "dhn:o:s:x:y:z:";
    (void )fprintf(stderr,
    "Usage: %s [-d] [-h] [-n#] [-o<file>] [-s#,#]\n"
    "\t\t[-x[#][,[#]]] [-y[#][,[#]]] [-z[#][,[#]]]\n"
    "Creates 3D domain object to test WlzBuildObj3().\n"
    "Options are:\n"
    "  -d  Print debug information to stderr (set to %d).\n"
    "  -h  Help, prints this usage message.\n"
    "  -n  NUmber of 2D rectangular buffers to add (set to %d).\n"
    "  -o  Output file (set to %s).\n"
    "  -s  Rectangular buffer size (side lengths) (set to %d,%d).\n"
    "  -x  Bounding column limits (set to %d,%d).\n"
    "  -y  Bounding line limits (set to %d,%d).\n"
    "  -z  Bounding plane limits (set to %d,%d).\n",
    argv[0], debug, nT, outFileS, tSz.vtX, tSz.vtY,
    bBox.xMin, bBox.xMax, bBox.yMin, bBox.yMax, bBox.zMin, bBox.zMax);
  }
  return(!ok);
}

/*!
* \return	Non zero if there's a parse error.
* \ingroup	WlzTst
* \brief	Parses the given string for [<min>][,[<max>]] setting
* 		the destination min and max values given.
* \param	dstMin			Destination pointer for minimum value.
* \param	dstMax			Destination pointer for maximum value.
* \param	str			Given string to parse.
*/
static int	ParseDimLim(int *dstMin, int *dstMax, char *str)
{
  char		*s0,
  		*s1;
  int		err = 0;

  str = WlzStringWhiteSpSkip(str);
  if(*str == ',')
  {
    s0 = NULL;
    s1 = strtok(str, ",");
  }
  else
  {
    s0 = strtok(str, ",");
    s1 = strtok(NULL, ",");
  }
  if(s0 && (sscanf(s0, "%d", dstMin) != 1))
  {
    err = 1;
  }
  if((err ==  0) && s1 && (sscanf(s1, "%d", dstMax) != 1))
  {
    err = 2;
  }
  return(err);
}
