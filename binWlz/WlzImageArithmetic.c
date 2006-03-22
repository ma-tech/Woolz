#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzImageArithmetic_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzImageArithmetic.c
* \author       Bill Hill
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
* \brief	Binary image arithmetic on a pair of domain objects.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzimagearithmetic "WlzImageArithmetic"
*/

/*!
\ingroup BinWlz
\defgroup wlzimagearithmetic WlzImageArithmetic
\par Name
WlzImageArithmetic - binary image arithmetic on a pair of domain objects.
\par Synopsis
\verbatim
WlzImageArithmetic [-o<out file>] [-a] [-d] [-l] [-g] [-m] [-s] [-n] [-N] [-h]
                   [-O#] [<in object 0>] [<in object 1>]
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
  <tr> 
    <td><b>-a</b></td>
    <td>Add the object's grey values.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Divide the grey values of the 1st object by those of the 2nd.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Vector magnitude of horizontal and vertical component objects.</td>
  </tr>
  <tr>
    <td><b>-l</b></td>
    <td>Compute the modulus of the grey values of the 1st object wrt
        those of the 2nd.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Multiply the object's grey values.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Subtract the grey values of the 2nd object from those of the 1st.</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>Normalises the output object to the range of the imput objects.</td>
  </tr>
  <tr>
    <td><b>-N</b></td>
    <td>Normalises the output objects to the range [0-255].</td>
  </tr>
  <tr>
    <td><b>-O</b></td>
    <td>Overwrite option (only useful for debugging).</td>
  </tr>
</table>
\par Description
Computes an arithmetic binary (two objects) operation on two domain
objects. The default operator is add.
The input objects are read from stdin and values are written to stdout
unless the filenames are given.
\par Examples
\verbatim
cat obj1.wlz | WlzImageArithmetic -o obj3.wlz -a - obj2.wlz
\endverbatim
A new object obj3.wlz is formed by adding the grey values of obj1 and
obj2.wlz.
\par File
\ref WlzImageArithmetic.c "WlzImageArithmetic.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzImageArithmetic "WlzImageArithmetic(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

typedef enum
{
  WLZ_IMARNORM_NONE,
  WLZ_IMARNORM_INPUT,
  WLZ_IMARNORM_256
} WlzImArNorm;

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		idx,
		option,
		overwrite = 1,
		ok = 1,
		usage = 0;
  WlzImArNorm	norm = WLZ_IMARNORM_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*outObj = NULL;
  WlzObject	*inObj[2];
  WlzBinaryOperatorType operator = WLZ_BO_ADD;
  char 		*outObjFileStr;
  char  	*inObjFileStr[2];
  WlzPixelV	gMin[3],
  		gMax[3];
  const char	*errMsg;
  static char	optList[] = "O:o:ab:dglmsnNh",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObj[0] = NULL;
  inObj[1] = NULL;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr[0] = inObjFileStrDef;
  inObjFileStr[1] = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
    case 'O':
      if((sscanf(optarg, "%d", &overwrite) != 1) ||
	 (overwrite < 0) || (overwrite > 2))
      {
	usage = 1;
	ok = 0;
      }
      break;

    case 'o':
      outObjFileStr = optarg;
      break;

    case 'a':
      operator = WLZ_BO_ADD;
      break;

    case 'b':
      operator = atoi(optarg);
      switch( operator ){
      case WLZ_BO_ADD:
      case WLZ_BO_SUBTRACT:
      case WLZ_BO_MULTIPLY:
      case WLZ_BO_DIVIDE:
      case WLZ_BO_MODULUS:
      case WLZ_BO_EQ:
      case WLZ_BO_NE:
      case WLZ_BO_GT:
      case WLZ_BO_GE:
      case WLZ_BO_LT:
      case WLZ_BO_LE:
      case WLZ_BO_AND:
      case WLZ_BO_OR:
      case WLZ_BO_XOR:
      case WLZ_BO_MAX:
      case WLZ_BO_MIN:
      case WLZ_BO_MAGNITUDE:
	break;

      default:
	usage = 1;
	ok = 0;
	break;
      }
      break;

    case 'd':
      operator = WLZ_BO_DIVIDE;
      break;

    case 'g':
      operator = WLZ_BO_MAGNITUDE;
      break;

    case 'l':
      operator = WLZ_BO_MODULUS;
      break;

    case 'm':
      operator = WLZ_BO_MULTIPLY;
      break;

    case 's':
      operator = WLZ_BO_SUBTRACT;
      break;

    case 'n':
      norm = WLZ_IMARNORM_INPUT;
      break;

    case 'N':
      norm = WLZ_IMARNORM_256;
      break;

    case 'h':
    default:
      usage = 1;
      ok = 0;
      break;
    }
  }
  if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
     (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    idx = 0;
    while((idx < 2) && (optind < argc))
    {
      inObjFileStr[idx] = *(argv + optind);
      ++optind;
      ++idx;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
    ok = 0;
  }
  if(ok)
  {
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx]= WlzAssignObject(WlzReadObj(fP,
	  					   &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to read object %d from file %s (%s)\n",
		       *argv, idx, inObjFileStr[idx], errMsg);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
      }
      ++idx;
    }
  }
  if(ok)
  {
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      if((inObj[idx]->type != WLZ_2D_DOMAINOBJ) &&
         (inObj[idx]->type != WLZ_3D_DOMAINOBJ))
      {
        ok = 0;
	(void )fprintf(stderr, "%s: input object %d is not a domain object\n",
		       *argv, idx);
      }
      ++idx;
    }
  }
  if(ok && (norm == WLZ_IMARNORM_INPUT))
  {
    if(((errNum = WlzGreyRange(inObj[0], gMin, gMax)) != WLZ_ERR_NONE) ||
       ((errNum = WlzGreyRange(inObj[1], gMin + 1, gMax + 1)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to find input objects grey range (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if((outObj = WlzImageArithmetic(inObj[0], inObj[1], operator, 1,
				    &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute new object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok && (norm != WLZ_IMARNORM_NONE))
  {
    if((errNum = WlzGreyRange(outObj, gMin+ 2, gMax + 2)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to find output object's grey range (%s)\n",
		     *argv, errMsg);
    }
    if(ok)
    {
      if(norm == WLZ_IMARNORM_256)
      {
        gMin[0].type = WLZ_GREY_UBYTE;
	gMin[0].v.ubv = 0;
        gMax[0].type = WLZ_GREY_UBYTE;
	gMax[0].v.ubv = 255;
	(void )WlzValueConvertPixel(gMin, gMin[0], gMin[2].type);
	(void )WlzValueConvertPixel(gMax, gMax[0], gMax[2].type);
      }
      else
      {
	(void )WlzValueConvertPixel(gMin, gMin[0], gMin[2].type);
	(void )WlzValueConvertPixel(gMax, gMax[0], gMax[2].type);
	(void )WlzValueConvertPixel(gMin + 1, gMin[1], gMin[2].type);
	(void )WlzValueConvertPixel(gMax + 1, gMax[1], gMax[2].type);
	switch(gMin[2].type)
	{
	  case WLZ_GREY_INT:
	    if(gMin[1].v.inv < gMin[0].v.inv)
	    {
	      gMin[0].v.inv = gMin[1].v.inv;
	    }
	    if(gMax[1].v.inv > gMax[0].v.inv)
	    {
	      gMax[0].v.inv = gMax[1].v.inv;
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    if(gMin[1].v.shv < gMin[0].v.shv)
	    {
	      gMin[0].v.shv = gMin[1].v.shv;
	    }
	    if(gMax[1].v.shv > gMax[0].v.shv)
	    {
	      gMax[0].v.shv = gMax[1].v.shv;
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    if(gMin[1].v.ubv < gMin[0].v.ubv)
	    {
	      gMin[0].v.ubv = gMin[1].v.ubv;
	    }
	    if(gMax[1].v.ubv > gMax[0].v.ubv)
	    {
	      gMax[0].v.ubv = gMax[1].v.ubv;
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    if(gMin[1].v.flv < gMin[0].v.flv)
	    {
	      gMin[0].v.flv = gMin[1].v.flv;
	    }
	    if(gMax[1].v.flv > gMax[0].v.flv)
	    {
	      gMax[0].v.flv = gMax[1].v.flv;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    if(gMin[1].v.dbv < gMin[0].v.dbv)
	    {
	      gMin[0].v.dbv = gMin[1].v.dbv;
	    }
	    if(gMax[1].v.dbv > gMax[0].v.dbv)
	    {
	      gMax[0].v.dbv = gMax[1].v.dbv;
	    }
	    break;
	}
      }
    }
    if((errNum = WlzGreySetRange(outObj, gMin[2], gMax[2],
    				 gMin[0], gMax[0])) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to set output object's grey range (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?  fopen(outObjFileStr, "w"):
	      				    stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  WlzFreeObj(inObj[0]);
  WlzFreeObj(inObj[1]);
  WlzFreeObj(outObj);
  if(usage)
  {

    fprintf(stderr,
	    "Usage: %s"
	    " [-O#] [-o<out file>] [-a] [-b <op>] [-d] [-l] [-m] [-s] [-h] "
	    "[<in object 0>] [<in object 1>]\n"
	    "Options:\n"
	    "  -O        Overwrite option (only useful for debugging)\n"
	    "  -o        Output file name.\n"
	    "  -a        Add the object's grey values.\n"
	    "  -b <op>   Apply binary operation to the grey-values.\n"
	    "            op =  %d - Add\n"   
	    "               =  %d - SUBTRACT\n"
	    "               =  %d - MULTIPLY\n"
	    "               =  %d - DIVIDE\n"
	    "               =  %d - MODULUS\n"
	    "               =  %d - EQ\n"
	    "               =  %d - NE\n"
	    "               =  %d - GT\n"
	    "               =  %d - GE\n"
	    "               =  %d - LT\n"
	    "               =  %d - LE\n"
	    "               =  %d - AND\n"
	    "               =  %d - OR\n"
	    "               =  %d - XOR\n"
	    "               =  %d - MAX\n"
	    "               =  %d - MIN\n"
	    "               =  %d - MAGNITUDE\n"
	    "  -d        Divide the grey values of the 1st object by those of the 2nd.\n"
	    "  -g        Vector magnitude of horizontal and vertical component objects\n"
	    "  -l        Compute the modulus of the grey values of the 1st object wrt\n"
	    "            those of the 2nd.\n"
	    "  -m        Multiply the object's grey values.\n"
	    "  -s        Subtract the grey values of the 2nd object from those of the 1st.\n"
	    "  -n        Normalises the output object to the range of the input objects.\n"
	    "  -N        Normalises the output objects to the range [0-255].\n"
	    "  -h        Help, prints this usage message.\n"
	    "Computes an arithmetic binary (two objects) operation on two domain\n"
	    "objects. The default operator is add.\n"
	    "The input objects are read from stdin and values are written to stdout\n"
	    "unless the filenames are given.\n"
	    "Example:\n%s%s%s",
	    *argv,
	    WLZ_BO_ADD,
	    WLZ_BO_SUBTRACT,
	    WLZ_BO_MULTIPLY,
	    WLZ_BO_DIVIDE,
	    WLZ_BO_MODULUS,
	    WLZ_BO_EQ,
	    WLZ_BO_NE,
	    WLZ_BO_GT,
	    WLZ_BO_GE,
	    WLZ_BO_LT,
	    WLZ_BO_LE,
	    WLZ_BO_AND,
	    WLZ_BO_OR,
	    WLZ_BO_XOR,
	    WLZ_BO_MAX,
	    WLZ_BO_MIN,
	    WLZ_BO_MAGNITUDE,
	    "cat obj1.wlz | ",
	    *argv,
	    " -o obj3.wlz -a - obj2.wlz\n"
	    "A new object 'obj3.wlz is formed by adding the grey values of obj1 and\n"
	    "obj2.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
