#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzLUTGreyTransformObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzLUTGreyTransformObj.c
* \author       Bill Hill
* \date         November 2011
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
* \brief	Applies an look up table grey transform.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzlutgreytransformobj "WlzLUTGreyTransformObj"
*/
 
/*!
\ingroup      BinWlz
\defgroup     wlzlutgreytransformobj WlzLUTGreyTransformObj
\par Name
WlzLUTGreyTransformObj -
Reads, writes, computes and applies grey value look up tables.
\par Synopsis
\verbatim
WlzLUTGreyTransformObj [-o<output object file>] [-h]
           [-d] [-i] [-A] [-N] [-T] [-f E|G|I|L|S] [-g i|r] [-G i|s|u|r]
	   [-l#] [-u#] [-L#] [-U#] [-p#] [-q#] [-t <input LUT>]
           [<input object>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints a usage message.</td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Dither grey values.</td>
  </tr>
  <tr>
    <td><b>-i</b></td>
    <td>Apply LUT grey value transform in place.</td>
  </tr>
  <tr>
    <td><b>-A</b></td>
    <td>Output the LUT as ASCII text to the output file.</td>
  </tr>
  <tr>
    <td><b>-N</b></td>
    <td>No LUT transformation (input object is not required).</td>
  </tr>
  <tr>
    <td><b>-T</b></td>
    <td>Output the LUT as a Woolz object to the output file.</td>
  </tr>
  <tr>
    <td><b>-f</b></td>
    <td>Grey transform function with which to set the look up table
        specified by a single letter (case insensitive and rest of
	word is ignored):</td>
    <br> <b>g</b>amma,
    <br> <b>i</b>dentity,
    <br> <b>l</b>inear,
    <br> <b>s</b>igmoid.
  </tr>
  <tr>
    <td><b>-g</b></td>
    <td>LUT grey type specified by a single letter:</td>
    <br> <b>i</b> (32 bit) integer,
    <br> <b>r</b> (32 bit) red,green,blue,alpha colour,
  </tr>
  <tr>
    <td><b>-G</b></td>
    <td>Output grey type specified by a single letter:</td>
    <br> <b>i</b> (32 bit) integer,
    <br> <b>s</b> (16 bit) short integer,
    <br> <b>u</b> (8  bit) unsigned byte,
    <br> <b>r</b> (32 bit) red,green,blue,alpha colour.
  </tr>
  <tr>
    <td><b>-l</b></td>
    <td>Lower input grey value.</td>
  </tr>
  <tr>
    <td><b>-u</b></td>
    <td>Upper input grey value.</td>
  </tr>
  <tr>
    <td><b>-L</b></td>
    <td>Lower output grey value.</td>
  </tr>
  <tr>
    <td><b>-H</b></td>
    <td>Upper output grey value.</td>
  </tr>
  <tr>
    <td><b>-p</b></td>
    <td>First grey transform parameter. This is unused for identity
        and linear transforms, but is \f$\gamma\f$ or sigmoid \f$\mu\f$.</td>
  </tr>
  <tr>
    <td><b>-q</b></td>
    <td>Second grey transform parameter. This is unused for identity,
        linear and gamma transforms, but is \f$\sigma\f$ for sigmoid
	transforms.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Read LUT from a Woolz object file.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>output object file name.</td>
  </tr>
</table>

\par Description
Reads or computes a look up table object and then either applies
this to a grey value object or writes out the lookup table object.
If a look up table object is specified on the command line then none
of the command line grey transform primatives are used.
By default files are read from the standard input and written to
the standard output.

\par Examples
  # Remaps the grey values of the given input object so that they
  # are unsigned bytes, mapped from the input range 0-4095 to the
  # output range 0-255, using a Gamma function with a Gamma value
  # of 0.6.
  WlzLUTGreyTransformObj -l 0 -u 4095 -L 0 -U 255 -f G -g i -G u
                         -p 0.6 -o out.wlz in.wlz
\par File
\ref WlzLUTGreyTransformObj.c "WlzLUTGreyTransformObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzfacts "WlzFacts(1)"
\ref WlzLUTGreyTransformNew "WlzLUTGreyTransformNew(3)"
\ref WlzLUTTransformObj "WlzLUTTransformObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

static int			WlzLGTOParsePixel(
				  char *str,
				  WlzPixelV *pix);

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		iIL,
  		iIU,
		option,
		dither = 0,
		inplace = 0,
		outputLUTObj = 0,
		outputLUTAscii = 0,
		noLUTTransform = 0,
                setLUTTransform = 0,
		ok = 1,
		usage = 0;
  double	p = 1.0,
  		q = 1.0;
  WlzGreyType	lGType = WLZ_GREY_INT,
  		oGType = WLZ_GREY_UBYTE;
  WlzPixelV	gOL,
		gOU;
  WlzGreyTransformType gTrType = WLZ_GREYTRANSFORMTYPE_LINEAR;
  WlzObject	*inObj = NULL,
  		*outObj = NULL,
		*lutObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  const char	*errMsg;
  char 		*outFileStr,
  		*inObjFileStr,
		*inLUTFileStr = NULL;
  static char	optList[] = "hdiANTf:g:G:o:l:u:L:U:p:q:t:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
  iIL = gOL.v.inv = 0;
  iIU = gOU.v.inv = 255;
  gOL.type = gOU.type = WLZ_GREY_INT;
  /* Parse the command line. */
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'd':
        dither = 1;
	break;
      case 'i':
        inplace = 1;
	break;
      case 'f':
	setLUTTransform = 1;
        switch(*optarg)
	{
	  case 'g': /* FALLTHROUGH */
	  case 'G':
	    gTrType = WLZ_GREYTRANSFORMTYPE_GAMMA;
	    break;
	  case 'i': /* FALLTHROUGH */
	  case 'I':
	    gTrType = WLZ_GREYTRANSFORMTYPE_IDENTITY;
	    break;
	  case 'l': /* FALLTHROUGH */
	  case 'L':
	    gTrType = WLZ_GREYTRANSFORMTYPE_LINEAR;
	    break;
	  case 's': /* FALLTHROUGH */
	  case 'S':
	    gTrType = WLZ_GREYTRANSFORMTYPE_SIGMOID;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'g':
        switch(*optarg)
	{
	  case 'i': /* FALLTHROUGH */
	  case 'I':
	    lGType = WLZ_GREY_INT;
	    break;
	  case 'r': /* FALLTHROUGH */
	  case 'R':
	    lGType = WLZ_GREY_RGBA;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'G':
        switch(*optarg)
	{
	  case 'i': /* FALLTHROUGH */
	  case 'I':
	    oGType = WLZ_GREY_INT;
	    break;
	  case 's': /* FALLTHROUGH */
	  case 'S':
	    oGType = WLZ_GREY_SHORT;
	    break;
	  case 'u': /* FALLTHROUGH */
	  case 'U':
	    oGType = WLZ_GREY_UBYTE;
	    break;
	  case 'r': /* FALLTHROUGH */
	  case 'R':
	    oGType = WLZ_GREY_RGBA;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'l':
	setLUTTransform = 1;
	if(sscanf(optarg, "%d", &iIL) != 1)
	{
	  usage = 1;
	}
	break;
      case 'u':
	setLUTTransform = 1;
	if(sscanf(optarg, "%d", &iIU) != 1)
	{
	  usage = 1;
	}
	break;
      case 'L':
	setLUTTransform = 1;
	usage = WlzLGTOParsePixel(optarg, &gOL);
	break;
      case 'U':
	setLUTTransform = 1;
	usage = WlzLGTOParsePixel(optarg, &gOU);
	break;
      case 'p':
	setLUTTransform = 1;
	if(sscanf(optarg, "%lg", &p) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'q':
	setLUTTransform = 1;
	if(sscanf(optarg, "%lg", &q) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'N':
      	noLUTTransform = 1;
	break;
      case 'A':
      	outputLUTAscii = 1;
      	outputLUTObj = 0;
	break;
      case 'T':
      	outputLUTAscii = 0;
      	outputLUTObj = 1;
	break;
      case 't':
        inLUTFileStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	break;
    }
  }
  if((usage == 0) && (noLUTTransform == 0))
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
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
	inObjFileStr = *(argv + optind);
      }
    }
  }
  ok = (usage == 0);
  /* Read if input LUT file given. */
  if(ok)
  {
    if(inLUTFileStr)
    {
      errNum = WLZ_ERR_READ_EOF;
      if(((fP = (strcmp(inLUTFileStr, "-")?
		fopen(inLUTFileStr, "r"): stdin)) == NULL) ||
	 ((lutObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to read LUT object from file %s (%s)\n",
		       *argv, inLUTFileStr, errMsg);
      }
      if(fP && strcmp(inLUTFileStr, "-"))
      {
	(void )fclose(fP);
        fP = NULL;
      }
      if(lutObj)
      {
        if(lutObj->type != WLZ_LUT)
	{
          errNum = WLZ_ERR_OBJECT_TYPE;
	}
	else if(lutObj->domain.core == NULL)
	{
          errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(lutObj->values.core == NULL)
	{
          errNum = WLZ_ERR_VALUES_NULL;
	}
	else if(lutObj->domain.lut->type != WLZ_LUT)
	{
          errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if(lutObj->values.lut->type != WLZ_LUT)
	{
          errNum = WLZ_ERR_VALUES_TYPE;
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr, "%s: Inappropriate LUT object (%s)\n",
	                 *argv, errMsg);
	}
      }
    }
  }
  /* Create new LUT object if one wasn't read from file otherwise set
   * the LUT object values within given range. */
  if(ok)
  {
    if(gOL.type != lGType)
    {
      errNum = WlzValueConvertPixel(&gOL, gOL, lGType);
    }
    if((errNum == WLZ_ERR_NONE) && (gOU.type != lGType))
    {
      errNum = WlzValueConvertPixel(&gOU, gOU, lGType);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if( lutObj )
      {
	if( setLUTTransform )
	{
	  errNum = WlzLUTGreyTransformSet(lutObj, gTrType,
					  lGType, iIL, iIU, gOL.v, gOU.v,
					  p, q);
	}
      }
      else 
      {
	lutObj = WlzLUTGreyTransformNew(gTrType, lGType,
					iIL, iIU, gOL.v, gOU.v,
					p, q, &errNum);
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;   
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to create/set LUT object (%s)\n",
		     *argv, errMsg);
    }
  }
  /* Optionaly output the LUT object either as ASCII or a Woolz object. */
  if(ok && outputLUTAscii)
  {
    if((fP = (strcmp(outFileStr, "-")?
	     fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output ASCII LUT to file %s\n",
		     *argv, outFileStr);
    }
    else
    {
      int	idx,
	      bin1,
	      nBin;
      WlzGreyP gP;

      gP = lutObj->values.lut->val;
      bin1 = lutObj->domain.lut->bin1;
      nBin = lutObj->domain.lut->lastbin - bin1 + 1;
      switch(lutObj->values.lut->vType)
      {
	case WLZ_GREY_INT:
	  (void )fprintf(fP, "WLZ_GREY_INT\n");
	  for(idx = 0; idx < nBin; ++idx)
	  {
	    if(fprintf(fP, "% 8d % 8d\n", bin1 + idx, gP.inp[idx]) <= 0) 
	    {
	      ok = 0;
	      break;
	    }
	  }
	  break;
	case WLZ_GREY_SHORT:
	  (void )fprintf(fP, "WLZ_GREY_SHORT\n");
	  for(idx = 0; idx < nBin; ++idx)
	  {
	    if(fprintf(fP, "% 8d % 8d\n", bin1 + idx, gP.shp[idx]) <= 0) 
	    {
	      ok = 0;
	      break;
	    }
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  (void )fprintf(fP, "WLZ_GREY_UBYTE\n");
	  for(idx = 0; idx < nBin; ++idx)
	  {
	    if(fprintf(fP, "% 8d % 3d\n", bin1 + idx, gP.ubp[idx]) <= 0) 
	    {
	      ok = 0;
	      break;
	    }
	  }
	  break;
	case WLZ_GREY_RGBA:
	  (void )fprintf(fP, "WLZ_GREY_RGBA\n");
	  for(idx = 0; idx < nBin; ++idx)
	  {
	    if(fprintf(fP, "% 8d % 3d % 3d % 3d % 3d\n",
		       bin1 + idx,
		       WLZ_RGBA_RED_GET(gP.inp[idx]),
		       WLZ_RGBA_GREEN_GET(gP.inp[idx]),
		       WLZ_RGBA_BLUE_GET(gP.inp[idx]),
		       WLZ_RGBA_ALPHA_GET(gP.inp[idx])) <= 0) 
	    {
	      ok = 0;
	      break;
	    }
	  }
	  break;
	default:
	  ok = 0;
	  break;
      }
      if(ok == 0)
      {
	(void )fprintf(stderr,
		       "%s: Failed to write ASCII LUT to file %s\n",
		       *argv, outFileStr);
      }
    }
  }
  if(ok && outputLUTObj)
  {
    if(((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL) ||
	(WlzWriteObj(fP, lutObj) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write LUT object to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
      fP = NULL;
    }
  }
  if(ok && (noLUTTransform == 0) &&
     (outputLUTAscii == 0) && (outputLUTObj == 0))
  {
    errNum = WLZ_ERR_READ_EOF;
    if(((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s)\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP);
      fP = NULL;
    }
    if(ok)
    {
      outObj = WlzLUTTransformObj(inObj, lutObj, oGType, inplace, dither,
				  &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
        (void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to use LUT to transform grey values (%s)\n",
		       *argv, errMsg);

      }
    }
    if(ok)
    {
      if(((fP = (strcmp(outFileStr, "-")?
		fopen(outFileStr, "w"): stdout)) == NULL) ||
	  (WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write output object to file %s (%s).\n",
		       *argv, outFileStr, errMsg);
      }
      if(fP && strcmp(outFileStr, "-"))
      {
	fclose(fP);
	fP = NULL;
      }
    }
  }
  WlzFreeObj(inObj);
  WlzFreeObj(lutObj);
  WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-o<output object>] [-h]\n"
    "         [-d] [-i] [-A] [-N] [-T] [-f E|G|I|L|S] [-g i|r] [-G i|s|u|r]\n"
    "         [-l#] [-u#] [-L#] [-U#] [-p#] [-q#] [-t <input LUT>]\n"
    "         [<input object>]\n" 
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -d  Dither grey values.\n"
    "  -i  Apply LUT grey value transform in place.\n"
    "  -A  Output the LUT as ASCII text to the output file.\n"
    "  -N  No LUT transformation (input object is not required).\n"
    "  -T  Output the LUT as a Woolz object to the output file.\n"
    "  -f  Grey transform function with which to set the look up table\n"
    "      specified by a single letter (case insensitive and rest of the\n"
    "      word is ignored):\n"
    "        *g*amma,\n"
    "        *i*dentity,\n"
    "        *l*inear,\n"
    "        *s*igmoid.\n"
    "  -g  LUT grey type specified by a single letter:\n"
    "        i (32 bit) integer,\n"
    "        r (32 bit) red,green,blue,alpha colour.\n"
    "  -G  Output grey type specified by a single letter:\n"
    "        i (32 bit) integer,\n"
    "        s (16 bit) short integer,\n"
    "        u (8  bit) unsigned byte,\n"
    "        r (32 bit) red,green,blue,alpha colour.\n"
    "  -l  Lower input grey value.\n"
    "  -u  Upper input grey value.\n"
    "  -L  Lower output grey value.\n"
    "  -U  Upper output grey value.\n"
    "  -p  First grey transform parameter. This is unused for identity\n"
    "      and linear transforms, but is gamma, or sigmoid mu.\n"
    "  -q  Second grey transform parameter. This is unused for identity,\n"
    "      linear and gamma transforms, but is the sigmoid sigma.\n"
    "  -t  Read LUT from a Woolz object file.\n"
    "  -o  Output object file name.\n"
    "Reads or computes a look up table object and then either applies\n"
    "this to a grey value object or writes out the lookup table object.\n"
    "If a look up table object is specified on the command line then none\n"
    "of the command line grey transform primatives are used.\n"
    "By default files are read from the standard input and written to\n"
    "the standard output.\n",
    *argv,
    " -l 0 -u 4095 -L 0 -U 255 -f G -g i -G u\n"
    "                                -p 0.6 -o out.wlz in.wlz\n"
    "Remaps the grey values of the given input object so that they\n"
    "are unsigned bytes, mapped from the input range 0-4095 to the\n"
    "output range 0-255, using a Gamma function with a Gamma value\n"
    "of 0.6.\n");
  }
  return(!ok);
}

/*!
* \return	Non-zero on parse error.
* \ingroup
* \brief	Parses the given string for a grey value. This is set
* 		to WLZ_GREY_RGBA if the given string contains a comma otherwise
* 		the value is WLZ_GREY_INT.
* \param	str			Given string to parse.
* \param	pix			Pixel value to be set.
*/
static int	WlzLGTOParsePixel(char *str, WlzPixelV *pix)
{
  int		usage = 1;
  int		t[4];

  if(str && pix)
  {
    if(strchr(optarg, ','))
    {
      t[3] = 255;
      t[0] = t[1] = t[2] = 0;
      if((sscanf(optarg, "%d,%d,%d,%d",
	      &(t[0]), &(t[1]), &(t[2]), &(t[3])) >= 3) &&
	  (t[0] >= 0) && (t[1] >= 0) && (t[2] >= 0) && (t[3] >= 0) && 
	  (t[0] <= 255) && (t[1] <= 255) && (t[2] <= 255) && (t[3] <= 255))
      {
	usage = 0;
	pix->type = WLZ_GREY_RGBA;
	WLZ_RGBA_RGBA_SET(pix->v.rgbv, t[0], t[1], t[2], t[3]);
      }
    }
    else
    {
      if(sscanf(optarg, "%d", &(pix->v.inv)) == 1)
      {
	usage = 0;
        pix->type = WLZ_GREY_INT;
      }
    }
  }
  return(usage);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
