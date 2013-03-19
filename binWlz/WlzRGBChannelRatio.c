#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRGBChannelRatio_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzRGBChannelRatio.c
* \author       Bill Hill
* \date         February 2008
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
* \brief	Computes log ratio of RGB channels in a RGBA object.
* \ingroup	BinWlz
*/

/*!
\ingroup      BinWlz
\defgroup     wlzrgbchannelratio WlzRGBChannelRatio
\par Name
WlzRGBChannelRatio - computes log ratio of RGB channels in a RGBA object.
\par Synopsis
\verbatim
WlzRGBChannelRatio [-h] [-o<output object file>]
                   [-d<channel>] [-n<channel>] [-m<channel>]
		   [-N] [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Denominator, must be a valid channel letter. </td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>Numerator, must be a valid channel letter. </td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Multiplier, must be a valid channel letter. </td>
  </tr>
  <tr>
    <td><b>-N</b></td>
    <td>Normalise the values and convert to unsigned byte. </td>
  </tr>
</table>
\par Description
Computes log ratio of RGB channels in a RGBA object for each
pixel using ratio \f$r\f$ with
\f[
r = \log(1 + \frac{n}{1 + d}).
\f]
or if a multiplier channel is specified
\f[
r = m \log(1 + \frac{n}{1 + d}).
\f]

This results in an object with either float values or values normalised
to the range 0-255 and converted to ubyte.
Colour channels are specified by the ** letter of the colour channel
strings *r*ed, *g*reen, *b*lue, *y*ellow, *m*agenta, "*c*yan, *h*ue,
*s*aturation, *v*alue and gr*e*y, ie r, g, b, y, m, c, h, s, v and e.
No multiplier is used unless given on the command line.
By default the input object is read from the standard input and the
output object is written to standard the output.
\par Examples
\verbatim
# The following commad computes a new object from a RGBA object in.wlz.
# In the resulting object the pixel values are the ratio of the input
#  red and green components evaluated at each pixel. The ratio object
# is written to out.obj.
WlzRGBChannelRatio -n r -d g -o out.wlz in.wlz

\endverbatim

\par File
\ref WlzRGBChannelRatio.c "WlzRGBChannelRatio.c"
\par See Also
\ref WlzRGBChanRatio "WlzRGBChanRatio(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

static int			WlzRGBRatioChan(
				  WlzRGBAColorChannel *chan,
				  char *str);

int		main(int argc, char *argv[])
{
  int		option,
  		usage = 0,
		norm = 0,
		useMul = 0;
  WlzRGBAColorChannel denC = WLZ_RGBA_CHANNEL_GREEN,
  		mulC = WLZ_RGBA_CHANNEL_GREY,
		numC = WLZ_RGBA_CHANNEL_RED;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE          *fP = NULL;
  char		*inObjStr,
  		*outObjStr;
  const char    *errMsg;
  static char   optList[] = "d:m:n:o:Nh",
  		outObjStrDef[] = "-",
		inObjStrDef[] = "-";

  outObjStr = outObjStrDef;
  inObjStr = inObjStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'd':
	usage = WlzRGBRatioChan(&denC, optarg);
	break;
      case 'n':
	usage = WlzRGBRatioChan(&numC, optarg);
	break;
      case 'm':
	useMul = 1;
	usage = WlzRGBRatioChan(&mulC, optarg);
	break;
      case 'o':
        outObjStr = optarg;
	break;
      case 'N':
        norm = 1;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if((*inObjStr == '\0') || (*outObjStr == '\0'))
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
        inObjStr = *(argv + optind);
      }
    }
  }
  if(usage == 0)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjStr == NULL) ||
       (*inObjStr == '\0') ||
       ((fP = (strcmp(inObjStr, "-")?
	      fopen(inObjStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s\n",
		     *argv, inObjStr);
    }
    if(fP && strcmp(inObjStr, "-"))
    {
      fclose(fP);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      outObj = WlzRGBChanRatio(inObj, numC, denC, mulC, useMul, norm, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to write output object (%s).\n",
		       argv[0], errMsg);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WLZ_ERR_WRITE_EOF;
      if(((fP = (strcmp(outObjStr, "-")?  fopen(outObjStr, "w"):
					      stdout)) == NULL) ||
	 ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
      {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to write output object.\n",
		       argv[0]);
      }
      if(fP && strcmp(outObjStr, "-"))
      {
	fclose(fP);
      }
    }
  }
  else
  {
    (void )fprintf(stderr,
      "Usage: %s [-h] [-N] [-o<out>] [-d<r|g|b>] [-n<r|g|b>] [<in>]\n"
      "Computes log ratio of RGB channels in an RGBA object for each\n"
      "pixel using: ratio = ln(1 + n/(1 + d)). This results in an\n"
      "object normalised to the range 0-255.\n"
      "Colour channels are specified by the ** letter of the colour\n"
      "channel strings *r*ed, *g*reen, *b*lue, *y*ellow, *m*agenta\n"
      "*c*yan, *h*ue, *s*aturation, *v*alue and gr*e*y, ie r, g, b,\n"
      "y, m, c, h, s, v and e.\n"
      "No multiplier is used unless given on the command line.\n"
      "By default the input object is read from the standard input and the\n"
      "output object is written to standard the output.\n"
      "Version: %s\n"
      "Options:\n"
      "  -N  Output is normalised to range 0-255.\n"
      "  -d  Denominator colour channel.\n"
      "  -n  Numerator colour channel.\n"
      "  -m  Multiplier colour channel.\n"
      "  -o  Output file.\n"
      "  -h  Help, prints this usage message.\n",
      *argv,
      WlzVersion());
  }
  exit(errNum);
}

static int	WlzRGBRatioChan(WlzRGBAColorChannel *chan, char *str)
{
  int		val,
  		usage = 0;

  usage = (str == NULL) || (strlen(str) != 1);
  if(usage == 0)
  {
    val = toupper(*str);
    switch(val)
    {
      case 'E':
        *chan = WLZ_RGBA_CHANNEL_GREY;
	break;
      case 'R':
        *chan = WLZ_RGBA_CHANNEL_RED;
	break;
      case 'G':
        *chan = WLZ_RGBA_CHANNEL_GREEN;
	break;
      case 'B':
        *chan = WLZ_RGBA_CHANNEL_BLUE;
	break;
      case 'C':
        *chan = WLZ_RGBA_CHANNEL_CYAN;
	break;
      case 'M':
        *chan = WLZ_RGBA_CHANNEL_MAGENTA;
	break;
      case 'Y':
        *chan = WLZ_RGBA_CHANNEL_YELLOW;
	break;
      case 'H':
        *chan = WLZ_RGBA_CHANNEL_HUE;
	break;
      case 'S':
        *chan = WLZ_RGBA_CHANNEL_SATURATION;
	break;
      case 'V':
        *chan = WLZ_RGBA_CHANNEL_BRIGHTNESS;
	break;
      default:
        usage = 1;
	break;
    }
  }
  return(usage);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
