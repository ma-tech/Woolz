#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRGBAConvert_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzRGBAConvert.c
* \author       Richard Baldock
* \date         June 2004
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
* \brief	Colour conversion for an RGBA colour object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzrgbaconvert "WlzRGBAConvert"
*/

/*!
\ingroup BinWlz
\defgroup wlzrgbaconvert WlzRGBAConvert
\par Name
WlzRGBAConvert - colour conversion for an RGBA colour object.
\par Synopsis
\verbatim
WlzRGBAConvert [-c] [-m] [-s#] [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-c</b></td>
    <td>Convert to compound, default.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Convert to modulus.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Select colour space
    <table width="500" border="0">
    <tr> <td>1</td> <td>RGB, default.</td> </tr>
    <tr> <td>2</td> <td>HSB.</td> </tr>
    <tr> <td>3</td> <td>CMY.</td> </tr>
    </table>
    </td>
  </tr>
</table>
\par Description
Converts a RGBA (colour) object to either a compound object
or a modulus of the RGB values,
writing the new object to standard output.
The input object <b>must</b> have grey-value type RGBA.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzRGBAConvert.c "WlzRGBAConvert.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzconvertpix "WlzConvertPix(1)"
\ref WlzRGBAToCompound "WlzRGBAToCompound(3)"
\ref WlzRGBAToModulus "WlzRGBAToModulus(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <Wlz.h>
				
/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-c] [-C <channel>] [-m] [-s#][-h] [<input file>]\n"
	  "\tConvert the RGBA woolz object to a compound object\n"
	  "\tor to a specific channel, including modulus, of the rgb\n"
	  "\tvalues, writing the new object to standard output.\n"
	  "\tNote input object MUST have grey-value type RGBA\n."
	  "\tOptions are:\n"
	  "\t  -c       convert to compound (default)\n"
	  "\t  -C #     convert to a given channel\n"
	  "\t           # = %d: red\n"
	  "\t               %d: green\n"
	  "\t               %d: blue\n"
	  "\t               %d: hue\n"
	  "\t               %d: saturation\n"
	  "\t               %d: brightness\n"
	  "\t               %d: cyan\n"
	  "\t               %d: magenta\n"
	  "\t               %d: yellow\n"
	  "\t  -m       convert to modulus\n"
	  "\t  -s#	select colour space:\n"
	  "\t           # = %d: RGB (default)\n"
	  "\t               %d: HSB\n"
	  "\t               %d: CMY\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str,
	  WLZ_RGBA_CHANNEL_RED,
	  WLZ_RGBA_CHANNEL_GREEN,
	  WLZ_RGBA_CHANNEL_BLUE,
	  WLZ_RGBA_CHANNEL_HUE,
	  WLZ_RGBA_CHANNEL_SATURATION,
	  WLZ_RGBA_CHANNEL_BRIGHTNESS,
	  WLZ_RGBA_CHANNEL_CYAN,
	  WLZ_RGBA_CHANNEL_MAGENTA,
	  WLZ_RGBA_CHANNEL_YELLOW,
	  WLZ_RGBA_SPACE_RGB,
	  WLZ_RGBA_SPACE_HSB,
	  WLZ_RGBA_SPACE_CMY);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *newobj;
  FILE		*inFile;
  char 		optList[] = "cC:ms:hv";
  int		option;
  int		compoundFlg=1;
  WlzRGBAColorSpace	colSpc=WLZ_RGBA_SPACE_RGB;
  WlzRGBAColorChannel	chan=WLZ_RGBA_CHANNEL_GREY;
  WlzErrorNum	errNum;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'c':
      compoundFlg = 1;
      break;

    case 'C':
      switch( chan = (WlzRGBAColorChannel) atoi(optarg) ){
      case WLZ_RGBA_CHANNEL_RED:
      case WLZ_RGBA_CHANNEL_GREEN:
      case WLZ_RGBA_CHANNEL_BLUE:
      case WLZ_RGBA_CHANNEL_HUE:
      case WLZ_RGBA_CHANNEL_SATURATION:
      case WLZ_RGBA_CHANNEL_BRIGHTNESS:
      case WLZ_RGBA_CHANNEL_CYAN:
      case WLZ_RGBA_CHANNEL_MAGENTA:
      case WLZ_RGBA_CHANNEL_YELLOW:
	break;

      default:
        fprintf(stderr, "%s: colour channel = %d is invalid\n",
		argv[0], chan );
        usage(argv[0]);
        return 1;

      }
      compoundFlg = 0;
      break;

    case 'm':
      compoundFlg = 0;
      break;

    case 's':
      switch( colSpc = (WlzRGBAColorSpace) atoi(optarg) ){

      case WLZ_RGBA_SPACE_RGB:
      case WLZ_RGBA_SPACE_HSB:
      case WLZ_RGBA_SPACE_CMY:
	break;

      default:
        fprintf(stderr, "%s: colour space = %d is invalid\n",
		argv[0], colSpc );
        usage(argv[0]);
        return 1;

      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* read objects and convert if possible */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL ){
    if( obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ ){
      if( compoundFlg ){
	if( (newobj = (WlzObject *) 
	     WlzRGBAToCompound(obj, colSpc, &errNum)) != NULL ){
	  (void )WlzWriteObj(stdout, newobj);
	  WlzFreeObj(newobj);
	}
	else {
	  if( errNum == WLZ_ERR_VALUES_TYPE ){
	    fprintf(stderr, "%s: wrong values type, object not converted\n",
		    argv[0]);
	    (void )WlzWriteObj(stdout, obj);
	  }
	  else {
	    fprintf(stderr, "%s: something wrong, time to quit.\n",
		    argv[0]);
	    usage(argv[0]);
	    return 1;
	  }
	}
      }
      else {
	if( chan == WLZ_RGBA_CHANNEL_GREY ){
	  if( (newobj = WlzRGBAToModulus(obj, &errNum)) != NULL ){
	    (void )WlzWriteObj(stdout, newobj);
	    WlzFreeObj(newobj);
	  }
	  else {
	    if( errNum == WLZ_ERR_VALUES_TYPE ){
	      fprintf(stderr, "%s: wrong values type, object not converted\n",
		      argv[0]);
	      (void )WlzWriteObj(stdout, obj);
	    }
	    else {
	      fprintf(stderr, "%s: something wrong, time to quit.\n",
		      argv[0]);
	      usage(argv[0]);
	      return 1;
	    }
	  }
	}
	else {
	  if( (newobj = WlzRGBAToChannel(obj, chan, &errNum)) != NULL ){
	    (void )WlzWriteObj(stdout, newobj);
	    WlzFreeObj(newobj);
	  }
	  else {
	    if( errNum == WLZ_ERR_VALUES_TYPE ){
	      fprintf(stderr, "%s: wrong values type, object not converted\n",
		      argv[0]);
	      (void )WlzWriteObj(stdout, obj);
	    }
	    else {
	      fprintf(stderr, "%s: something wrong, time to quit.\n",
		      argv[0]);
	      usage(argv[0]);
	      return 1;
	    }
	  }
	}
      }
    }
    else {
      (void )WlzWriteObj(stdout, obj);
    }

    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
