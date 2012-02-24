#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyMask_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzGreyMask.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Sets all object values within a mask domain to a mask value.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzgreymask "WlzGreyMask"
*/

/*!
\ingroup BinWlz
\defgroup wlzgreymask WlzGreyMask
\par Name
WlzGreyMask - sets all object values within a mask domain to a mask value.
\par Synopsis
\verbatim
WlzGreyMask [-m#] [-h] [-v] [<input mask> [<input obj>]]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Mask value, default 0.</td>
  </tr>
</table>
\par Description
Apply the mask to the given object setting all pixels
within the mask to the mask value. The object must be
a woolz grey-level object, the mask can be a domain object
polyline or boundary list which will be filled to make
a domain. The mask is always the second object if read
from the standard input.
\par Examples
\verbatim
WlzGreyMask -m 64 dom.wlz obj.wlz >masked.wlz
\endverbatim
Creates a new object mask.wlz which has the same domain and values of
obj.wlz, except where the domain of obj.wlz and mask.wlz intersect
where the values are set to 64.
\par File
\ref WlzGreyMask.c "WlzGreyMask.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgreytemplate "WlzGreyTemplate(1)"
\ref WlzGreyMask "WlzGreyMask(3)"
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
	  "Usage:\t%s [-m#] [-h] [-v] [<input mask> [<input obj>]]\n"
	  "\tApply the mask to the given object setting all pixels\n"
	  "\twithin the mask to the mask value. The object must be\n"
	  "\t a woolz grey-level object, the mask can be a domain object,\n"
	  "\tpolyline or boundary list which will be filled to make\n"
	  "\t a domain. The mask is always the second object if read\n"
	  "\tfrom the standard input\n"
	  "\tOptions are:\n"
	  "\t  -m#       the mask value - default 0\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *mask, *newobj;
  FILE		*inFile = NULL;
  char 		optList[] = "m:hv";
  int		option;
  WlzPixelV	maskVal;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  maskVal.type = WLZ_GREY_FLOAT;
  maskVal.v.flv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'm':
      maskVal.v.flv = atof(optarg);
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* check for read from file, note if one file is defined then
     it is the mask and the object should be on stdin */
  if( (argc - optind) >= 2 ){
    /* read the mask */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (mask = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read mask from file %s\n", argv[0],
	      *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
    optind++;

    /* read the object */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read obj from file %s\n", argv[0],
	      *(argv+optind));

      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
  }
  else if( (argc - optind) == 1 ){
    /* read the mask */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (mask = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read mask from file %s\n", argv[0],
	      *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
    inFile = stdin;
    if( (obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read obj from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    /* else read objects from stdin */
    if( (obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read obj from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
    if( (mask = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read mask from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
    
  /* apply mask and write resultant object */
  if((newobj = WlzGreyMask(obj, mask, maskVal, &errNum)) != NULL){
    if((errNum = WlzWriteObj(stdout, newobj)) != WLZ_ERR_NONE) {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write object (%s).\n",
		     argv[0], errMsg);
      return(1);
    }
    WlzFreeObj(newobj);
  }
  WlzFreeObj(mask);
  WlzFreeObj(obj);

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
