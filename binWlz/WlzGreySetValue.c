#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreySetValue_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzGreySetValue.c
* \author       Richard Baldock
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
* \brief	Sets the grey values of the object to a specified value.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzgreysetvalue "WlzGreySetValue"
*/

/*!
\ingroup BinWlz
\defgroup wlzgreysetvalue WlzGreySetValue
\par Name
WlzGreySetValue - sets the grey values of the object to a specified value.
\par Synopsis
\verbatim
WlzGreySetValue [-g#] [-h] [-v] [<input mask> [<input obj>]]
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
    <td><b>-g</b></td>
    <td>The specified grey value - default 0.</td>
  </tr>
</table>
\par Description
Set the grey values of the object to the input value.
A valuetable will be attached if required.
\par Examples
\verbatim
WlzGreySetValue -g 128 dom.wlz >out.wlz
\endverbatim
Creates a new object with the domain of the object read from dom.wlz
but with all grey values having value 128.
The resulting object is written to out.wlz.
\par File
\ref WlzGreySetValue.c "WlzGreySetValue.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgreymask "WlzGreyMask(1)"
\ref WlzGreySetValue "WlzGreySetValue(3)"
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
	  "Usage:\t%s [-g#] [-h] [-v] [<input mask> [<input obj>]]\n"
	  "\tSet the grey values of the object to the input value.\n"
	  "\tA valuetable will be attached if required.\n"
	  "\tOptions are:\n"
	  "\t  -g#       the new grey value - default 0\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *tmpObj;
  WlzDomain	*domains;
  WlzValues	values, *valuess;
  FILE		*inFile;
  char 		optList[] = "g:hv";
  int		option;
  WlzPixelV	greyVal;
  WlzPixelV	bckgrnd;
  WlzObjectType	type;
  int		verboseFlg=0, p;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  greyVal.type = WLZ_GREY_FLOAT;
  greyVal.v.flv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'g':
      greyVal.v.flv = atof(optarg);
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

  /* check for read from file */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* set up type and background */
  type = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_UBYTE, NULL);
  bckgrnd.type = WLZ_GREY_UBYTE;
  bckgrnd.v.ubv = 0;

  /* read objects and set value if possible */
  while(((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) &&
        (errNum == WLZ_ERR_NONE))
  {
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
      if( obj->values.core == NULL ){
	values.v = WlzNewValueTb(obj, type, bckgrnd, &errNum);
	obj->values = WlzAssignValues(values, NULL);
      }
      errNum = WlzGreySetValue(obj, greyVal);
      break;

    case WLZ_3D_DOMAINOBJ:
      switch( obj->domain.p->type ){
      case WLZ_PLANEDOMAIN_DOMAIN:
	if( obj->values.core == NULL ){
	  if( values.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					       obj->domain.p->plane1,
					       obj->domain.p->lastpl,
					       bckgrnd, NULL, &errNum) ){
	    obj->values = WlzAssignValues(values, NULL);
	    domains = obj->domain.p->domains;
	    valuess = obj->values.vox->values;
	    for(p=0; p < (obj->domain.p->lastpl - obj->domain.p->plane1 + 1);
		p++, domains++, valuess++){
	      if( (*domains).core ){
		values.core = NULL;
		tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, values,
				     NULL, NULL, NULL);
		values.v = WlzNewValueTb(tmpObj, type, bckgrnd, &errNum);
		*valuess = WlzAssignValues(values, NULL);
		WlzFreeObj(tmpObj);
	      }
	    }
	  }
	}
	errNum = WlzGreySetValue(obj, greyVal);
	break;

      default:
	break;
      }
      break;

    default:
      break;
    }

    if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write object (%s).\n",
		     argv[0], errMsg);
      return(1);
    }
    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
