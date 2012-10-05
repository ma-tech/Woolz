#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreySetValue_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzGreySetValue.c
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
* \brief	Sets the grey values of the object to a specified value.
* \ingroup	BinWlz
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
WlzGreySetValue [-c #,#,#] [-g#] [-h] [<input mask> [<input obj>]]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-c</b></td>
    <td>New colour value r,g,b - default 0,0,0.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>The specified grey value - default 0.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
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
The grey type of the resulting object are chosen to be the minimum
which can be used to encode the specified value,
with the possible types (in increasing order) being:
unsigned byte, short, int and double.
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
#include <limits.h>
#include <float.h>
#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  (void )fprintf(stderr,
	  "Usage:\t%s [-c#,#,#] [-g#] [-h]\n"
	  "\t  [<input mask> [<input obj>]]\n"
	  "\tSet the grey values of the object to the input value.\n"
	  "\tA valuetable will be attached if required.\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -c#,#,#   the new colour value r,g,b - default 0,0,0\n"
	  "\t  -g#       the new grey value - default 0\n"
	  "\t  -h        help - prints this usage message\n",
	  proc_str,
	  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *tmpObj;
  WlzDomain	*domains;
  WlzValues	values, *valuess;
  FILE		*inFile;
  char 		optList[] = "c:g:h";
  int		igv,
		pCnt,
  		option;
  WlzPixelV	greyVal;
  WlzPixelV	bckgrnd;
  WlzObjectType	type;
  int		p, colFlg=0;
  int		red=0, green=0, blue=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  greyVal.type = WLZ_GREY_DOUBLE;
  greyVal.v.dbv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'c':
      if(sscanf(optarg, "%d,%d,%d", &red, &green, &blue) != 3)
      {
        usage(argv[0]);
	return 1;
      }
      WLZ_RGBA_RGBA_SET(greyVal.v.rgbv, red, green, blue, 255);
      greyVal.type = WLZ_GREY_RGBA;
      colFlg = 1;
      break;

    case 'g':
      if(sscanf(optarg, "%lg", &(greyVal.v.dbv)) != 1)
      {
        usage(argv[0]);
	return 1;
      }
      colFlg = 0;
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

  /* Set up the type and background: Type will be the minimum required from
     unsigned byte, short, int and double. */
  if( colFlg ){
    bckgrnd.v.rgbv = 0;
    bckgrnd.type = WLZ_GREY_RGBA;
  }
  else {
    igv = WLZ_NINT(greyVal.v.dbv);
    bckgrnd.type = WLZ_GREY_DOUBLE;
    bckgrnd.v.dbv = 0.0;
    if((fabs(greyVal.v.dbv - igv) < DBL_EPSILON) &&
       (igv >= INT_MIN) && (igv <= INT_MAX))
    {
      if((igv >= 0) && (igv <= 255))
      {
	bckgrnd.v.ubv = 0;
	bckgrnd.type = WLZ_GREY_UBYTE;
      }
      else if((igv >= SHRT_MIN) && (igv <= SHRT_MAX))
      {
	bckgrnd.v.shv = 0;
	bckgrnd.type = WLZ_GREY_SHORT;
      }
      else
      {
	bckgrnd.v.inv = 0;
	bckgrnd.type = WLZ_GREY_INT;
      }
      /* WlzValueConvertPixel() sets greyVal.type. */
      (void )WlzValueConvertPixel(&greyVal, greyVal, bckgrnd.type);
    }
  }
  type = WlzGreyTableType(WLZ_GREY_TAB_RAGR, greyVal.type, NULL);

  /* read objects and set value if possible */
  while(((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) &&
        (errNum == WLZ_ERR_NONE))
  {
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
      if((obj->values.core == NULL) || (obj->values.core->type != type))
      {
	if(obj->values.core)
	{
	  (void )WlzFreeValueTb(obj->values.v);
	  obj->values.core = NULL;
	}
	values.v = WlzNewValueTb(obj, type, bckgrnd, &errNum);
	obj->values = WlzAssignValues(values, NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzGreySetValue(obj, greyVal);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      switch( obj->domain.p->type ){
      case WLZ_PLANEDOMAIN_DOMAIN:
	pCnt = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
        if(obj->values.core)
	{
	  domains = obj->domain.p->domains;
	  valuess = obj->values.vox->values;
	  p = 0;
	  while((errNum == WLZ_ERR_NONE) && (p < pCnt))
	  {
	    if((*domains).core)
	    {
	      if(((*valuess).core == NULL) || (valuess->core->type != type))
	      {
		values.core = NULL;
		tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, values,
				     NULL, NULL, NULL);
		values.v = WlzNewValueTb(tmpObj, type, bckgrnd, &errNum);
		if((*valuess).core)
		{
		  (void )WlzFreeValueTb(valuess->v);
		}
		*valuess = WlzAssignValues(values, NULL);
		WlzFreeObj(tmpObj);
	      }
	    }
	    ++p;
	    ++valuess;
	    ++domains;
	  }
	}
	else
	{
	  if((values.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					       obj->domain.p->plane1,
					       obj->domain.p->lastpl,
					       bckgrnd, NULL,
					       &errNum)) != NULL){
	    obj->values = WlzAssignValues(values, NULL);
	    domains = obj->domain.p->domains;
	    valuess = obj->values.vox->values;
	    p = 0;
	    while((errNum == WLZ_ERR_NONE) && (p < pCnt))
	    {
	      if((*domains).core)
	      {
		values.core = NULL;
		tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, values,
				     NULL, NULL, NULL);
		values.v = WlzNewValueTb(tmpObj, type, bckgrnd, &errNum);
		*valuess = WlzAssignValues(values, NULL);
		WlzFreeObj(tmpObj);
	      }
	      ++p;
	      ++valuess;
	      ++domains;
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzGreySetValue(obj, greyVal);
	}
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
