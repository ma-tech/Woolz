#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzIntersect_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzIntersect.c
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
* \brief	Computes the intersection of the input objects.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzintersect "WlzIntersect"
*/

/*!
\ingroup BinWlz
\defgroup wlzintersect WlzIntersect
\par Name
WlzIntersect - computes the intersection of the input objects.
\par Synopsis
\verbatim
WlzIntersect [-h] [-n] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Maximum number of objects, default 100.</td>
  </tr>
</table>
\par Description
Computes the intersection of the input objects.
\par Examples
\verbatim
cat obj0.wlz obj1.wlz obj2.wlz | WlzIntersect >int.wlz
\endverbatim
Computes the intersection of the domain of the objects read from the files
obj0.wlz, obj1.wlz and obj2.wlz. The resulting object is written to
int.wlz.
\par File
\ref WlzIntersect.c "WlzIntersect.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzdiffdomain "WlzDiffDomain(1)"
\ref wlzunion "WlzUnion(1)"
\ref WlzIntersectN "WlzIntersectN(3)"
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
	  "Usage:\t%s [-h] [-n] [<input file>]\n"
	  "\tCalculate the intersection of the input woolz objects\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -n#       Maximum number of objects -default=100\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj1, *obj, **objlist;
  WlzObjectType	type = (WlzObjectType) -1;
  int 		n, nmax;
  FILE		*inFile;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		optList[] = "n:h";
  int		option;
  int		uvt = 1;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  nmax = 100;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'n':
      nmax = atoi(optarg);
      if( nmax < 1 ){
	fprintf(stderr, "%s: nmax = %d is invalid\n", argv[0], nmax);
	usage(argv[0]);
	return( 1 );
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
      return( 1 );
    }
  }

  /* allocate space for the object pointers */
  if( (objlist = (WlzObject **)
       AlcMalloc(sizeof(WlzObject *) * nmax)) == NULL ){
    (void )fprintf(stderr, "%s: Memory allocation failed - nmax too high?\n",
    		   argv[0]);
    return( 1 );
  }

  /* read objects accumulating compatible types */
  n = 0;
  while( ((obj = WlzReadObj(inFile, NULL)) != NULL) && (n < nmax) ) {

    if( type == -1 &&
	(obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ) ){
      type = obj->type;
    }

    if( obj->type == type ){
      objlist[n++] = obj;
      if( obj->values.core == NULL ){
	uvt = 0;
      }
    } else {
      WlzFreeObj( obj );
    }
  }

  if( type == (WlzObjectType) -1 ){
    return( WLZ_ERR_NONE );
  }

  if(((obj1 = WlzIntersectN(n, objlist, uvt, &errNum)) != NULL) &&
     (errNum == WLZ_ERR_NONE)) {
    errNum = WlzWriteObj(stdout, obj1);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr, "%s: failed to find intersection (%s).\n",
    		   argv[0], errMsg);
  }
  /* freespace so purify can check for leaks */
  WlzFreeObj(obj1);
  while( n-- ){
    WlzFreeObj(objlist[n]);
  }
  AlcFree((void *) objlist);

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
