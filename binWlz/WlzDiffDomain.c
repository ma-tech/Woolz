#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzDiffDomain_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzDiffDomain.c
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
* \brief	Computes the difference between two domains.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzdiffdomain "WlzDiffDomain"
*/

/*!
\ingroup BinWlz
\defgroup wlzdiffdomain WlzDiffDomain
\par Name
WlzDiffDomain  -  computes the difference between two domains.
\par Synopsis
\verbatim
WlzDiffDomain  [-h]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Calculate the difference between the domains of the two input objects.
\par Examples
\verbatim
cat obj1.wlz obj2.wlz | WlzDiffDomain >diff.wlz
\endverbatim
Computes the difference between the domains of objects obj1.wlz
and obj2.wlz, with the resulting object being written to diff.wlz.
\par File
\ref WlzDiffDomain.c "WlzDiffDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzunion "WlzUnion(1)"
\ref wlzintersect "WlzIntersect(1)"
\ref WlzDiffDomain "WlzDiffDomain(3)"
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
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tCalculate the difference of the input woolz objects\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *obj1, *obj2;
  WlzObjectType	type = (WlzObjectType) -1;
  int 		n;
  FILE		*inFile;
  char 		optList[] = "h";
  int		option;
  WlzErrorNum	errNum;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case '~': /* dummy to avoid compiler warning */
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

  /* read objects accepting the first two compatible types */
  n = 0;
  while( ((obj = WlzReadObj(inFile, NULL)) != NULL) && (n < 2) ) {

    if( type == -1 &&
	(obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ) ){
      type = obj->type;
    }

    if( obj->type != type ){
      WlzFreeObj( obj );
      continue;
    }

    n++;
    if( n == 1 ){
      obj1 = obj;
    }
    else {
      obj2 = obj;
    }

  }

  if( type == (WlzObjectType) -1 ){
    fprintf(stderr, "%s: no domain objects input\n", argv[0]);
    return( WLZ_ERR_NONE );
  }

  if( n < 2 ){
    fprintf(stderr, "%s: insufficient objects of type %d\n", argv[0], type);
    return( WLZ_ERR_NONE );
  }
    
  obj = WlzDiffDomain(obj1, obj2, &errNum);
  errNum = WlzWriteObj(stdout, obj);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    fprintf(stderr, "%s: failed to write output object (%s)\n",
    	    argv[0],
	    errMsg);

  }

  /* freespace so purify can check for leaks */
  WlzFreeObj(obj1);
  WlzFreeObj(obj2);
  WlzFreeObj(obj);

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
