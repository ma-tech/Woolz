#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCompound_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCompound.c
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
* \brief        Generates a Woolz compound object from objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcompound "WlzCompound"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzcompound WlzCompound
\par Name
WlzCompound - Create a Woolx compounf object.
\par Synopsis
\verbatim
WlzCompound [-n <num objects>] [-h]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-n</b></td>
    <td>Maximum number of objects, default 1024.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>
Reads woolz objects from standard input and writes a compound
 object to standard out.

\par Description
Create a woolz compound object from a given input stream. Only
 objects of the same type as the first object read in are included
 and the compound object will be of type WLZ_COMPOUND_ARR_2.

\par Examples
\verbatim
#
# read up to 150 woolz objects in the current directory
# and create a compound object 
#
cat *.wlz | WlzCompound -n 150 > compound_obj.wlz

\endverbatim

\par File
\ref WlzCompound.c "WlzCompound.c"
\par See Also
\ref wlzexplode "WlzExplode(1)"
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
	  "Usage:\t%s [-h] [-n#]\n"
	  "\tGenerate a woolz compound object from the given\n"
	  "\tinput objects. Objects are read from stdin and written\n"
	  "\tto stdout.\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -n#       Maximum number of objects -default=1024\n",
	  proc_str,
	  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj1, *obj, **objlist;
  WlzCompoundArray	*cobj;
  WlzObjectType	type = (WlzObjectType) -1;
  int 		n, nmax;
  FILE		*inFile;
  char 		optList[] = "n:h";
  int		option;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

    
  /* read the argument list and check for an input file */
  opterr = 0;
  nmax = 1024;
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

  /* reading from stdin */
  inFile = stdin;

  /* allocate space for the object pointers */
  if( (objlist = (WlzObject **)
       AlcMalloc(sizeof(WlzObject *) * nmax)) == NULL ){
    (void )fprintf(stderr, "%s: memory allocation failed.\n",
    		   argv[0]);
    return( 1 );
  }

  /* read objects accumulating compatible types */
  n = 0;
  while( ((obj = WlzAssignObject(
                 WlzReadObj(inFile, NULL), NULL)) != NULL) && (n < nmax) ) {

    if( type == -1 &&
	(obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ) ){
      type = obj->type;
    }

    if( obj->type == type || obj->type == WLZ_EMPTY_OBJ ){
      objlist[n++] = WlzAssignObject(obj, NULL);
    }
    (void )WlzFreeObj( obj );
  }

  if( type == WLZ_EMPTY_OBJ ){
    return( 0 );
  }

  if( (cobj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_2,
				   3, n, objlist, WLZ_NULL,
				   &errNum)) == NULL ){
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr, "%s: failed to create compon object (%s).\n",
    		   argv[0], errMsg);
    return(1);
  }
  else {
    obj1 = (WlzObject *) cobj;
    if((errNum = WlzWriteObj(stdout, obj1)) != WLZ_ERR_NONE) {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to write compound object (%s).\n",
      		     argv[0], errMsg);
    }
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
