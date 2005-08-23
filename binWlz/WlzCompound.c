#pragma ident "MRC HGU $Id$"
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
    <td>Maximum number of objects, default 100.</td>
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

\par See Also
\ref wlzexplode "WlzExplode(1)"

\par Bugs
None known
\author       richard <Richard.Baldock@hgu.mrc.ac.uk>
\date         Wed Jul 27 23:10:10 2005
\version      MRC HGU $Id$
              $Revision$
              $Name$
\par Copyright:
             1994-2003 Medical Research Council, UK.
              All rights reserved.
\par Address:
              MRC Human Genetics Unit,
              Western General Hospital,
              Edinburgh, EH4 2XU, UK.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/***********************************************************************
* Project:      Woolz
* Title:        WlzCompound.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Generates a Woolz compound object from the given
*		input objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
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
	  "\tto stdout. Options are:\n"
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
  while( ((obj = WlzAssignObject(WlzReadObj(inFile, NULL),
  			         NULL)) != NULL) && (n < nmax) ) {

    if( type == -1 &&
	(obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ) ){
      type = obj->type;
    }

    if( obj->type == type ){
      objlist[n++] = WlzAssignObject(obj, NULL);
    } else {
      WlzFreeObj( obj );
    }
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
