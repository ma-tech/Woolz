#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzGreyTransfer.c
* \author       Richard Baldock
* \date         January 2002
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
* \brief	Copies grey values from a source object to a destination
* 		object within the intersection of the object's domains.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzgreytransfer "WlzGreyTransfer"
*/

/*!
\ingroup BinWlz
\defgroup wlzgreytransfer WlzGreyTransfer
\par Name
WlzGreyTransfer - copies grey values from a source object to a destination
object within the intersection of the object's domains.
\par Synopsis
\verbatim
WlzGreyTransfer [-h] [-v] [<source image> [<destination obj>]]
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
</table>
\par Description
Copy grey values from the source object to the destination
object within the domain of intersection. The output object
has the domain and grey-values of the destination object
otherwise. If both objects are read from stdin then the
first object is the destination and the second the source.
\par Examples
\verbatim
WlzGreyTransfer obj0.wlz obj1.wlz >out.wlz
Creates a new object with the domain of obj1.wlz
in which the values are those of obj0.wlz within the intersection of the
two domains and obj1.wlz elsewhere.
\endverbatim
\par File
\ref WlzGreyTransfer.c "WlzGreyTransfer.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgreytemplate "WlzGreyTemplate(1)"
\ref WlzGreyTransfer "WlzGreyTransfer(1)"
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
	  "Usage:\t%s [-h] [-v] [<source image> [<destination obj>]]\n"
	  "\tCopy grey values from the source object to the destination\n"
	  "\tobject within the domain of intersection. The output object\n"
	  "\thas the domain and grey-values of the destination object \n"
	  "\totherwise. If both objects are read from stdin then the\n"
	  "\tfirst object is the destination and the second the source.\n"
	  "\tOptions are:\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *tmpl, *newobj;
  FILE		*inFile;
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* check for read from file, note if only one file is defined then
     it is the template and the object should be on stdin */
  if( (argc - optind) >= 2 ){
    /* read the template */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (tmpl = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read template from file %s\n", argv[0],
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
    /* read the template */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (tmpl = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read template from file %s\n", argv[0],
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
    inFile = stdin;
    if( (obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read obj from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
    if( (tmpl = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read template from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
    
  /* apply template and write resultant object */
  if( newobj = WlzGreyTransfer(obj, tmpl, &errNum) ){
    if((errNum = WlzWriteObj(stdout, newobj)) != WLZ_ERR_NONE) {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write object (%s).\n",
		     argv[0], errMsg);
      return(1);
    }
    WlzFreeObj(newobj);
  }
  WlzFreeObj(tmpl);
  WlzFreeObj(obj);

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
