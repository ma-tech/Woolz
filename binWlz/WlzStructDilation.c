#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzStructDilation_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzStructDilation.c
* \author       Richard Baldock
* \date         May 2001
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
* \brief	Dilates an object using a structuring element object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzstructdilation "WlzStructDilation"
*/

/*!
\ingroup BinWlz
\defgroup wlzstructdilation WlzStructDilation
\par Name
WlzStructDilation - dilates an object using a structuring element object.
\par Synopsis
\verbatim
WlzStructDilation [-h] [-s <struct-elem file>] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Structuring element object.</td>
  </tr>
</table>
\par Description
Dilates a domain object using a structuring element writing the new object
to the standard output.
\par Examples
\verbatim
WlzStructDilation -s str.wlz in.wlz >out.wlz
\endverbatim
Dilates the object read from in.wlz using a structuring object read from
str.wlz. The resulting object is written to out.wlz.
\par File
\ref WlzStructDilation.c "WlzStructDilation.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzerosion "WlzErosion(1)"
\ref wlzdilation "WlzDilation(1)"
\ref wlzdomainfill "WlzDomainFill(1)"
\ref wlzstructerosion WlzStructErosion(1)"
\ref WlzStructDilation "WlzStructDilation(3)"
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
	  "Usage:\t%s [-s <struct-elem file>] [-h] [<input file>]\n"
	  "\tDilate a domain woolz object(s) using given\n"
	  "\tstructuring element writing the new object(s)\n"
	  "\tto standard output\n"
	  "\tOptions are:\n"
	  "\t  -s <file> file holding the structuring\n"
	  "\t            element object. Default is the\n"
	  "\t            first object read in from the input\n"
	  "\t            file or standard input.\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *seObj, *nobj;
  FILE		*inFile;
  FILE		*seFile=NULL;
  char 		optList[] = "hs:";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 's':
      if( (seFile = fopen(optarg, "r")) == NULL ){
	fprintf(stderr, "%s: can't open structuring element file: %s\n",
		argv[0], optarg);
	usage(argv[0]);
	return WLZ_ERR_UNSPECIFIED;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }
  if( !seFile ){
    seFile = inFile;
  }
  
  if((seObj = WlzAssignObject(WlzReadObj(seFile, &errNum), NULL)) != NULL){

    /* read objects */
    while((obj = WlzAssignObject(WlzReadObj(inFile, &errNum), NULL)) != NULL) 
    {
      if((nobj = WlzStructDilation(obj, seObj, &errNum)) != NULL){
	nobj = WlzAssignObject(nobj, NULL);
	WlzWriteObj(stdout, nobj);
	WlzFreeObj(nobj);
      }
      WlzFreeObj(obj);
    }
    WlzFreeObj(seObj);
  }
  
  /* trap the WLZ_ERR_READ_EOF since this is a legal way of indicating
     the end of objects in a file */
  if( errNum == WLZ_ERR_READ_EOF ){
    errNum = WLZ_ERR_NONE;
  }
  return errNum;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
