#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzFacts_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzFacts.c
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
* \brief	Gives information about Woolz objects.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzfacts "WlzFacts"
*/

/*!
\ingroup BinWlz
\defgroup wlzfacts WlzFacts
\par Name
WlzFacts  -  print  out  information about the input woolz
       objects.
\par Synopsis
\verbatim
WlzFacts  [-h] [-f] [-m] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-m</b></td>
    <td>many facts printed out.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Few facts printed out (default).</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
WlzFacts reads woolz objects from the  standard  input  or
       the  given  input  file  and  writes out information about
       their structure (types,  domain  bounds,  area,  etc.)  to
       stderr.
\par Examples
\verbatim
WlzFacts toucan.wlz
\endverbatim
Prints information about the object read from toucan.wlz
to the standard error output.
\par File
\ref WlzFacts.c "WlzFacts.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgreystats "WlzGreyStats(3)"
\ref WlzObjectFacts "WlzObjectFacts(3)"
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
	  "Usage:\t%s [-h] [-f] [-m] [<input file>]\n"
	  "\tPrint facts about the input woolz objects\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -f        Few facts (default)\n"
	  "\t  -m        Many facts - print out more detailed information\n"
	  "",
	  proc_str);
  return;
}
 
main(
     int		argc,
     char**	argv)
{
  WlzObject	*obj;
  FILE	*inFile;
  char 	optList[] = "fmh";
  int		option;
  int		manyFactsFlg=0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  inFile = stdin;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'f':
      manyFactsFlg=0;
      break;

    case 'm':
      manyFactsFlg=1;
      break;

    case 'h':
      usage(argv[0]);
      return( WLZ_ERR_NONE );

    default:
      usage(argv[0]);
      return( WLZ_ERR_PARAM_TYPE );

    }
  }
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( WLZ_ERR_PARAM_DATA );
    }
  }

  /* read objects until EOF or an error */
  while( obj = WlzReadObj(inFile, &errNum) ){
    (void ) WlzObjectFacts(obj, stderr, NULL, manyFactsFlg);
    WlzFreeObj(obj);
    printf("\n");
  }


  /* trap the WLZ_ERR_READ_EOF since this is a legal way of indicating
     the end of objects in a file */
  if(errNum == WLZ_ERR_READ_EOF){
    errNum = WLZ_ERR_NONE;
  }
  return errNum;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
