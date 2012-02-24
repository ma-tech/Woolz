#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzEMAPDomainTransform_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzEMAPDomainTransform.c
* \author       Richard Baldock
* \date         December 2006
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
* \ingroup      BinWlzApp
* \brief        Apply a warp transform to an EMAP model domain to
*               convert to a linked model.
*               
* \par Binary
* \ref wlzemapdomaintransform "WlzEMAPDomainTransform"
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzemapdomaintransform WlzEMAPDomainTransform
\par Name
WlzEMAPDomainTransform - 
\par Synopsis
\verbatim
WlzEMAPDomainTransform

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-</b></td>
    <td> </td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
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

\par Description

\par Examples
\verbatim
\endverbatim

\par See Also
\par Bugs
None known
*/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>

#include <Wlz.h>
#include <WlzEMAP.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */


static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s -m <src> -M <dest> [-h] [-v]"
	  " <object input file>\n"
	  "\t  -h                 Help - this message\n"
	  "\t  -v                 verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{
  char 		optList[] = "m:M:hv";
  int		option;
  FILE		*inFP;
  WlzObject	*inObj, *outObj;
  int		verboseFlg=0;
  char	      	*srcModel, *dstModel;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* defaults */
  srcModel = NULL;
  dstModel = NULL;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'm':
      srcModel = AlcStrDup(optarg);
      break;

    case 'M':
      dstModel = AlcStrDup(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 0;

    case 'v':
      verboseFlg = 1;
      break;

    }
  }
  if((srcModel == NULL) || (dstModel == NULL)){
    fprintf(stderr, "%s: source and destination models must be defined\n",
	    argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* read the domain to be transformed */
  inFP = stdin;
  if( optind < argc ){
    if( (inFP = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }
  if( ((inObj = WlzReadObj(inFP, &errNum)) != NULL) ){
    outObj = WlzEMAPDomainTransform(srcModel, dstModel, inObj, &errNum);
    if( outObj ){
      WlzWriteObj(stdout, outObj);
      WlzFreeObj(outObj);
    }
    else {
      fprintf(stderr, "%s: failed to transform object\n", argv[0]);
      return 1;
    }
    WlzFreeObj(inObj);
  }
  

  return 0;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
