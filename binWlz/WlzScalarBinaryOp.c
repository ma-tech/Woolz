#if defined(__GNUC__)
#ident "MRC HGU $Id:"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id:"
#else static char _Wlz_calarBinaryOp_c[] = "MRC HGU $Id:";
#endif
#endif
/*!
* \file         WlzScalarBinaryOp.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Oct 20 09:57:05 2006
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par Copyright:
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
* \ingroup      WlzBinaryOps
* \brief        Apply a scalar operation to a woolz object.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/
 
/*!
\ingroup      WlzBinaryOps
\defgroup     wlzscalarbinaryop WlzScalarBinaryOp
\par Name
WlzScalarBinaryOp - 
\par Synopsis
\verbatim
WlzScalarBinaryOp

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

/* insert code here */
/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-b <op>] [-h] [-v] [<value>] [<input file>]\n"
	  "\tApply a scalar binary operation to the grey values of\n"
	  "\tthe given image with respect to the input value. Note\n"
	  "\tmost options require an input value.\n"
	  "\tOptions are:\n"
	  "\t  -b <op>  convert to a given channel\n"
	  "\t           op =  %d - Add\n"   
	  "\t              =  %d - SUBTRACT\n"
	  "\t              =  %d - MULTIPLY\n"
	  "\t              =  %d - DIVIDE\n"
	  "\t              =  %d - MODULUS\n"
	  "\t              =  %d - EQ\n"
	  "\t              =  %d - NE\n"
	  "\t              =  %d - GT\n"
	  "\t              =  %d - GE\n"
	  "\t              =  %d - LT\n"
	  "\t              =  %d - LE\n"
	  "\t              =  %d - AND\n"
	  "\t              =  %d - OR\n"
	  "\t              =  %d - XOR\n"
	  "\t              =  %d - MAX\n"
	  "\t              =  %d - MIN\n"
	  "\t              =  %d - MAGNITUDE\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str,
	  WLZ_BO_ADD,
	  WLZ_BO_SUBTRACT,
	  WLZ_BO_MULTIPLY,
	  WLZ_BO_DIVIDE,
	  WLZ_BO_MODULUS,
	  WLZ_BO_EQ,
	  WLZ_BO_NE,
	  WLZ_BO_GT,
	  WLZ_BO_GE,
	  WLZ_BO_LT,
	  WLZ_BO_LE,
	  WLZ_BO_AND,
	  WLZ_BO_OR,
	  WLZ_BO_XOR,
	  WLZ_BO_MAX,
	  WLZ_BO_MIN,
	  WLZ_BO_MAGNITUDE);
    return;
}

int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *newobj;
  FILE		*inFile;
  char 		optList[] = "b:hv";
  int		option;
  int		verboseFlg=0;
  WlzBinaryOperatorType operator = WLZ_BO_ADD;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'b':
      operator = atoi(optarg);
      switch( operator ){
      case WLZ_BO_ADD:
      case WLZ_BO_SUBTRACT:
      case WLZ_BO_MULTIPLY:
      case WLZ_BO_DIVIDE:
      case WLZ_BO_MODULUS:
      case WLZ_BO_EQ:
      case WLZ_BO_NE:
      case WLZ_BO_GT:
      case WLZ_BO_GE:
      case WLZ_BO_LT:
      case WLZ_BO_LE:
      case WLZ_BO_AND:
      case WLZ_BO_OR:
      case WLZ_BO_XOR:
      case WLZ_BO_MAX:
      case WLZ_BO_MIN:
      case WLZ_BO_MAGNITUDE:
	break;

      default:
	usage(argv[0]);
	return( 1 );
      }
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );
    }
  }

  return 0;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
