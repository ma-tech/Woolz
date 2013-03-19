#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstGrayCode_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstGrayCode.c
* \author       Bill Hill
* \date         September 2011
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
* \brief	Simple test for AlgGreyCode(), AlgGreyCodeInv()
* 		and and their long long equivalents.
* \ingroup	binAlgTst
*/
#include <float.h>
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

typedef union _AlgTstValue
{
  unsigned int	     i;
  unsigned long long l;
} AlgTstValue;

int             main(int argc, char **argv)
{
  int           n,
		option,
  		ok = 1,
		usage = 0,
		flgInv = 0,
  		flgLL = 0;
  AlgTstValue   v;
  const char	*optList = "hiln:";

  n = 16;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'i':
	flgInv = 1;
	break;
      case 'l':
	flgLL = 1;
	break;
      case 'n':
	if((sscanf(optarg, "%d", &n) != 1) || (n < 8) ||
	   ((flgLL == 0) && (n > 30)) || (n > 60))
	{
	  usage = 1;
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-i] [-h] [-l] [-n #]\n%s\n",
    *argv,
    "A test for AlgGreyCode(), AlgGreyCodeInv() and their long long\n"
    "equivalents.\n"
    "Options are:\n"
    "  i  Compute inverse.\n"
    "  l  Use unsigned long long functions instead of the unsigned int ones.\n"
    "  n  Number of bits.\n"
    "Values are read from the stadard input and written to the standard\n"
    "output.\n");
  }
  else
  {
    while(ok && (fscanf(stdin, "%Lu", &(v.l)) == 1))
    {
      if(flgInv)
      {
	if(flgLL)
	{
	  v.l = AlgGrayCodeInvLL(v.l);
	}
	else
	{
	  v.i = (unsigned int )(v.l);
	  v.l = AlgGrayCodeInv(v.i);
	}
      }
      else
      {
	if(flgLL)
	{
	  v.l = AlgGrayCodeLL(v.l);
	}
	else
	{
	  v.i = (unsigned int )(v.l);
	  v.l = AlgGrayCode(v.i);
	}
      }
      (void )fprintf(stdout, "0x%08llx %lld\n", v.l, v.l);
    }
  }
  return(!ok);
}
