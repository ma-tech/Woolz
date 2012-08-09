#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyValueMixing_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzGreyValueMixing.c
* \author       Jianguo Rao
* \date         April 2003
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
* \brief	Produces a new object the grey values of which are a
* 		blend of the input objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzgreyvaluemixing "WlzGreyValueMixing"
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str);  

 
int main(int	argc,
	 char	**argv)
{

  FILE	       *inFile = NULL;   /* used to read Woolz object */ 
  int		option;

  double        xmidle         = 0;
  
  char                *inFileStr, *outFileStr, *targetfileStr;
  WlzErrorNum	       errNum = WLZ_ERR_NONE;
  WlzObject           *wObjS          = NULL; /* source Woolz object */
  WlzObject           *wObjT          = NULL; /* target Woolz object */
  WlzObject           *wObjW          = NULL; /* mixed grey value Woolz object
                                                 using source Obj as basis */

  /* read the argument list and check for an input file */

  static char	optList[] = "i:t:o:x:h";

   while( (option = getopt(argc, argv, optList)) != EOF )
   {
      switch( option )
      {
        case 'h':
             usage(argv[0]);
	     return(0);
        case 'i':
	    inFileStr = optarg;
	    break;
        case 't':
	    targetfileStr = optarg;
	    break;
        case 'o':
	    outFileStr = optarg;
	    break;
	case 'x':
	    if(sscanf(optarg, "%lg", &xmidle) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        default:
              return(0);
      }
   }

  /* read Woolz object */
  if((inFile = fopen(inFileStr, "r")) == NULL )
  {
    printf("cannot open the input woolz file.\n");
    exit(1);
  }

  if( !(wObjS = WlzReadObj(inFile, &errNum) ) )
  {
    printf("input Woolz Object Error.\n");
    fclose(inFile); 
    exit(1);
  }
  fclose(inFile); 
  inFile = NULL;

  /* read target  Woolz object */
  if((inFile = fopen( targetfileStr, "r")) == NULL )
  {
    printf("cannot open the input woolz file.\n");
    exit(1);
  }

  if( !(wObjT = WlzReadObj(inFile, &errNum) ) )
  {
    printf("input Woolz Object Error.\n");
    fclose(inFile); 
    exit(1);
  }
  fclose(inFile); 
  inFile = NULL;

  /*--- mixing grey values  ----*/
  if(errNum == WLZ_ERR_NONE)
  {
    wObjW = WlzGreyValueMixing_s(wObjS, wObjT, xmidle, &errNum);
    
    /* output the woolz object */
    if((inFile = fopen(outFileStr, "w")) == NULL )
    {
       printf("cannot open output file.\n");
       exit(1);
    } 
  
    if( (errNum =	WlzWriteObj(inFile, wObjW) ) !=WLZ_ERR_NONE )
    {
       printf("can not out put the transformed Woolz Object Error.\n");
       fclose(inFile);
       exit(1);
    }
    fclose(inFile);
    inFile = NULL;
    AlcFree(wObjW);
  }  
  
  AlcFree(wObjS);
  AlcFree(wObjT);
   return ( 0 );
}

static void usage(char *proc_str)
{
  (void )fprintf(stderr,
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tAutomatically produce a Woolz object from two Woolz objects. \n"
	  "\tthe output woolz object has the grey values  mixed from the\n"
	  "\tthe grey value of the two input Woolz objects. \n"
	  "Version: %s\n"
	  "\tOptions:\n"
	  "\t  -h        Help - prints this usage message\n"
          "\t  -i        input file name for the Woolz Object\n"
	  "\t  -t        input file name for the target Woolz Object\n"
          "\t  -o        output file name for the mixed grey value Woolz Object\n"
	  "\t  -x        input the mixing parameter x\n"
	  "\t            The grey value will be (1-x) Sobj + x Tobj \n"
	  "\t            where x is between 0 and 1\n",
	  proc_str,
	  WlzVersion());
  return;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
