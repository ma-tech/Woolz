/***********************************************************************
* \Project:      Woolz
* \Title:        Wlz3DGreyValueMixing.c
* \Date:         23 April 2003 
* \Author:       J. Rao
* \Copyright:	 2003 Medical Research Council, UK.
*		 All rights reserved.
* \Address: 	 MRC Human Genetics Unit,
*		 Western General Hospital,
*		 Edinburgh, EH4 2XU, UK.
* \Purpose:      Woolz 3D Warping using MQ 
*		
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
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
  int		option,i;

  double        xmidle         = 0;
  
  char                *inFileStr, *outFileStr, *targetfileStr;
  const char	      *errMsg;
  WlzErrorNum	       errNum = WLZ_ERR_NONE;
  WlzObject           *wObjS          = NULL; /* source Woolz object */
  WlzObject           *wObjT          = NULL; /* target Woolz object */
  WlzObject           *wObjW          = NULL; /* mixed grey value Woolz object using source Obj as basis */
  char                *cstr;

  /* read the argument list and check for an input file */

  static char	optList[] = "i:t:o:x:h",

  opterr = 0;
 
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
  fprintf(stderr,
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tAutomatically produce a tetrahedron mesh from a woolz object. But\n"
	  "\tthe user should input some parameters to give a large cuboid which\n"
	  "\tcover the Woolz object. Sure, you can also just give a small cuboid\n"
	  "\tby transfer only part of your Woolz object. \n"
	  "\n"
	  "\tPlease report bugs to RAO.JIANGUO@hgu.mrc.ac.uk\n"
	  "\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
          "\t  -i        input file name for the Woolz Object\n"
          "\t  -t        input file name for the target Woolz Object\n"
          "\t  -o        output file name for the mixed grey value Woolz Object\n"
	  "\t  -x        input the mixing parameter x\n"
	  "\t            The grey value will be (1-x) Sobj + x Tobj \n"
	  "\t            where x is between 0 and 1\n"
	  "\t                                                      \n"
					  "",
	  proc_str);
  return;
}
