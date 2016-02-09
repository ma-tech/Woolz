#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzLUTGreyTransformFromTxt_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzLUTGreyTransformFromTxt.c
* \author       Richard Baldock
* \date         May 2013
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
* \brief	Creates a grey LUT transform from text input
* 		to form a grey LUT transform object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzlutgreytransformfromtxt "WlzLUTGreyTransformFromTxt"
*/
 
/*!
\ingroup      BinWlz
\defgroup     wlzlutgreytransformfromtxt WlzLUTGreyTransformFromTxt
\par Name
WlzLUTGreyTransformFromTxt -
Creates a grey LUT transform from text input to form a grey LUT transform object.
Current version assumes index values from 0-255.
\par Synopsis
\verbatim
WlzLUTGreyTransformFromTxt [-f] [-o<output file>] [-h] [<LUT text file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints a usage message.</td>
  </tr>
  <tr>
    <td><b>-f</b></td>
    <td>Fills gaps with linear interpolation.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>output object file name.</td>
  </tr>
</table>

\par Description
Reads the LUT values from a text file or from standard input.
Assumes value-pairs (index, value)
and in this version a LUT will be created with
\f$0 \le \f$ index \f$\le  255\f$.
If the "-f" flag is used then values wil be
interpolated from 0-255 using the values that are set.
Values below the minimum index value and above
the maximum index  value are set to the min and max values respectively.

\par Examples

\par File
\ref WlzLUTGreyTransformFromTxt.c "WlzLUTGreyTransformFromTxt.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzfacts "WlzFacts(1)"
\ref WlzLUTGreyTransformNew "WlzLUTGreyTransformNew(3)"
\ref WlzLUTMergeToRGBA "WlzLUTMergeToRGBA(3)"
\ref wlzlutgreytransformmerge "WlzLUTTransformMerge(1)"
\ref WlzMakeEmpty "WlzMakeEmpty(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

typedef struct _LUTEntry{
  int	index;
  int	value;
} LUTEntry;

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

static int ComparLUTEntry(
  const void 	*p1,
  const void	*p2)
{
  LUTEntry	*entry1 = (LUTEntry *) p1;
  LUTEntry	*entry2 = (LUTEntry *) p2;

  if( entry1->index < entry2->index ){
    return( -1 );
  } else if ( entry1->index > entry2->index ){
    return( 1 );
  } else {
    return( 0 );
  }
}

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		usage = 0;
  int		fillFlg = 0;
  LUTEntry	LUTEntries[256];
  int		LUTEntryCount;
  int		entryIdx, index;
  int		val = 0, indexGap, indexDist;
  WlzObject	*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  const char	*errMsg;
  char 		*outFileStr;
  char		*inFileStr;
  static char	optList[] = "fho:",
		outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inFileStr = inFileStrDef;
  /* Parse the command line. */
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
    case 'f':
      fillFlg = 1;
      break;

    case 'o':
      outFileStr = optarg;
      break;

    case 'h':
    default:
      usage = 1;
      break;
    }
  }
  if(usage == 0)
  {
    if((outFileStr == NULL) || (*outFileStr == '\0'))
    {
      usage = 1;
    }
    if(usage == 0)
    {
      if((argc - optind) > 0)
      {
	inFileStr = *(argv+optind);
      }
    }
  }
  ok = (usage == 0);
  /* Read  input LUT from infile or stdin. */
  if(ok)
  {
    if( *inFileStr == '-' ){
      fP = stdin;
    } else {
      if( (fP = fopen(inFileStr, "r")) == NULL ){
	 ok = 0;
	 (void) fprintf(stderr,
	                 "%s: failed to open input file %s\n",
			 *argv, inFileStr);
      }
    }
    if(ok){
      /* read the LUT index and values */
      for( LUTEntryCount=0; LUTEntryCount < 256; LUTEntryCount++){
	if( fscanf(fP, "%d,%d", &LUTEntries[LUTEntryCount].index, &LUTEntries[LUTEntryCount].value) < 2 ){
	  break;
	}
      }
      /* check for no data */
      if( LUTEntryCount == 0 ) {
	ok = 0;
	(void) fprintf(stderr,
		       "%s: No LUT entries found please check the input file\n",
		       *argv);
      } else {
	/* sort the LUTEntries */
	(void) qsort(LUTEntries, (size_t) LUTEntryCount, sizeof(LUTEntry), ComparLUTEntry);

	/* check highest index add endpoint if needed */
	if( LUTEntries[LUTEntryCount-1].index < 255 ){
	  LUTEntries[LUTEntryCount].index = 255;
	  LUTEntries[LUTEntryCount].value = LUTEntries[LUTEntryCount - 1].value;
	  LUTEntryCount++;
	}
	
	/* check lowest index add startpoint if needed */
	if( LUTEntries[0].index > 0 ){
	  LUTEntries[LUTEntryCount].index = 0;
	  LUTEntries[LUTEntryCount].value = LUTEntries[0].value;
	  LUTEntryCount++;
	}

	/* re-sort the LUTEntries */
	(void) qsort(LUTEntries, (size_t) LUTEntryCount, sizeof(LUTEntry), ComparLUTEntry);

      }
    }
  }

  /* Create a LUT object for the range 0-255 */
  if(ok)
  {
    outObj = WlzAssignObject(
             WlzMakeLUTObject(WLZ_GREY_INT, 0, 255, &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to create new LUT object (%s)\n",
		     *argv, errMsg);
    }
  }
  /* Calculate the new LUT from the LUTEntries */
  if(ok){
    entryIdx = 0;
    for(index=0; index < 256; index++){
      if( index == LUTEntries[entryIdx].index ){
	val = LUTEntries[entryIdx].value;
	entryIdx++;
      }
      else if( index < LUTEntries[entryIdx].index ){
	indexGap = LUTEntries[entryIdx].index - LUTEntries[entryIdx-1].index;
	indexDist = index - LUTEntries[entryIdx-1].index;
	if( indexGap == 0 ){
	  /* shouldn't get here - something weird */
	  val =  LUTEntries[entryIdx].value;
	} else {
	  val = ((indexGap - indexDist) * LUTEntries[entryIdx-1].value + \
		 indexDist * LUTEntries[entryIdx].value) / indexGap;
	}
      }
      *((outObj->values.lut->val.inp) + index) = val;
    }
  }

  /* write object as required */
  if(ok)
  {
    if(((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL) ||
	(WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write LUT object to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
      fP = NULL;
    }
  }
  WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-f] [-o<output object>] [-h]\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "Reads the LUT values from a text file or from standard input. Assumes\n"
    "value-pairs (index, value) and in this version a LUT will be created\n"
    "with 0 <= index <= 255. If the \"-f\" flag is used then values wil be\n"
    "interpolated from 0-255 using the values that are set. Values below\n"
    "the minimum index value and above the maximum index value are set to\n"
    "the min and max values respectively.\n");
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
