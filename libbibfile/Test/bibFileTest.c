#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        bibFileTest.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Short and simple test main for the bibtex based file
*		syntax library.
*		Records are read from the file named on the command
*		line and then written to the standard output.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <bibFile.h>

int		main(int argc, char **argv)
{
  BibFileError	errFlag = BIBFILE_ER_NONE;
  BibFileRecord *record;
  char		*errMsg = NULL;
  FILE		*fP;

  if(argc != 2)
    (void )fprintf(stderr, "Usage: %s <file>\n", *argv);
  else
  {
    if((fP = fopen(*(argv + 1), "r")) != NULL)
    {
      while((errFlag == BIBFILE_ER_NONE) &&
	    ((errFlag = BibFileRecordRead(&record, &errMsg,
					  fP)) == BIBFILE_ER_NONE))
      {
	errFlag = BibFileRecordWrite(stdout, &errMsg, record);
	BibFileRecordFree(&record);
      }
      (void )fclose(fP);
      if((errFlag != BIBFILE_ER_NONE) && (errFlag != BIBFILE_ER_EOF))
      {
	(void )fprintf(stderr, "%s: Failed to read/write", *argv);
	if(errMsg)
	  (void )fprintf(stderr, ", %s.\n", errMsg);
	else
	  (void )fprintf(stderr, ".\n");
      }
    }
  }
  return(errFlag);
}
