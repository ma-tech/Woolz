

#ifndef WLZIO_H
#define WLZIO_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        bibFile.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Types and constants for the bibtex based file syntax
*		used for serial section data, ....
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <bibFile.h>
#include <Wlz.h>


/* From	WlzBibFileIOUtils.c */
extern WlzErrorNum write_Wlz3DSectionViewParams_Record(
  FILE			*fp,
  WlzThreeDViewStruct	*wlzViewStr);
extern WlzErrorNum parse_Wlz3DSectionViewParams_Record(
  BibFileRecord		*bibfileRecord,
  WlzThreeDViewStruct	*wlzViewStr);
extern WlzErrorNum write_WlzWarpTransformParams_Record(
  FILE			*fp,
  WlzBasisFnType	basisFnType,
  WlzMeshGenMethod	meshMthd,
  int			meshMinDst,
  int	 		meshMaxDst);
extern WlzErrorNum parse_WlzWarpTransformParams_Record(
  BibFileRecord		*bibfileRecord,
  WlzBasisFnType	*basisFnType,
  WlzMeshGenMethod	*meshMthd,
  int			*meshMinDst,
  int	 		*meshMaxDst);
extern WlzErrorNum write_WlzTiePointVtxs_Record(
  FILE		*fp,
  int		index,
  WlzDVertex3	dstVtx,
  WlzDVertex3	srcVtx);
extern WlzErrorNum parse_WlzTiePointVtxs_Record(
  BibFileRecord		*bibfileRecord,
  int		*index,
  WlzDVertex3	*dstVtx,
  WlzDVertex3	*srcVtx);

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* WLZIO_H */
