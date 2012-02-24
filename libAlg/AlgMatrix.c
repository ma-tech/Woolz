#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrix_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrix.c
* \author       Bill Hill
* \date         October 2010
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
* \brief	Matrix allocation and maintenance functions.
* \ingroup	AlgMatrix
*/

#include <string.h>
#include <Alg.h>

static void 			AlgMatrixLLRInsert(
				  AlgMatrixLLR *mat,
				  size_t row,
				  AlgMatrixLLRE *p);
static AlgError			AlgMatrixRectWriteAscii(
				  AlgMatrixRect *mat,
				  FILE *fP);
static AlgError			AlgMatrixSymWriteAscii(
				  AlgMatrixSym *mat,
				  FILE *fP);
static AlgError			AlgMatrixLLRWriteAscii(
				  AlgMatrixLLR *mat,
				  FILE *fP);

/*!
* \return	New matrix or matrix with core NULL on error.
* \ingroup	AlgMatrix
* \brief	Allocates a new matrix or the requested type.
* \param	aType			Matrix type.
* \param	nR			Number of rows.
* \param	nC			Number of columns.
* \param	nE			Number of list entries to allocate,
* 					only used for linked list row matrices,
* 					may be zero.
* \param	tol			Matrix tollerance value, only used for
* 					linked list row matrices.
* \param	dstErr			Destination error pointer, may be NULL.
*/
AlgMatrix	AlgMatrixNew(AlgMatrixType aType, size_t nR, size_t nC,
			     size_t nE, double tol, AlgError *dstErr)
{
  AlgMatrix	mat;
  AlgError 	errNum = ALG_ERR_NONE;

  mat.core = NULL;
  switch(aType)
  {
    case ALG_MATRIX_RECT:
      mat.rect = AlgMatrixRectNew(nR, nC, &errNum);
      break;
    case ALG_MATRIX_SYM:
      mat.sym = AlgMatrixSymNew(nR, &errNum);
      break;
    case ALG_MATRIX_LLR:
      mat.llr = AlgMatrixLLRNew(nR, nC, nE, tol, &errNum);
      break;
    default:
      errNum = ALG_ERR_MATRIX_TYPE;
      break;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mat);
}

/*!
* \return	New rectangular matrix or NULL on error.
* \ingroup	AlgMatrix
* \brief	Allocates a new rectangular matrix.
* \param	nR			Number of rows.
* \param	nC			Number of columns.
* \param	dstErr			Destination error pointer, may be NULL.
*/
AlgMatrixRect	*AlgMatrixRectNew(size_t nR, size_t nC, AlgError *dstErr)
{
  AlgMatrixRect *mat;
  AlgError 	errNum = ALG_ERR_NONE;

  if((mat = (AlgMatrixRect *)AlcCalloc(1, sizeof(AlgMatrixRect))) == NULL)
  {
    errNum = ALG_ERR_MALLOC;
  }
  else
  {
    mat->type = ALG_MATRIX_RECT;
    mat->nR = nR;
    mat->nC = nC;
    if((errNum = AlcDouble2Calloc(&(mat->array), nR, nC)) != ALG_ERR_NONE)
    {
      AlcFree(mat);
      mat = NULL;
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mat);
}

/*!
* \ingroup	AlcMatrix
* \brief	Frees a matrix.
* \param	mat			Given matrix.
*/
void		AlgMatrixFree(AlgMatrix mat)
{
  if(mat.core)
  {
    switch(mat.core->type)
    {
      case ALG_MATRIX_RECT:
        AlgMatrixRectFree(mat.rect);
        break;
      case ALG_MATRIX_SYM:
        AlgMatrixSymFree(mat.sym);
        break;
      case ALG_MATRIX_LLR:
        AlgMatrixLLRFree(mat.llr);
        break;
      default:
        break;
    }
  }
}
/*!
* \ingroup	AlcMatrix
* \brief	Frees a rectangular matrix.
* \param	mat			Rectangular matrix.
*/
void		AlgMatrixRectFree(AlgMatrixRect *mat)
{
  if(mat)
  {
    AlcDouble2Free(mat->array);
    AlcFree(mat);
  }
}

/*!
* \return	New symmetric matrix or NULL on error.
* \ingroup	AlgMatrix
* \brief	Allocates a new symmetric matrix.
* \param	nN			Number of rows and columns.
* \param	dstErr			Destination error pointer, may be NULL.
*/
AlgMatrixSym	*AlgMatrixSymNew(size_t nN, AlgError *dstErr)
{
  size_t	nE;
  double	*p0;
  double	**p1;
  AlgMatrixSym  *mat;
  AlgError 	errNum = ALG_ERR_NONE;

  if((mat = (AlgMatrixSym *)AlcCalloc(1, sizeof(AlgMatrixSym))) == NULL)
  {
    errNum = ALG_ERR_MALLOC;
  }
  else
  {
    mat->type = ALG_MATRIX_SYM;
    mat->nR = nN;
    mat->nC = nN;
    nE = ((nN + 1) * nN) / 2;
    if(((p0 = (double *)AlcCalloc(nE * nE, sizeof(double))) == NULL) ||
       ((p1 = (double **)AlcMalloc(nN * sizeof(double *))) == NULL))
    {
      AlcFree(p0);
      errNum = ALG_ERR_MALLOC;
    }
    else
    {
      size_t idx;

      mat->array = p1;
      for(idx = 0; idx < nN; ++idx)
      {
	mat->array[idx] = p0;
	p0 += idx + 1;
      }
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mat);
}

/*!
* \ingroup	AlcMatrix
* \brief	Frees a symmetric matrix.
* \param	mat			symmetric matrix.
*/
void		AlgMatrixSymFree(AlgMatrixSym *mat)
{
  if(mat)
  {
    AlcDouble2Free(mat->array);
    AlcFree(mat);
  }
}

/*!
* \return	New linked list row matrix or NULL on error.
* \ingroup	AlgMatrix
* \brief	Allocates a new linked list row matrix.
* \param	nR			Number of rows.
* \param	nC			Number of columns.
* \param	nE			Number of list entries to allocate,
* 					may be zero.
* \param	tol			Matrix tollerance value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
AlgMatrixLLR	*AlgMatrixLLRNew(size_t nR, size_t nC, size_t nE, double tol,
				 AlgError *dstErr)
{
  AlgMatrixLLR	*mat;
  AlgError 	errNum = ALG_ERR_NONE;

  if(((mat = (AlgMatrixLLR *)AlcCalloc(1, sizeof(AlgMatrixLLR))) == NULL) ||
     ((mat->tbl = (AlgMatrixLLRE **)
                  AlcCalloc(nR, sizeof(AlgMatrixLLRE *)))== NULL))
  {
    errNum = ALG_ERR_MALLOC;
  }
  else if(nE > 0)
  {
    errNum = AlgMatrixLLRExpand(mat, nE);
  }
  if(errNum == ALG_ERR_NONE)
  {
    mat->type = ALG_MATRIX_LLR;
    mat->nR = nR;
    mat->nC = nC;
  }
  else
  {
    AlgMatrixLLRFree(mat);
    mat = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mat);
}

/*!
* \ingroup	AlcMatrix
* \brief	Frees a linked list row matrix.
* \param	mat		Linked list row matrix.
*/
void		AlgMatrixLLRFree(AlgMatrixLLR *mat)
{
  if(mat)
  {
    (void )AlcFreeStackFree(mat->blk);
    AlcFree(mat->tbl);
    AlcFree(mat);
  }
}

/*!
* \return	Linked list row matrix entry or NULL if none available.
* \ingroup	AlgMatrix
* \brief	Gets a linked list row matrix entry from the matrix. The
* 		return value may be NULL if none are available and the
* 		entries need to be expanded sing AlgMatrixLLRExpand().
* \note		Call this function rather than manipulating the free
* 		entries directly.
* \param	mat		Linked list row matrix.
*/
AlgMatrixLLRE	*AlgMatrixLLRENew(AlgMatrixLLR *mat)
{
  AlgMatrixLLRE *p;

  if((p = mat->freeStk) != NULL)
  {
    mat->freeStk = p->nxt;
    p->nxt = NULL;
    ++(mat->numEnt);
  }
  return(p);
}

/*!
* \ingroup	AlgMatrix
* \brief	Returns a linked list row matrix entry to the free
* 		stack of the matrix.
* \note		Call this function rather than manipulating the free
* 		entries directly.
* \param	mat		Linked list row matrix.
* \param	p		Linked list row matrix entry.
*/
void		AlgMatrixLLREFree(AlgMatrixLLR *mat, AlgMatrixLLRE *p)
{
  if(p != NULL)
  {
    p->nxt = mat->freeStk;
    mat->freeStk = p;
    --(mat->numEnt);
  }
}

/*!
* \return	Alg error code.
* \ingroup	AlgMatrix
* \brief	Writes a matrix in numeric ASCI format to the
* 		given file file. The rows are on separate lines and the
* 		columns of each row are white space seperated.
* \param	mat			Given matrix.
* \param	fP			Output file pointer.
*/
AlgError	AlgMatrixWriteAscii(AlgMatrix mat, FILE *fP)
{
  AlgError	errNum = ALG_ERR_NONE;

  if(mat.core != NULL)
  {
    switch(mat.core->type)
    {
      case ALG_MATRIX_RECT:
        errNum = AlgMatrixRectWriteAscii(mat.rect, fP);
	break;
      case ALG_MATRIX_SYM:
        errNum = AlgMatrixSymWriteAscii(mat.sym, fP);
	break;
      case ALG_MATRIX_LLR:
        errNum = AlgMatrixLLRWriteAscii(mat.llr, fP);
	break;
      default:
        errNum = ALG_ERR_MATRIX_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Alg error code.
* \ingroup	AlgMatrix
* \brief	Writes a rectangular matrix in numeric ASCI format to the
* 		given file file. The rows are on separate lines and the
* 		columns of each row are white space seperated.
* \param	mat			Given matrix.
* \param	fP			Output file pointer.
*/
static AlgError	AlgMatrixRectWriteAscii(AlgMatrixRect *mat, FILE *fP)
{
  int		idR,
  		idC;
  AlgError	errNum = ALG_ERR_NONE;

  for(idR = 0; idR < mat->nR; ++idR)
  {
    for(idC = 0; idC < mat->nC; ++idC)
    {
      (void )fprintf(fP, "%lg ", *(*(mat->array + idR) + idC));
    }
    if(fprintf(fP, "\n") != 1)
    {
      errNum = ALG_ERR_WRITE;
      break;
    }
  }
  return(errNum);
}

/*!
* \return	Alg error code.
* \ingroup	AlgMatrix
* \brief	Writes a symmetric matrix in numeric ASCI format to the
* 		given file file. The rows are on separate lines and the
* 		columns of each row are white space seperated.
* \param	mat			Given matrix.
* \param	fP			Output file pointer.
*/
static AlgError	AlgMatrixSymWriteAscii(AlgMatrixSym *mat, FILE *fP)
{
  int		idR,
  		idC;
  AlgError	errNum = ALG_ERR_NONE;

  for(idR = 0; idR < mat->nR; ++idR)
  {
    for(idC = 0; idC <= idR; ++idC)
    {
      (void )fprintf(fP, "%lg ", *(*(mat->array + idR) + idC));
    }
    for(idC = idR + 1; idC < mat->nC; ++idC)
    {
      (void )fprintf(fP, "%lg ", *(*(mat->array + idC) + idR));
    }
    if(fprintf(fP, "\n") != 1)
    {
      errNum = ALG_ERR_WRITE;
      break;
    }
  }
  return(errNum);
}

/*!
* \return	Alg error code.
* \ingroup	AlgMatrix
* \brief	Writes a linked list row matrix in numeric ASCI format to the
* 		given file file. The rows are on separate lines and the
* 		columns of each row are white space seperated.
* \param	mat			Given matrix.
* \param	fP			Output file pointer.
*/
static AlgError	AlgMatrixLLRWriteAscii(AlgMatrixLLR *mat, FILE *fP)
{
  int		idR,
  		idC;
  double	val;
  AlgMatrixLLRE *p;
  AlgError	errNum = ALG_ERR_NONE;

  for(idR = 0; idR < mat->nR; ++idR)
  {

    p = mat->tbl[idR];
    for(idC = 0; idC < mat->nC; ++idC)
    {
      if((p != NULL) && (p->col == idC))
      {
	val = p->val;
	p = p->nxt;
      }
      else
      {
        val = 0.0;
      }
      (void )fprintf(fP, "%lg ", val);
    }
    if(fprintf(fP, "\n") != 1)
    {
      errNum = ALG_ERR_WRITE;
      break;
    }
  }
  return(errNum);
}

/*!
* \return	New matrix read from the file. Matrix core member will be
* 		NULL on error.
* \ingroup	AlgMatrix
* \brief	Reads a matrix of the requested type from an ASCII file.
* 		If the type is ALG_MATRIX_SYM then only the first half
* 		of each row will be used, although all must be given.
* \param	mType			Required type of matrix.
* \param	tol			Tolerance value for ALG_MATRIX_LLR
* 					matrices (unused for other types).
* \param	fP			Input file pointer.
* \param	fSep			Field separator string containing
* 					possible field separator characters.
* \param	recMax			Maximum number of characters per
* 					input record.
* \param	dstErr			Destination error pointer, may be NULL.
*/
AlgMatrix	AlgMatrixReadAscii(AlgMatrixType mType, double tol,
				   FILE *fP, const char *fSep, size_t recMax,
				   AlgError *dstErr)
{
  size_t	nC,
  		nR,
		nV;
  AlgMatrix	mat;
  AlcVector	*vec = NULL;
  AlcErrno	alcErr = ALC_ER_NONE;
  AlgError	errNum = ALG_ERR_NONE;

  mat.core = NULL;
  if((vec = AlcVecReadDouble2Asci(fP, fSep, recMax, &nR, &nC,
                                  &alcErr)) == NULL)
  {
    switch(alcErr)
    {
      case ALC_ER_READ:
        errNum = ALG_ERR_READ;
	break;
      case ALC_ER_ALLOC:
        errNum = ALG_ERR_MALLOC;
	break;
      default:
        errNum = ALG_ERR_FUNC;
	break;
    }
  }
  else
  {
    switch(mType)
    {
      case ALG_MATRIX_RECT:
	if((mat.rect = (AlgMatrixRect *)
		       AlcCalloc(1, sizeof(AlgMatrixRect))) == NULL)
	{
	  errNum = ALG_ERR_MALLOC;
	}
	else
	{
	  mat.rect->type = ALG_MATRIX_RECT;
	  mat.rect->nR = nR;
	  mat.rect->nC = nC;
	  mat.rect->array = (double **)
		            AlcVectorToArray2D(vec, nR, nC, &alcErr);
	}
	break;
      case ALG_MATRIX_SYM:
	if(nR != nC)
	{
	  errNum = ALG_ERR_READ;
	}
	else if((mat.sym = AlgMatrixSymNew(nR, &errNum)) != NULL)
	{
	  size_t iR;

	  for(iR = 0; iR < nR; ++iR)
	  {
	    AlcVectorSetArray1D(vec, iR * nR, iR * (nR + 1),
	                        mat.sym->array[iR]);
	  }
	}
	break;
      case ALG_MATRIX_LLR:
	mat.llr = AlgMatrixLLRNew(nR, nC, 0, tol, &errNum);
	if(mat.llr != NULL)
	{
	  size_t iR,
	  	 iV;

	  iV = 0;
	  for(iR = 0; iR < nR; ++iR)
	  {
	    size_t iC;
	    AlgMatrixLLRE *p = NULL;

	    if((errNum = AlgMatrixLLRExpand(mat.llr, nC)) != ALG_ERR_NONE)
	    {
	      break;
	    }
	    for(iC = 0; iC < nC; ++iC)
	    {
	      double v;
	      
	      v = *(double *)AlcVectorItemGet(vec, iV);
	      if(fabs(v) > tol)
	      {
		if(p == NULL)
		{
		  p = mat.llr->tbl[iR] = AlgMatrixLLRENew(mat.llr);
		}
		else
		{
		  p->nxt = AlgMatrixLLRENew(mat.llr);
		  p = p->nxt;
		}
		p->col = iC;
		p->val = v;
	      }
	      ++iV;
	    }
	  }
	}
	break;
      default:
	errNum = ALG_ERR_MATRIX_TYPE;
	break;
    }
  }
  (void )AlcVectorFree(vec);
  if(errNum != ALG_ERR_NONE)
  {
    AlgMatrixFree(mat);
    mat.core = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mat);
}

/*!
* \return	Alg error code.
* \ingroup	AlgMatrix
* 		Errors can occur because of a memory allocation
* 		failure in ALG_MATRIX_LLR matrices or an invalid
* 		matrix type.
* \param	mat			Given matrix.
* \param	row			Row coordinate.
* \param	col			Column coordinate.
* \param	val			Matrix value.
*/
AlgError	AlgMatrixSet(AlgMatrix mat,
                             size_t row, size_t col, double val)
{
  AlgError 	errNum = ALG_ERR_NONE;

  if(mat.core)
  {
    switch(mat.core->type)
    {
      case ALG_MATRIX_LLR:
	errNum = AlgMatrixLLRSet(mat.llr, row, col, val);
	break;
      case ALG_MATRIX_SYM:
	if(col <= row)
	{
	  *(*(mat.sym->array + row) + col) = val;
	}
	else
	{
	  *(*(mat.sym->array + col) + row) = val;
	}
	break;
      case ALG_MATRIX_RECT:
	*(*(mat.sym->array + row) + col) = val;
	break;
      default:
	break;
    }
  }
  return(errNum);
}

AlgError	AlgMatrixLLRCopyInPlace(AlgMatrixLLR *aM, AlgMatrixLLR *bM)
{
  AlgError	errNum = ALG_ERR_NONE;                   

  AlgMatrixLLRZero(aM);
  errNum = AlgMatrixLLRExpand(aM, bM->numEnt);
  if(errNum == ALG_ERR_NONE)
  {
    size_t id0;

    for(id0 = 0; id0 < bM->nR; ++id0)
    {
      AlgMatrixLLRE *p,
      		    *q;

      p = NULL;
      q = bM->tbl[id0];
      while(q != NULL)
      {
	if(p == NULL)
	{
	  p = aM->tbl[id0] = AlgMatrixLLRENew(aM);
	}
	else
	{
	  p->nxt = AlgMatrixLLRENew(aM);
	  p = p->nxt;
	}
	p->col = q->col;
	p->val = q->val;
        q = q->nxt;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Alg error code.
* \ingroup	AlgMatrix
* \brief	This can have one of three actions:
* 		1. If no value exists at the given coordinates then
* 		the new value is added with a new entry.
* 		2. If a value exists at the given coordinate and the
* 		given value is non-zero then the entry value is
* 		replaced with the given value.
* 		3. If a value exists at the given coordinate and the
* 		given value is zero then the corresponding entry is
* 		removed.
* 		Errors can only occur because of a memory allocation
* 		failure. These can be avoided by preallocating
* 		sufficient entries using AlgMatrixLLRExpand().
* \param	mat			Linked list row matrix.
* \param	row			Row coordinate.
* \param	col			Column coordinate.
* \param	val			Matrix value.
*/
AlgError	AlgMatrixLLRSet(AlgMatrixLLR *mat,
                                size_t row, size_t col, double val)
{
  AlgMatrixLLRE *p;
  AlgError 	errNum = ALG_ERR_NONE;

  if(fabs(val) > mat->tol)
  {
    if(mat->freeStk == NULL)
    {
      errNum = AlgMatrixLLRExpand(mat, 0);
    }
    if(errNum == ALG_ERR_NONE)
    {
      p = AlgMatrixLLRENew(mat);
      p->col = col;
      p->val = val;
      AlgMatrixLLRInsert(mat, row, p);
    }
  }
  else
  {
    AlgMatrixLLRERemove(mat, row, col);
  }
  return(errNum);
}

/*!
* \return	Value in matrix at given coordinates.
* \ingroup	AlgMatrix
* \brief	Returns the value in the matrix at the given coordinates.
* \param	mat			Given matrix.
* \param	row			Given row.
* \param	col			Given column.
*/
double		AlgMatrixValue(AlgMatrix mat, size_t row, size_t col)
{
  double	val = 0.0;

  if(mat.core != NULL)
  {
    switch(mat.core->type)
    {
      case ALG_MATRIX_LLR:
        val = AlgMatrixLLRValue(mat.llr, row, col);
	break;
      case ALG_MATRIX_SYM:
	if(col <= row)
	{
	  val = *(*(mat.sym->array + row) + col);
	}
	else
	{
	  val = *(*(mat.sym->array + col) + row);
	}
	break;
      case ALG_MATRIX_RECT:
	val = *(*(mat.sym->array + row) + col);
	break;
      default:
	break;
    }
  }
  return(val);
}

/*!
* \return	Value in matrix at given coordinates.
* \ingroup	AlgMatrix
* \brief	Returns the value in the matrix at the given coordinates.
* \param	mat			Linked list row matrix.
* \param	row			Given row.
* \param	col			Given column.
*/
double		AlgMatrixLLRValue(AlgMatrixLLR *mat, size_t row, size_t col)
{
  double	val = 0.0;
  AlgMatrixLLRE *p;

  p = mat->tbl[row];
  if(p != NULL)
  {
    while((p->nxt != NULL) && (p->col < col))
    {
      p = p->nxt;
    }
    if(p->col == col)
    {
      val = p->val;
    }
  }
  return(val);
}

/*!
* \return	Alg error code.
* \ingroup	AlgMatrix
* \brief	Ensures that there are at least the requested number of free
* 		linked list row matrix entries available.
* \param	mat			Linked list row matrix.
* \param	nE			Requested minimum number of free entries
* 					or if zero a default number of entries
* 					are added.
*/
AlgError	AlgMatrixLLRExpand(AlgMatrixLLR *mat, size_t nE)
{
  size_t	idx;
  void		*b;
  AlgMatrixLLRE *p = NULL;
  AlgError 	errNum = ALG_ERR_NONE;
  const size_t	defNE = 4096;

  if(nE == 0)
  {
    nE = defNE;
  }
  else if(mat->maxEnt < mat->numEnt + nE)
  {
    nE = (mat->maxEnt + nE) - mat->numEnt;
  }
  else
  {
    nE = 0;
  }
  if(nE != 0)
  {
    if((p = (AlgMatrixLLRE *)AlcCalloc(nE, sizeof(AlgMatrixLLRE))) == NULL)
    {
      errNum = ALG_ERR_MALLOC;
    }
    else
    {
      if((b = AlcFreeStackPush(mat->blk, p, NULL)) == NULL)
      {
	errNum = ALG_ERR_MALLOC;
	AlcFree(p);
      }
      else
      {
	mat->blk = b;
	for(idx = 1; idx < nE; ++idx)
	{
	  p->nxt = mat->freeStk;
	  mat->freeStk = p;
	  ++p;
	}
	p->nxt = mat->freeStk;
	mat->freeStk = p;
	mat->maxEnt += nE;
      }
    }
  }
  return(errNum);
}

/*!
* \ingroup	AlgMatrix
* \brief	Inserts an entry at the given coordinates or changes
* 		it's value if it exists.
* \param	mat			Linked list row matrix.
* \param	row			Row coordinate.
* \param	p			Entry with column and value.
*/
static void 	AlgMatrixLLRInsert(AlgMatrixLLR *mat, size_t row,
				   AlgMatrixLLRE *g)
{
  AlgMatrixLLRE *q;

  q = mat->tbl[row];
  if(q == NULL)
  {
    mat->tbl[row] = g;
    g->nxt = NULL;
  }
  else if(g->col <= q->col)
  {
    if(g->col == q->col)
    {
      q->val = g->val;
      AlgMatrixLLREFree(mat, g);
    }
    else
    {
      mat->tbl[row] = g;
      g->nxt = q;
    }
  }
  else
  {
    AlgMatrixLLRE *p;

    while((q != NULL) && (q->col < g->col))
    {
      p = q;
      q = q->nxt;
    }
    if((q == NULL) || (q->col != g->col))
    {
      p->nxt = g;
      g->nxt = q;
    }
    else
    {
      q->val = g->val;
      AlgMatrixLLREFree(mat, g);
    }
  }
}

/*!
* \ingroup      AlgMatrix
* \brief        Sets all elements of the matrix to zero.
* \param        mat                     Given rectangular matrix.
*/
void            AlgMatrixRectZero(AlgMatrixRect *mat)
{
  (void )memset(mat->array[0], 0, mat->nR * mat->nC * sizeof(double));
}

/*!
* \ingroup      AlgMatrix
* \brief        Sets all elements of the matrix to zero.
* \param        mat                     Given symmetric matrix.
*/
void            AlgMatrixSymZero(AlgMatrixSym *mat)
{
  (void )memset(mat->array[0], 0,
                mat->nR * (mat->nR + 1) * sizeof(double) / 2);
}

/*!
* \ingroup      AlgMatrix
* \brief        Sets all elements of the matrix to zero.
* \param        mat                     Given linked list row matrix.
*/
void            AlgMatrixLLRZero(AlgMatrixLLR *mat)
{
  size_t        idR;

  for(idR = 0; idR < mat->nR; ++idR)
  {
    AlgMatrixLLRE *p,
                  *q;

    q = mat->tbl[idR];
    while(q != NULL)
    {
      p = q;
      q = q->nxt;

      p->nxt = mat->freeStk;
      mat->freeStk = p;
    }
    mat->tbl[idR] = NULL;
  }
  mat->numEnt = 0;
}

/*!
* \ingroup	AlgMatrix
* \brief	Removes the entry at the given coordinates.
* \param	mat			Linked list row matrix.
* \param	row			Row coordinate.
* \param	col			Column coordinate.
*/
void		AlgMatrixLLRERemove(AlgMatrixLLR *mat, size_t row, size_t col)
{
  AlgMatrixLLRE *p,
  		*q;

  if(mat->nR > 0)
  {
    q = mat->tbl[row];
    if(q != NULL)
    {
      p = NULL;
      while((q != NULL) && (q->col < col))
      {
        p = q;
	q = q->nxt;
      }
      if((q != NULL) && (q->col == col))
      {
        if(p == NULL)
	{
          mat->tbl[row] = q->nxt;
        }
	else
	{
	  p->nxt = q->nxt;
	}
	AlgMatrixLLREFree(mat, q);
      }
    }
  }
}
