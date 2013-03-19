#ifndef ALGTYPE_H
#define ALGTYPE_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgType_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgType.h
* \author       Bill Hill
* \date         March 1999
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
* \brief        Type definitions for the Woolz numerical algorithm
*               library.
* \ingroup	Alg
*/


#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
extern "C" {
#endif
#endif /* WLZ_EXT_BIND */

/* Standard min, max, absolute value and nearest integer macros */
#define	ALG_MAX(X,Y)	(((X)>(Y))?(X):(Y))
#define	ALG_MIN(X,Y)	(((X)<(Y))?(X):(Y))
#define	ALG_MAXIDX(X,Y)	(((X)>(Y))?(0):(1))
#define	ALG_MINIDX(X,Y)	(((X)<(Y))?(0):(1))
#define	ALG_ABS(X)	(((X)>0)?(X):(-(X)))
#define	ALG_NINT(X)	((int)(((X)<0)?((X)-(0.5)):((X)+(0.5))))
#define	ALG_SQR(X)	((X)*(X))
#define	ALG_MAX3(X,Y,Z)	(((X)>(Y))?(((X)>(Z))?(X):(Z)):(((Y)>(Z))?(Y):(Z)))
#define	ALG_MIN3(X,Y,Z)	(((X)<(Y))?(((X)<(Z))?(X):(Z)):(((Y)<(Z))?(Y):(Z)))
#define ALG_MAXIDX3(X,Y,Z) \
			(((X)>(Y))?(((X)>(Z))?(0):(3)):(((Y)>(Z))?(1):(3)))
#define	ALG_MININD3(X,Y,Z) \
                        (((X)<(Y))?(((X)<(Z))?(0):(3)):(((Y)<(Z))?(1):(3)))

/* Standard math constants */
#define	ALG_M_E		(2.7182818284590452354)
#define	ALG_M_LOG2E	(1.4426950408889634074)
#define	ALG_M_LOG10E	(0.43429448190325182765)
#define	ALG_M_LN2	(0.69314718055994530942)
#define	ALG_M_LN10	(2.30258509299404568402)
#define ALG_M_PI	(3.14159265358979323846)
#define ALG_M_PI_2	(1.57079632679489661923)
#define ALG_M_PI_4	(0.78539816339744830961)
#define	ALG_M_1_PI	(0.31830988618379067154)
#define	ALG_M_2_PI	(0.63661977236758134308)
#define	ALG_M_2_SQRTPI	(1.12837916709551257390)
#define ALG_M_SQRT2	(1.41421356237309504880)
#define ALG_M_SQRT3	(1.73205080756887729353)
#define	ALG_M_SQRT1_2	(0.70710678118654752440)

/* A tollerance value for double precission arithmetic. */
#define ALG_DBL_TOLLERANCE	(1.0E-9)

/*!
* \enum		_AlgDistribution
* \brief	Statistical distributions.
* 		Typedef: ::AlgDistribution.
*/
typedef enum _AlgDistribution
{
  ALG_DISTRIBUTION_NORMAL,
  ALG_DISTRIBUTION_EXP,   
  ALG_DISTRIBUTION_POISSON,
  ALG_DISTRIBUTION_BINOMIAL 
} AlgDistribution;

/*!
* \enum		_AlgMatrixType
* \brief	Matrix representations.
*		Typedef: ::AlgMatrixType
*/
typedef enum _AlgMatrixType
{
  ALG_MATRIX_NULL = 0,		/*!< A NULL matric with no elements. */
  ALG_MATRIX_RECT,		/*!< Rectangular matrix, with storage
  				     for each element. These matrices should
				     be allocated using the libAlc array
				     allocation functions. */
  ALG_MATRIX_SYM,		/*!< Symmetric matrix, with storage
  				     for the upper triangle only. These
				     matrices should be allocated using the
				     libAlc symmetric array allocation
				     functions. */
  ALG_MATRIX_LLR		/*!< Sparse matrix stored in linked list
                                     row format. */
} AlgMatrixType;

/*!
* \struct	_AlgMatrix
* \brief	A union of all valid matrix types.
* 		Typedef: ::AlgMatrix..
*/
typedef union _AlgMatrix
{
  struct _AlgMatrixCore	*core;
  struct _AlgMatrixRect	*rect;
  struct _AlgMatrixSym	*sym;
  struct _AlgMatrixLLR	*llr;
} AlgMatrix;

/*!
* \struct	_AlgMatrixCore
* \brief	A core matrix type with members common to all matrix types.
* 		Typedef: ::AlgMatrixCore.
*/
typedef struct _AlgMatrixCore
{
  AlgMatrixType	type;		/*!< Matrix type. */
  size_t	nR;		/*!< Number of rows. */
  size_t	nC;		/*!< Number of columns. */
} AlgMatrixCore;

/*!
* \struct	_AlgMatrixRect
* \brief	Rectangular matrix.
* 		Typedef: ::AlgMatrixRect.
*/
typedef struct _AlgMatrixRect
{
  AlgMatrixType	type;		/*!< From AlgmatrixCore. */
  size_t	nR;		/*!< From AlgmatrixCore. */
  size_t	nC;		/*!< From AlgmatrixCore. */
  size_t	maxR;		/*!< Rows space allocated for. */
  size_t	maxC;		/*!< Columns space allocated for. */
  double	**array;	/* Array of elements. */
} AlgMatrixRect;

/*!
* \struct	_AlgMatrixSym
* \brief	Symmetric matrix.
* 		Typedef: ::AlgMatrixRect.
*/
typedef struct _AlgMatrixSym
{
  AlgMatrixType	type;		/*!< From AlgmatrixCore. */
  size_t	nR;		/*!< From AlgmatrixCore. */
  size_t	nC;		/*!< From AlgmatrixCore. */
  size_t	maxN;		/*!< Max rows/columns space allocated for. */
  double	**array;
} AlgMatrixSym;

/*!
* \struct	_AlgMatrixLLRE
* \brief	Entry in the linked list row matrix.
* 		Typedef: ::AlgMatrixLLRE.
*/
typedef struct _AlgMatrixLLRE
{
  size_t	col;		/*!< Column in matrix. */
  double	val;		/*!< Value in the row, column. */
  struct _AlgMatrixLLRE *nxt;   /*!< Next entry either in the value list or
                                     the free list. */
} AlgMatrixLLRE;

/*!
* \struct	_AlgMatrixLLRE
* \brief	Linked list row matrix, in which the values are stored
* 		in linked lists, with a linked list for each row of
* 		the matrix. This can be very efficient if the matrix
* 		is very sparse, but is very ineffiecient if the matrix
* 		is dense. The cross over is around 5 percent values
* 		being non-zero. At 1 percent of values being non-zero
* 		the linked list matrix is faster than a rectangular matrix
* 		(eg for matrix multiplication) by about a factor of ten.
*/
typedef struct _AlgMatrixLLR
{
  AlgMatrixType type;		/*!< From AlgmatrixCore. */
  size_t	nR;		/*!< From AlgmatrixCore. */
  size_t	nC;		/*!< From AlgmatrixCore. */
  size_t	numEnt;		/*!< Number of (no-zero) entries. */
  size_t	maxEnt;		/*!< Maximum number of entries. */
  double	tol;		/*!< Lowest absolute non-zero value. */
  void		*blk;	        /*!< Stack of blocks of triples allocated
  				     managed using AlcFreeStack. */
  AlgMatrixLLRE *freeStk;	/*!< Stack of free linked list entries. */
  AlgMatrixLLRE **tbl;		/*!< Table of matrix linked lists with a list
  				     for each row of the matrix. */
} AlgMatrixLLR;

typedef struct _AlgMatrixTriple
{
  size_t	row;		/*!< Row in matrix. */
  size_t	col;		/*!< Column in matrix. */
  double	val;		/*!< Value in the row, column. */
} AlgMatrixTriple;

/*!
* \enum		_AlgPadType
* \brief	Types of daat padding.
* 		Typedef: ::AlgPadType.
*/
typedef enum _AlgPadType
{
  ALG_PAD_NONE,              	/*!< No padding, same as padding with zeros. */
  ALG_PAD_ZERO,                 /*!< Pad data with zeros. */
  ALG_PAD_END                   /*!< Pad data with first/last data values. */
} AlgPadType;

/*!
* \struct	_ComplexD
* \brief	Complex number data type.
* 		Typedef: ::ComplexD.
*/
typedef struct _ComplexD
{
  double	re;
  double	im;
} ComplexD;


/*
* \enum		_AlgError
* \brief	Error codes.
* 		Typedef: ::AlgError.
*/
typedef enum _AlgError
{
  ALG_ERR_NONE		= (0),
  ALG_ERR_CONVERGENCE,		/*!< Failure to converge. */
  ALG_ERR_DIVZERO,   		/*!< Divide by zero. */
  ALG_ERR_FUNC,			/*!< Function parameters invalid. */
  ALG_ERR_MALLOC,		/*!< Memory allocation failure. */
  ALG_ERR_MATRIX_CONDITION,   	/*!< Matrix condition number out of range. */
  ALG_ERR_MATRIX_HOMOGENEOUS,	/*!< Homogeneous matrix. */
  ALG_ERR_MATRIX_SINGULAR,	/*!< Singular matrix. */
  ALG_ERR_MATRIX_TYPE,		/*!< Invalid matrix type given. */
  ALG_ERR_NONGLOBAL,   		/*!< Finds local solution, but fails to global
  				     solution. */
  ALG_ERR_READ,			/*!< Read failure. */
  ALG_ERR_WRITE,		/*!< Write failure. */
  ALG_ERR_MAX
} AlgError;

/*
* \enum		_AlgDbgMask
* \brief	Debug mask values.
* 		Typedef: ::AlgDbgMask.
*/
typedef enum _AlgDbgMask
{
  ALG_DBG_NONE          = (0),
  ALG_DBG_LVL_1         = (1),
  ALG_DBG_LVL_2         = (1<<1),
  ALG_DBG_LVL_3         = (1<<2),
  ALG_DBG_LVL_FN        = (1<<3)
} AlgDbgMask;

typedef AlgError        (*AlgDbgFn)(char *, ...);
 
extern AlgDbgFn		algDbgOutFn;
#define ALG_DBG_FN      (*algDbgOutFn)
#define ALG_DBG(F,M)    ((((F)&(algDbgMask))==(F))?ALG_DBG_FN M:ALG_ERR_NONE)
 
#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus 
}
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#endif /* ! ALGTYPE_H */
