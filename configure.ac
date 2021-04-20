##
# \file         configure.in
# \author       Bill Hill
# \date         September 2007
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2012],
# The University Court of the University of Edinburgh,
# Old College, Edinburgh, UK.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief        Process this file with autoconf to produce a configure script.
##

m4_define([WOOLZ_VERSION], [1.8.1])

AC_INIT(Woolz, WOOLZ_VERSION, ma-tech@igmm.ed.ac.uk)
AC_CANONICAL_TARGET

AC_GNU_SOURCE

AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_SRCDIR([libWlz/Wlz.h])

# Set defaults
AC_PREFIX_DEFAULT(/opt/MouseAtlas)
AC_DISABLE_SHARED
CP="cp"
RM="rm -f"
RMRF="rm -rf"
SIMD=''
CDEF_ARCH=""
USE_CTAGS=1
USE_ETAGS=0
HAVE_CJSON="no"
CJSON_DIR="no"
HAVE_JPEG="no"
JPEG_DIR="no"
HAVE_NIFTI="no"
NIFTI_DIR="no"
HAVE_TIFF="no"
TIFF_DIR="no"
HAVE_BZLIB="no"
HAVE_LZMALIB="no"
HAVE_ZLIB="no"
AC_SUBST(CP)
AC_SUBST(RM)
AC_SUBST(RMRF)

# Check CPU type
echo "Checking target cpu/system type"
case ${target_cpu} in
  i386|i486|i586|i686)
    CDEF_ARCH="-D__x86"
    ;;
  x86_64|amd64)
    CDEF_ARCH="-D__x86"
    ;;
  ia64)
    CDEF_ARCH="-D__x86"
    ;;
  alpha)
    CDEF_ARCH="-D__alpha"
    ;;
  sparc|sparc64)
    CDEF_ARCH="-D__sparc"
    ;;
  ppc|powerpc)
    CDEF_ARCH="-D__ppc"
    ;;
  mips)
    CDEF_ARCH="-D__mpis"
    ;;
  *)
    echo "Warning: CPU type ${target_cpu} unknown."
    ;;
esac

# Check for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL
LT_INIT
AC_PROG_MKDIR_P

# Check for system functions
AC_FUNC_MMAP

# Check for libraries.
AC_CHECK_LIB(m, pow)

# Check for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([	assert.h \
			ctype.h \
			errno.h \
			float.h \
			limits.h \
			stddef.h \
			stdlib.h \
			string.h \
			strings.h \
                  	sys/time.h \
                  	sys/types.h \
			time.h \
			unistd.h])

# Check for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM
AC_C_RESTRICT

# Check for library functions.
AC_FUNC_MKTIME
AC_FUNC_STAT
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([	floor \
			getcwd \
			gethostname \
			gettimeofday \
		        malloc \
			memset \
			mkdir \
			modf \
			pow \
			rand_r \
			realloc \
			regcmp \
			regcomp \
			sqrt \
			strcasecmp \
			strchr \
			strerror \
			strrchr \
			strsep \
                        strstr \
                        strtok \
                        strtok_r \
			timersub])


# Check for doxygen
AC_PATH_PROG(DOXYGEN, doxygen)
AC_CHECK_PROG(HAVE_DOXYGEN, doxygen, "yes", "no")
AM_CONDITIONAL(HAVE_DOXYGEN, test x$HAVE_DOXYGEN = xyes)

# Compilation and linking flags
CFLAGS=""
INCL_EXTFF=""
AC_SUBST(LIBS_EXTFF)
AC_SUBST(INCL_EXTFF)

AC_PROG_CC

AC_ARG_ENABLE(icc,
  [  --enable-icc            build using the Intel C compiler])
AS_IF([test x"$enable_icc" = "xyes"], [
  CFLAGS="${CFLAGS} -w2"
  CC="icc"
  CXX="icpc"
  LD="xild"
  AR="xiar"
], [
  CFLAGS="${CFLAGS} -pipe -Wall -Wno-deprecated"
  AS_IF([test -z "$CC" ], [CC="gcc"], [CC="$CC"])
  AS_IF([test -z "$CXX" ], [CXX="g++"], [CXX="$CXX"])
  AS_IF([test -z "$LD" ], [LD="ld"], [LD="$LD"])
])

CFLAGS="${CFLAGS} ${CDEF_ARCH}"

# Set target architecture flags
if test "$CC" = "gcc"
then
  case ${target_cpu} in
    i386|i486|i586)
      CFLAGS="${CFLAGS} -m32 -march=i586"
      ;;
    i686)
      CFLAGS="${CFLAGS} -m32 -march=i686"
      ;;
    x86_64|amd64)
      CFLAGS="${CFLAGS} -m64"
      ;;
    *)
      ;;
  esac
else
  #Intel CC
  if test "$CC" = "icc"
  then
    case ${target_cpu} in
      i386|i486|i586)
	CFLAGS="${CFLAGS} -m32 -xK"
	;;
      i686)
	CFLAGS="${CFLAGS} -xN"
	;;
      x86_64|amd64)
	CFLAGS="${CFLAGS} -m64"
	;;
      *)
	;;
    esac
  fi
fi

CPP="cpp"

AC_ARG_ENABLE(debug,
  [  --enable-debug	  build with debugging flags set])
AS_IF([test "x$enable_debug" = "xyes"], [
    CFLAGS="${CFLAGS} -g -fPIC -DDEBUG"
    LDFLAGS="${LDFLAGS} -g -DDEBUG"
])

AC_ARG_ENABLE(optimise,
  [  --enable-optimise	  build with optimisation flags set])
AS_IF([test "x$enable_optimise" = "xyes"], [
    # GNU CC
    if test "$CC" = "gcc"
    then
      case ${target_cpu} in
	i386|i486|i586|i686)
	  CFLAGS="${CFLAGS} -fPIC -O2 -fomit-frame-pointer"
	  ;;
	x86_64|amd64)
	  if test "x$SIMD" = "x"
	  then
	    SIMD="-msse3"
	  fi
          CFLAGS="${CFLAGS} -fPIC -O3"
	  ;;
        *)
          CFLAGS="${CFLAGS} -fPIC -O2"
	  ;;
      esac
    else
      #Intel CC
      if test "$CC" = "icc"
      then
	CFLAGS="${CFLAGS} -fPIC -O3 -ip"
      fi
    fi
    ])

AC_ARG_ENABLE(avx2,
  [AC_HELP_STRING([--enable-avx2],[enable AVX2 SIMD optimizations])],
  SIMD='-mavx2')

AC_ARG_ENABLE(openmp,
  [  --enable-openmp	  build with openmp flags set])
AS_IF([test "x$enable_openmp" = "xyes"], [
    # GNU CC
    if test "$CC" = "gcc"
    then
      CFLAGS="${CFLAGS} -fopenmp"
    fi
    # Intel CC
    if test "$CC" = "icc"
    then
      CFLAGS="${CFLAGS} -openmp"
    fi
])

AC_ARG_ENABLE(lto,
  [  --enable-lto	          build with link time optimisation flags set])
AS_IF([test "x$enable_lto" = "xyes"], [
    # GNU CC
    if test "$CC" = "gcc"
    then
      CFLAGS="${CFLAGS} -flto"
    fi
    # Intel CC
    if test "$CC" = "icc"
    then
      CFLAGS="${CFLAGS} -ipo"
    fi
])

AC_ARG_ENABLE(profile,
  [  --enable-profile        build with profiling set])
AS_IF([test "x$enable_profile" = "xyes"], [
    CFLAGS="${CFLAGS} -pg"
    LDFLAGS="${LDFLAGS} -pg"
])

AC_ARG_ENABLE(extff,
  [  --enable-extff	  build with external file format support])
AM_CONDITIONAL([BUILD_EXTFF], [test x"$enable_extff" = x"yes"])
AC_SUBST(ENABLE_EXTFF, ${enable_extff})

AC_ARG_ENABLE(test,
  [  --enable-test		  build test programs])
AM_CONDITIONAL([BUILD_TEST], [test x"$enable_test" = x"yes"])

if test ${USE_CTAGS} = 1
then
  AC_CHECK_PROG(GNC_CTAGS_FILE, ctags, tags)
fi

AC_ARG_ENABLE(bzip2,
              AS_HELP_STRING([--disable-bzip2], [build without bzip2 support]))

AC_ARG_ENABLE(lzma,
              AS_HELP_STRING([--disable-lzma], [build without lzma support]))

AC_ARG_ENABLE(zlib,
              AS_HELP_STRING([--disable-zlib], [build without zlib support]))

AC_CONFIG_FILES([	Makefile \
			libAlc/Makefile \
			libAlg/Makefile \
			libWlz/Makefile \
			binWlz/Makefile \
			binAlgTst/Makefile \
			binWlzTst/Makefile ])

if test x${enable_extff} = "xyes"
then
  AC_ARG_WITH(jpeg, [  --with-jpeg=DIR         Directory containing the jpeg
                          lib and include directories.],
	      JPEG_DIR=$withval, JPEG_DIR=no)
  if test ${JPEG_DIR} = "no"
  then
    AC_CHECK_HEADER(jpeglib.h,
		    AC_CHECK_LIB(jpeg, jpeg_CreateCompress,
				 [
				   HAVE_JPEG=yes
				   LIBS_EXTFF="-ljpeg ${LIBS_EXTFF}"
				 ]))
  else
    HAVE_JPEG=yes
    INCL_EXTFF="${INCL_EXTFF} -I${JPEG_DIR}/include"
    LIBS_EXTFF="-L${JPEG_DIR}/lib -ljpeg ${LIBS_EXTFF}"
  fi

  AC_ARG_WITH(tiff, [  --with-tiff=DIR         Directory containing the tiff
                          lib and include directories.],
	      TIFF_DIR=$withval, TIFF_DIR=no)
  if test ${TIFF_DIR} = "no"
  then
    AC_CHECK_HEADER(tiff.h,
		    AC_CHECK_LIB(tiff, TIFFOpen,
				 [
				   HAVE_TIFF=yes
				   LIBS_EXTFF="-ltiff  ${LIBS_EXTFF}"
				 ]))
  else
    HAVE_TIFF=yes
    INCL_EXTFF="${INCL_EXTFF} -I${TIFF_DIR}/include"
    LIBS_EXTFF="-L${TIFF_DIR}/lib -ltiff ${LIBS_EXTFF}"
  fi

  AC_ARG_WITH(nifti, [  --with-nifti=DIR         Directory containing NIfTI
                          lib and include directories.],
	      NIFTI_DIR=$withval, NIFTI_DIR=no)
  if test ${NIFTI_DIR} = "no"
  then
    AC_CHECK_HEADER(nifti1_io.h,
		    AC_CHECK_LIB(nifti, nifti_simple_init_nim,
				 [
				   HAVE_NIFTI=yes
				   LIBS_EXTFF=" -lniftiio -lznz ${LIBS_EXTFF}"
				 ]))
  else
    HAVE_NIFTI=yes
    INCL_EXTFF="${INCL_EXTFF} -I${NIFTI_DIR}/include -I${NIFTI_DIR}/include/nifti"
    LIBS_EXTFF=" -L${NIFTI_DIR}/lib -lniftiio -lznz ${LIBS_EXTFF}"
  fi
  AC_ARG_WITH(cjson, [  --with-cjson=DIR        Directory containing the cJSON
                          lib and include directories.],
	      CJSON_DIR=$withval, CJSON_DIR=no)
  if test ${CJSON_DIR} = "no"
  then
    AC_CHECK_HEADER(cJSON.h,
		    AC_CHECK_LIB(cjson, cJSON_Version,
				 [
				   HAVE_CJSON=yes
				   LIBS_EXTFF=" -lcjson ${LIBS_EXTFF}"
				 ]))
  else
    HAVE_CJSON=yes
    INCL_EXTFF="${INCL_EXTFF} -I${CJSON_DIR}/include"
    LIBS_EXTFF=" -L${CJSON_DIR}/lib -lcjson ${LIBS_EXTFF}"
  fi

  if test ${ENABLE_EXTFF} = "yes"
  then
    AS_IF([test x"$enable_bzip2" != x"no"], [
      AC_SEARCH_LIBS([BZ2_bzCompress], [bz2], [HAVE_BZLIB=yes])
    ])
    AS_IF([test x"$enable_zlib" != x"no"], [
      AC_SEARCH_LIBS([zlibVersion], [z], [HAVE_ZLIB=yes])
    ])
    if test ${HAVE_BZLIB} = "no"
    then
      CFLAGS="${CFLAGS} -DHAVE_BZLIB=0"
    else
      CFLAGS="${CFLAGS} -DHAVE_BZLIB=1"
    fi
    AS_IF([test x"$enable_lzma" != x"no"], [
      AC_SEARCH_LIBS([lzma_version_string], [lzma], [HAVE_LZMALIB=yes])
    ])
    if test ${HAVE_BZLIB} = "no"
    then
      CFLAGS="${CFLAGS} -DHAVE_BZLIB=0"
    else
      CFLAGS="${CFLAGS} -DHAVE_BZLIB=1"
    fi
    if test ${HAVE_LZMALIB} = "no"
    then
      CFLAGS="${CFLAGS} -DHAVE_LZMALIB=0"
    else
      CFLAGS="${CFLAGS} -DHAVE_LZMALIB=1"
    fi
    if test ${HAVE_ZLIB} = "no"
    then
      CFLAGS="${CFLAGS} -DHAVE_ZLIB=0"
    else
      CFLAGS="${CFLAGS} -DHAVE_ZLIB=1"
    fi
    if test ${NIFTI_DIR} = "no"
    then
      CFLAGS="${CFLAGS} -DHAVE_NIFTI=0"
    else
      CFLAGS="${CFLAGS} -DHAVE_NIFTI=1"
    fi
  fi

  AC_CONFIG_FILES([	libbibfile/Makefile \
			libhguDlpList/Makefile \
			libReconstruct/Makefile \
  			libWlzExtFF/Makefile \
  		        libWlzBnd/Makefile \
			binWlzApp/Makefile \
			binWlzExtFF/Makefile])
fi

case $host in
  *mingw*)
    LIBS="-lregex $LIBS"
    ;;
  *cygwin*)
    LIBS="-lregex $LIBS"
    ;;
esac

AS_IF([ test "$CC" = "gcc" -a "x$enable_optimise" = "xyes" ],
      CFLAGS="$CFLAGS $SIMD")

AC_OUTPUT

AC_MSG_RESULT([
Woolz core has been configured on ${host} with

  Source directory:                 ${srcdir}
  Installation directory:           ${prefix}
  Target CPU:                       ${target_cpu}
  C compiler:                       ${CC}
  C++ compiler:                     ${CXX}
  Linker:                           ${LD}
  Archiver:                         ${AR}
  CFLAGS:                          ${CFLAGS}
  LDFLAGS:                         ${LDFLAGS}
  LIBS needed for configure:       ${LIBS_EXTFF} ${LIBS}
  Woolz external file format code:  ${enable_extff}
  Test programs:                    ${enable_test}

To build Woolz:                     make
To install Woolz:                   make install
To build Woolz documentation:       make doc
To install Woolz documentation:     make install-doc
])

