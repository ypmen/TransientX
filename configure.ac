#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([TransientX], [0.0], [ypmen@mpifr-bonn.mpg.de])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([include/AVL.h])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
CXXFLAGS='-std=c++11 -O3'

LT_INIT

# Checks for libraries.
AX_BOOST_BASE([1.56], [], [echo "Boost library was not found!";exit -1])
AC_CHECK_LIB([cfitsio], [ffgky], [], [echo "cfitsio library was not found!";exit -1])
AC_CHECK_LIB([png], [png_set_flush], [], [echo "library not found!";exit -1])
AC_CHECK_LIB([fftw3], [fftw_malloc], [], [echo "library not found!";exit -1])
AC_CHECK_LIB([fftw3_threads], [fftw_init_threads], [], [echo "library not found!";exit -1])
AC_CHECK_LIB([fftw3f], [fftwf_malloc], [], [echo "library not found!";exit -1])
AC_CHECK_LIB([fftw3f_threads], [fftwf_init_threads], [], [echo "library not found!";exit -1])

BEAR_LIB_PGPLOT([echo "library found"], [echo "library not found!";exit -1])
X_LIB_PLOTX([echo "library found"], [echo "library not found!";exit -1])
X_LIB_SOFA([echo "library found"], [echo "library not found!"])

AC_OPENMP
AX_CHECK_X86_FEATURES

CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T


# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_CHECK_FUNCS([floor memset pow sqrt strchr])

AC_CONFIG_FILES([Makefile
		 src/Makefile
                 src/container/Makefile
                 src/formats/Makefile
                 src/module/Makefile
                 src/singlepulse/Makefile
                 src/utils/Makefile
                 src/ymw16/Makefile])
AC_OUTPUT
