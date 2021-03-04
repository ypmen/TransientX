# X_LIB_YMW16([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
AC_DEFUN([X_LIB_YMW16],
[
  AC_PROVIDE([X_LIB_YMW16])
  AC_MSG_CHECKING([YMW16 installation])
  AC_TRY_LINK([#include "cn.h"], [dmdtau;],
                  have_ymw16=yes, have_ymw16=no)
  if test "$have_ymw16" = "yes"; then
    AC_DEFINE([HAVE_YMW16],[1],[Define if the ymw16 library is present])
    LIBS="-lymw16 $LIBS"
    [$1]
  else
    [$2]
  fi
])
