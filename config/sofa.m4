# X_LIB_SOFA([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
AC_DEFUN([X_LIB_SOFA],
[
  AC_PROVIDE([X_LIB_SOFA])
  AC_MSG_CHECKING([SOFA installation])
  AC_TRY_LINK([#include "sofa.h"], [iauIcrs2g;],
                  have_sofa=yes, have_sofa=no)
  if test "$have_sofa" = "yes"; then
    AC_DEFINE([HAVE_SOFA],[1],[Define if the sofa library is present])
    LIBS="-lsofa_c $LIBS"
    [$1]
  else
    [$2]
  fi
])
