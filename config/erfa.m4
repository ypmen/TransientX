# X_LIB_ERFA([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
AC_DEFUN([X_LIB_ERFA],
[
  AC_PROVIDE([X_LIB_ERFA])
  AC_MSG_CHECKING([ERFA installation])
  AC_TRY_LINK([#include "erfa.h"], [eraIcrs2g;],
                  have_erfa=yes, have_erfa=no)
  if test "$have_erfa" = "yes"; then
    AC_DEFINE([HAVE_ERFA],[1],[Define if the erfa library is present])
    LIBS="-lerfa $LIBS"
    [$1]
  else
    [$2]
  fi
])
