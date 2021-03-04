# BEAR_LIB_PGPLOT([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
AC_DEFUN([BEAR_LIB_PGPLOT],
[
  AC_PROVIDE([BEAR_LIB_PGPLOT])
  AC_MSG_CHECKING([PGPLOT installation])

  save_LIBS="$LIBS"
  LIBS="-lcpgplot -lpgplot $LIBS -lgfortran -lm -lquadmath -lX11"

  AC_TRY_LINK([#include "cpgplot.h"], [cpgopen("");cpgend();],
                  have_pgplot=yes, have_pgplot=no)
  if test "$have_pgplot" = "yes"; then
    AC_DEFINE([HAVE_PGPLOT],[1],[Define if the PGPLOT library is present])
    [$1]
  else
    LIBS="$save_LIBS"
    [$2]
  fi

  AM_CONDITIONAL(HAVE_PGPLOT,[test "$have_pgplot" = "yes"])

])
