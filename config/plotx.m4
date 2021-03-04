# X_LIB_PLOTX([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
AC_DEFUN([X_LIB_PLOTX],
[
  AC_PROVIDE([X_LIB_PLOTX])
  AC_MSG_CHECKING([PlotX installation])
  AC_LANG_PUSH([C++])
  AC_TRY_LINK([#include "plotx.h"], [PlotX::get_color;],
                  have_plotx=yes, have_plotx=no)
  AC_LANG_POP([C++])
  if test "$have_plotx" = "yes"; then
    AC_DEFINE([HAVE_PLOTX],[1],[Define if the plotx library is present])
    LIBS="-lplotx $LIBS"
    [$1]
  else
    [$2]
  fi
])
