bin_PROGRAMS=transientx_fil replot_fil

AM_CPPFLAGS=-I$(top_srcdir)/include
LDFLAGS=-L$(top_srcdir)/src/container -L$(top_srcdir)/src/formats -L$(top_srcdir)/src/utils -L$(top_srcdir)/src/module -L$(top_srcdir)/src/ymw16
LDADD=-lmodule -lutils -lcontainer -lformats -lymw16

transientx_fil_SOURCES=dedisperse_fil.cpp singlepulse.cpp
replot_fil_SOURCES=replot_fil.cpp candidate.cpp
