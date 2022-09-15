noinst_LTLIBRARIES=libmodule.la

AM_CPPFLAGS=-I$(top_srcdir)/include
LDFLAGS=-L$(top_srcdir)/src/container -L$(top_srcdir)/src/formats -L$(top_srcdir)/src/utils
LDADD=-lcontainer -lformats -lutils

libmodule_la_SOURCES=patch.cpp preprocess.cpp preprocesslite.cpp downsample.cpp equalize.cpp baseline.cpp rfi.cpp subdedispersion.cpp boxcar.cpp cluster.cpp candplot.cpp archivewriter_tx.cpp psrfitsreader.cpp filterbankreader.cpp
