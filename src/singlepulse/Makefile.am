bin_PROGRAMS=transientx transientx_fil

AM_CPPFLAGS=-I$(top_srcdir)/include
LDFLAGS=-L$(top_srcdir)/src/container -L$(top_srcdir)/src/formats -L$(top_srcdir)/src/utils -L$(top_srcdir)/src/module -L$(top_srcdir)/src/ymw16
LDADD=-lboost_program_options -lmodule -lutils -lcontainer -lformats -lymw16

transientx_SOURCES=dedisperse.cpp singlepulse.cpp
transientx_fil_SOURCES=dedisperse_fil.cpp singlepulse.cpp
