bin_PROGRAMS=dada_benchmarkdb dada_fildb dada_dbdedispdb dada_dbsearch dada_dbtim dada_dbnull

AM_CPPFLAGS=-I$(top_srcdir)/include
LDADD= -lxdada -lpsrdada -lxmodule -lxutils -lxcontainer -lxformats -lymw16

dada_benchmarkdb_SOURCES=dada_benchmarkdb.cpp
dada_fildb_SOURCES=dada_fildb.cpp
dada_dbdedispdb_SOURCES=dada_dbdedispdb.cpp
dada_dbsearch_SOURCES=dada_dbsearch.cpp ../singlepulse/boxcar.cpp ../singlepulse/cluster.cpp ../singlepulse/candplot.cpp
dada_dbtim_SOURCES=dada_dbtim.cpp
dada_dbnull_SOURCES=dada_dbnull.cpp
