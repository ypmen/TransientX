lib_LTLIBRARIES=libymw16.la
noinst_PROGRAMS=ymw16 ymw16_ne

AM_CPPFLAGS=-I$(top_srcdir)/include
LDFLAGS+=-L$(top_srcdir)/src/ymw16
LDADD=-lymw16
libymw16_la_SOURCES=dmdtau.c dora.c fermibubble.c frb_d.c galcen.c gum.c lmc.c localbubble.c ne_crd.c nps.c smc.c spiral.c thick.c thin.c ymw16par.c

ymw16_DEPENDENCIES=libymw16.la
ymw16_ne_DEPENDENCIES=libymw16.la

ymw16_SOURCES=ymw16.c
ymw16_ne_SOURCES=ymw16_ne.c

