## Process this file with automake to produce Makefile.in .
#
# Author: Paul_Kunz@slac.stanford.edu
#

# The following is set; otherwise, it would follow GNU conventions.
# AUTOMAKE_OPTIONS = foreign

check_PROGRAMS = \
	test_Minuit2_Quad4FMain \
	test_Minuit2_Quad8FMain \
	test_Minuit2_Quad12FMain \
	test_Minuit2_Quad1FMain

test_Minuit2_Quad4FMain_SOURCES = \
	Quad4F.h \
	Quad4FMain.cxx

test_Minuit2_Quad1FMain_SOURCES = \
	Quad1F.h \
	Quad1FMain.cxx

test_Minuit2_Quad8FMain_SOURCES = \
	Quad8F.h \
	Quad8FMain.cxx

test_Minuit2_Quad12FMain_SOURCES = \
	Quad12F.h \
	Quad12FMain.cxx

AllSOURCES = \
	$(test_Minuit2_Quad4FMain_SOURCES) \
	$(test_Minuit2_Quad8FMain_SOURCES) \
	$(test_Minuit2_Quad12FMain_SOURCES) \
	$(test_Minuit2_Quad1FMain_SOURCES)

AM_CPPFLAGS = -I$(top_srcdir)/inc

##--- to build executable before installation on Linux
LDFLAGS_local = \
  -L$(abs_top_builddir)/src/.libs \
  -R$(abs_top_builddir)/src/.libs

test_Minuit2_Quad4FMain_LDADD = \
	$(top_builddir)/src/libMinuit2.la

test_Minuit2_Quad1FMain_LDADD = \
	$(top_builddir)/src/libMinuit2.la

test_Minuit2_Quad8FMain_LDADD = \
	$(top_builddir)/src/libMinuit2.la

test_Minuit2_Quad12FMain_LDADD = \
	$(top_builddir)/src/libMinuit2.la

test_Minuit2_Quad4FMain_LDFLAGS = $(LDFLAGS_local)

test_Minuit2_Quad1FMain_LDFLAGS = $(LDFLAGS_local)

test_Minuit2_Quad8FMain_LDFLAGS = $(LDFLAGS_local)

test_Minuit2_Quad12FMain_LDFLAGS = $(LDFLAGS_local)

# For platforms where automake dependency generation doesn't work.
depend: $(AllSOURCES)
	$(top_srcdir)/config/makedepend $(DEFS) $(INCLUDES) $(CPPFLAGS) $?

clean-local:
