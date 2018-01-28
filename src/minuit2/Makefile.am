## Process this file with automake to produce Makefile.in
#
# Author: Paul_Kunz@slac.stanford.edu
#

# The following is set, otherwise it would have to follow GNU conventions.
AUTOMAKE_OPTIONS = foreign

# --- get aclocal m4 files in config/m4
ACLOCAL_AMFLAGS = -I config/m4

EXTRA_DIST = \
  aclocal.m4      config.h.in    configure.ac  LICENSE      README.md \
  autom4te.cache  config.log     config        config.status  stamp-h1 \
  config.h        configure      libtool       Makefile.in

# --- The subdirectories.
SUBDIRS = \
  inc \
  src \
  test \
  doc

# --- The following is for platforms where automake dependencies don't work.
depend:
	for dir in $(SUBDIRS); \
	do ( test -d $$dir && cd $$dir && $(GMAKE) depend ) done

# --- Make the documentation in the build directory
docs:
	( cd doc; $(MAKE) docs )

# --- extra clean targets for Sun CC compiler and generated distribution file
clean-local:
	-rm -rf *.tar.gz
