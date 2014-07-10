# -- src/Makefile.am --


# public header files that needs to be installed in order for people to use the library
include_HEADERS = TG.hpp TGUtils.hpp

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# Use libtool to create a library in the lib directory
lib_LTLIBRARIES = libtplgene.la

## Core Library ##

# Specify the library sources
libtplgene_la_LIBADD = @BOOST_LIBS@
libtplgene_la_SOURCES =  SDATA.cpp TGUtils.cpp\
						 TGTypes.cpp


# Add extra libraries to this. Using -x c++ to force all files to be compile as c++
libtplgene_la_CPPFLAGS= $(BOOST_CPPFLAGS) -x c++

# Make this a shared library
libtplgene_la_LDFLAGS=-shared -lboost_program_options -version-info $(VERSION_INFO_CORE)


###############################################################################
