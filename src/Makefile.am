bin_PROGRAMS= mat
mat_SOURCES = mat.cxx matrix.cxx matrix.h compare.cxx diff.cxx format.cxx grep.cxx nj.cxx mantel.cxx
mat_CPPFLAGS = -Wall -Wextra  -std=c++14
mat_CXXFLAGS = -ggdb

format:
	clang-format -i *.cxx *.h
