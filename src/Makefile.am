bin_PROGRAMS= mat
mat_SOURCES = mat.cxx matrix.cxx matrix.h mat-combine.cxx mat-compare.cxx mat-diff.cxx mat-format.cxx mat-grep.cxx mat-nj.cxx
mat_CPPFLAGS = -Wall -Wextra  -std=c++14
mat_CXXFLAGS = -ggdb

format:
	clang-format -i *.cxx *.h
