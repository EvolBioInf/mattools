AC_INIT([mattools], [1])
AM_INIT_AUTOMAKE([-Wall foreign ])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CPP
AC_PROG_CXX
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX([14], [], [mandatory])

AC_CHECK_HEADERS([err.h errno.h])

AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T
AC_CHECK_HEADER_STDBOOL

AC_CHECK_FUNCS([memmove sqrt])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 docs/mat.1
 docs/Makefile
])
AC_OUTPUT
