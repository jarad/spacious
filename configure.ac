AC_PREREQ([2.59])
AC_LANG(C++)

AC_INIT([spacious], [1.0], [rjparker@ncsu.edu])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_AUX_DIR([config])
#AC_CONFIG_HEADERS([config.h:config-h.in])
#AC_CONFIG_SRCDIR([sic/eval.c])
AC_CONFIG_MACRO_DIR([m4])

#AM_INIT_AUTOMAKE([1.9 foreign])
ACLOCAL_AMFLAGS=-Im4 

# needed for R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
	echo "could not determine R_HOME"
	exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
LAPACK_LIBS=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
INC="$INC -I${R_HOME}/include"

AC_HEADER_STDC

AC_PROG_CC
AC_PROG_CXX

AC_CHECK_HEADERS([stdio.h math.h])

#AC_CHECK_FUNCS([bzero memset], [break])
#AC_CHECK_FUNCS([calloc strchr strrchr])
#AC_REPLACE_FUNCS([basename strspn strcspn strerror strsignal waitpid])

# use pthreads if available
AC_CHECK_LIB(pthread, pthread_create, [
	PTHREAD_LIBS="$PTHREAD_LIBS -lpthread"
	FEATURES="$FEATURES -DPTHREAD"
])

# cuda support
AC_ARG_WITH([cuda],
	AC_HELP_STRING([--with-cuda=DIR],
		[cuda root [/usr/local/cuda]]),
		[
			CUDA="$withval"
			if test x"$CUDA" = xyes; then
				CUDA="/usr/local/cuda"
			fi

			LDFLAGS="${LDFLAGS} -L${CUDA}/lib64 -lcublas -lcudart"
			INC="${INC} -I${CUDA}/include"
			FEATURES="$FEATURES -DCUDA"
			CUDA_OBJS="utils_cuda.o utils_cuda_kernels.o"
		])

AC_SUBST(CC)
AC_SUBST(CXX)
AC_SUBST(R_HOME)
AC_SUBST(CFLAGS)
AC_SUBST(CPICFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(INC)
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(PTHREAD_LIBS)
AC_SUBST(CUDA)
AC_SUBST(CUDA_OBJS)
AC_SUBST(FEATURES)

AC_OUTPUT
