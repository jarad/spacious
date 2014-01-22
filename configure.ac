
# pthreads support
#AC_ARG_WITH([pthreads],
#	AC_HELP_STRING([--with-pthreads],
#		[use pthreads for parallelization]),
#		[LD_PTHREADS="-lpthreads"],
#		[LD_PTHREADS=""])

# cuda support
#AC_ARG_WITH([cuda],
#	AC_HELP_STRING([--with-cuda=DIR],
#		[cuda root [/usr/local/cuda]]),
#		[
#			CUDA="$withval"
#			LD_CUDA="-L${CUDA}/lib64 -lcublas -lcudart"
#		],
#		[CUDA=""; LD_CUDA=""])

AC_PREREQ([2.59])
AC_LANG(C++)

AC_INIT([spacious], [1.0], [rjparker@ncsu.edu])
#AC_CONFIG_HEADERS([config.h:config-h.in])
#AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
#AC_CONFIG_SRCDIR([sic/eval.c])
AC_CONFIG_AUX_DIR([config])
#AC_CONFIG_MACRO_DIR([m4])

#AM_INIT_AUTOMAKE([1.9 foreign])

# needed for R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
	echo "could not determine R_HOME"
	exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
R_CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
R_CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
R_CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
R_CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
R_LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
LAPACK_LIBS=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`

#CUDA="oh_yeah"

AC_HEADER_STDC

AC_PROG_CC
AC_PROG_CXX

AC_CHECK_HEADERS([stdio.h math.h])

#SIC_WITH_READLINE
#SIC_VAR_SYS_ERRLIST
#SIC_VAR_SYS_SIGLIST

#AC_CHECK_FUNCS([bzero memset], [break])
#AC_CHECK_FUNCS([calloc strchr strrchr])
#AC_REPLACE_FUNCS([basename strspn strcspn strerror strsignal waitpid])
#http://www.gnu.org/software/autoconf-archive/ax_pthread.html
#AX_PTHREAD([PTHREADS=1], [PTHREADS=0])

#save_LIBS="$LIBS"
#test x"${with_readline-no}" != xno && LIBS="$LIBS -lreadline"
#AC_REPLACE_FUNCS([readline])
#AC_CHECK_FUNCS([add_history])
#LIBS="$save_LIBS"

AC_SUBST(CC)
AC_SUBST(CXX)
AC_SUBST(R_HOME)
AC_SUBST(R_CFLAGS)
AC_SUBST(R_CPICFLAGS)
AC_SUBST(R_CPPFLAGS)
AC_SUBST(R_CXXFLAGS)
AC_SUBST(R_LDFLAGS)
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(CUDA)

#,[test -z "$CONFIG_HEADERS" || echo timestamp > src/stamp-h.in])
#AC_CONFIG_FILES([src/Makefile])
AC_OUTPUT
