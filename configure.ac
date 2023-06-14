#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([NCPA Propagation Modeling Suite], [2.2.0], [claus@olemiss.edu], ncpaprop)
#AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/common/include/parameterset.h])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG([C++])

# set default variables
default_petsc_build_dir="${abs_builddir}/extern/petsc"
default_petsc_branch="release"
default_petsc_version=""
petsc_url="https://gitlab.com/petsc/petsc.git"
default_slepc_build_dir="${abs_builddir}/extern/slepc"
default_slepc_version="3.18.0"
slepc_url_base="https://slepc.upv.es/download/distrib"
warning_compiler_flags="-Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas"

# Need OS name for automatic default architecture name
os_type=$(uname -s | awk '{print tolower($0)}')
default_petsc_arch_name_real="arch-${os_type}-c-real"
default_petsc_arch_name_complex="arch-${os_type}-c-complex"

# Set up options for automatic prerequisite building
AC_ARG_WITH([localpetsc],
	AS_HELP_STRING([--with-localpetsc],[Attempt to build a PETSc instance local to the $PACKAGE_NAME installation])
)
AC_ARG_WITH([localslepc],
	AS_HELP_STRING([--with-localslepc],[Attempt to build a SLEPc instance local to the $PACKAGE_NAME installation.  Automatically set to 'yes' if --with-localpetsc is specified])
)
AS_IF([test "x$with_localpetsc" = "xyes"],
	[with_localslepc="yes"] )

AC_ARG_WITH([localpetscdebug],
	AS_HELP_STRING([--with-localpetscdebug],[Build the debug version of PETSc/SLEPc instead of the default optimized version]),
	[],
	[with_localpetscdebug=no]
)
AC_ARG_WITH([insecuressl],
	AS_HELP_STRING([--with-insecuressl],[Ignore SSL certificate errors when downloading SLEPc]),
	[],
	[with_insecuressl=no]
)
AC_ARG_ENABLE([localpetscmpi],
	AS_HELP_STRING([--enable-localpetscmpi],[Build PETSc/SLEPc without MPI support]),
	[],
	[enable_localpetscmpi=no]
)
AC_ARG_ENABLE([autodependencies],
	AS_HELP_STRING([--enable-autodependencies],[Attempt to install dependencies using the default package manager (Linux only).  Currently supported package managers: apt, yum])
)
AC_ARG_ENABLE([compilerwarnings],
	AS_HELP_STRING([--enable-compilerwarnings],[Don't suppress additional compiler warnings (development only)])
)


# Environmental variables
AC_ARG_VAR(local_petsc_dir,[Override default local PETSc installation directory])
AC_ARG_VAR(local_slepc_dir,[Override default local SLEPc installation directory])
AC_ARG_VAR(local_petsc_arch_real,[Override automatic PETSc real architecture name])
AC_ARG_VAR(local_petsc_arch_complex,[Override automatic PETSc complex architecture name])
AC_ARG_VAR(local_petsc_version,Override automatic PETSc version number)
AC_ARG_VAR(local_slepc_version,Override automatic (or outdated) SLEPc version number)


AC_CANONICAL_HOST
AC_MSG_CHECKING([operating system])
AC_MSG_RESULT([$host_os])
# check host OS for support
case $host_os in
	darwin* )
		AS_IF([test "x$enable_autodependencies" = "xyes"],[
			AC_MSG_ERROR([Option --enable-autodependencies is not supported on MacOS.])
		])
		;;
	linux* )
		# For Linux, need to know if apt-get or yum is installed
		AC_CHECK_PROGS([default_package_manager],[apt-get yum],[])
		AC_MSG_NOTICE([Detected package manager = ${default_package_manager}])
		
		AS_IF([test "x${default_package_manager}" = "x"],[
			AS_IF([test "x$enable_autodependencies" = "xyes"],[
				AC_MSG_ERROR([No supported package manager detected.  Option --enable-autodependencies supports apt-get and yum.])
			])
		])
		AS_IF([test "x${enable_compilerwarnings}" != "xyes"],[
			warning_compiler_flags="${warning_compiler_flags} -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-result"
		])
		
		;;
	*)
		AC_MSG_ERROR([Your platform is not currently supported])
		;;
esac


# Look for bash and perl.  If not found, check with the user, cause
# that's unusual for a Linux-like system
AC_PATH_PROG([BASHPATH],[bash],)
AS_IF([test "x$BASHPATH" = "x"],[
	AC_MSG_FAILURE(["No bash interpreter found on your path!  Please check this."])
])
AC_PATH_PROG([PERLPATH],[perl],)
AS_IF([test "x$PERLPATH" = "x"],[
	AC_MSG_FAILURE(["No perl interpreter found on your path!  Please check this."])
])

AC_PROG_RANLIB
AC_PROG_CXX
AS_IF([test "x$CXX" = "x"],[
	AS_IF([test "x$enable_autodependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				AC_MSG_ERROR([Auto dependencies not available for MacOS.  Install g++ using your package manager of choice.])
			],
			[linux*],[
				AS_CASE([$default_package_manager],
					[*apt-get],[
						sudo $default_package_manager -y install g++
					],
					[*yum],[
						sudo $default_package_manager -y install gcc-c++
					],[
						AC_MSG_ERROR([No supported package manager found (apt and yum currently supported)])
					]
				)
		])
		unset ac_cv_prog_CXX
		unset ac_cv_cxx_compiler_gnu
		unset ac_cv_prog_cxx_g
		AC_PROG_CXX
		AS_IF([test "x${CXX}" = "x"],[
			AC_MSG_ERROR("Unable to install C++ compiler.  Please install manually and retry.")
		])
	],[
		AC_MSG_ERROR([Failed to detect a C++ compiler.  Please install the appropriate package or compile from source.])
	])
])
AC_PROG_CC
AS_IF([test "x$CC" = "x"],[
	AS_IF([test "x$enable_autodependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				AC_MSG_ERROR([Auto dependencies not available for MacOS.  Install gcc using your package manager of choice.])
			],
			[linux*],[
				AS_CASE([$default_package_manager],
					[*apt-get],[
						sudo $default_package_manager -y install gcc
					],
					[*yum],[
						sudo $default_package_manager -y install gcc
					],[
						AC_MSG_ERROR([No supported package manager found (apt and yum currently supported)])
					]
				)
		])
		unset ac_cv_prog_CC
		unset ac_cv_cc_compiler_gnu
		unset ac_cv_prog_cc_g
		AC_PROG_CC
		AS_IF([test "x${CC}" = "x"],[
			AC_MSG_ERROR("Unable to install C compiler.  Please install manually and retry.")
		])
	],[
		AC_MSG_ERROR([Failed to detect a C compiler.  Please install the appropriate package or compile from source.])
	])
])

AC_PROG_CPP

# Checks for libraries.
AC_CHECK_LIB([m],[cos],,[AC_MSG_ERROR([Missing M library])])

AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_autodependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				AC_MSG_ERROR([Auto dependencies not available for MacOS.  Install FFTW using your package manager of choice.])
			],
			[linux*],[
				AS_CASE([$default_package_manager],
					[*apt-get],[
						sudo $default_package_manager -y install libfftw3-dev
					],
					[*yum],[
						sudo $default_package_manager -y install fftw-devel
					],[
						AC_MSG_ERROR([No supported package manager found (apt and yum currently supported)])
					]
				)
		])
		unset ac_cv_lib_fftw3_fftw_plan_dft_1d
		AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],,[AC_MSG_ERROR([Could not install FFTW using ${default_package_manager}, please install manually])])
	],[
		AC_MSG_ERROR([Failed to detect FFTW.  Please install the appropriate fftw3 dev package or compile from source.])
	])
])

AC_CHECK_LIB([gslcblas], [cblas_dgemm],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_autodependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				AC_MSG_ERROR([Auto dependencies not available for MacOS.  Install GSL using your package manager of choice.])
			],
			[linux*],[
				AS_CASE([$default_package_manager],
					[*apt-get],[
						sudo $default_package_manager -y install libgsl-dev
					],
					[*yum],[
						sudo $default_package_manager -y install gsl-devel
					],[
						AC_MSG_ERROR([No supported package manager found (apt and yum currently supported)])
					]
				)
		])
		unset ac_cv_lib_gslcblas_cblas_dgemm
		AC_CHECK_LIB([gslcblas], [cblas_dgemm],,[AC_MSG_ERROR([Could not install GSL using ${default_package_manager}, please install manually])])
	],[
		AC_MSG_ERROR([Failed to detect GSL.  Please install the appropriate gsl dev package or compile from source.])
	])
])

AC_CHECK_LIB([gsl], [gsl_blas_dgemm],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				AC_MSG_ERROR([Auto dependencies not available for MacOS.  Install GSL using your package manager of choice.])
			],
			[linux*],[
				AS_CASE([$default_package_manager],
					[*apt-get],[
						sudo $default_package_manager -y install libgsl-dev
					],
					[*yum],[
						sudo $default_package_manager install gsl-devel
					],[
						AC_MSG_ERROR([No supported package manager found (apt and yum currently supported)])
					]
				)
		])
		unset ac_cv_lib_gsl_gsl_blas_dgemm
		AC_CHECK_LIB([gsl], [gsl_blas_dgemm],,[AC_MSG_ERROR([Could not install GSL using ${default_package_manager}, please install manually])])
	],[
		AC_MSG_ERROR([Failed to detect GSL.  Please install the appropriate gsl dev package or compile from source.])
	])
])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset mkdir modf pow sqrt strpbrk])

# Try to install PETSc and/or SLEPc if requested
AC_MSG_CHECKING([automatic local PETSc installation])
AS_IF(
	[test "x${with_localpetsc}" = "xyes"],
	[
		AC_MSG_RESULT([yes])
		AC_LANG_PUSH([C])

		# See if default PETSc directory is overridden
		AS_IF([test "x${local_petsc_dir}" != "x"],
			[petsc_build_dir=${local_petsc_dir}],
			[petsc_build_dir="$(pwd)/extern/petsc"])

		# See if default PETSc version is overridden
		AS_IF([test "x$local_petsc_version" = "x"],
			[
				AC_MSG_NOTICE([Using default PETSc])
				petsc_version=""
			],[
				AC_MSG_NOTICE([Using PETSc version ${local_petsc_version}])
				petsc_version="v${local_petsc_version}"
			]
		)

		AC_MSG_NOTICE([petsc_build_dir = ${petsc_build_dir}])

		# check for extra configuration options.
		petsc_extra_config="--with-fc=0"
	
		# See if default PETSc real arch name is overridden
		AS_IF([test "x${local_petsc_arch_real}" = "x"],
			[
				AS_IF([test "x${with_localpetscdebug}" = "xyes"],
					[
						petsc_local_arch_name_real="${default_petsc_arch_name_real}-debug"
					],
					[
						petsc_local_arch_name_real="${default_petsc_arch_name_real}"
					]
				)
			],[petsc_local_arch_name_real=$local_petsc_arch_real])
		
		AC_MSG_NOTICE([petsc_local_arch_name_real = ${petsc_local_arch_name_real}])
		
		# See if default PETSc complex arch name is overridden
		AS_IF([test "x${local_petsc_arch_complex}" = "x"],
			[
				AS_IF([test "x${with_localpetscdebug}" = "xyes"],
					[
						petsc_local_arch_name_complex="${default_petsc_arch_name_complex}-debug"
					],
					[
						petsc_local_arch_name_complex="${default_petsc_arch_name_complex}"
					]
				)
			],[petsc_local_arch_name_complex=$local_petsc_arch_complex])

		AC_MSG_NOTICE([petsc_local_arch_name_complex = ${petsc_local_arch_name_complex}])

		# debug version?
		AS_IF([test "x${with_localpetscdebug}" = "xyes"],
			[petsc_extra_config="${petsc_extra_config} --with-debugging=1"],
			[petsc_extra_config="${petsc_extra_config} --with-debugging=0"]
		)

		# Use MPI?
		AS_IF([test "x${enable_localpetscmpi}" = "xyes"],[
			
			AC_CHECK_PROGS(MPI_FOUND, [mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++])
			AS_IF([test "x$MPI_FOUND" = "x"],[
				petsc_extra_config="${petsc_extra_config} --download-openmpi"
			])
		],[
			petsc_extra_config="${petsc_extra_config} --with-mpi=0"
		])
	
		# Go get F2CBLAS-LAPACK?
		AX_BLAS([
			AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.])
		],[
			petsc_extra_config="${petsc_extra_config} --download-f2cblaslapack=1"
		])

		dnl AS_IF([test "x${with_localpetscfblas}" = "xyes"],
		dnl 	[
		dnl 		petsc_extra_config="${petsc_extra_config} --download-f2cblaslapack=1"
		dnl 	],[
		dnl 		petsc_extra_config="${petsc_extra_config} --with-fc=0"
		dnl 	])

		AC_MSG_NOTICE([petsc_extra_config = ${petsc_extra_config}])
		
		# make the install dir
		AC_MSG_CHECKING([PETSc build directory])
		AS_IF([test -d "$petsc_build_dir"],[
			AC_MSG_RESULT([yes])
		],[
			AC_MSG_RESULT([no, creating])
			mkdir -p ${petsc_build_dir}
			AS_IF([test -d "$petsc_build_dir"],,AC_MSG_ERROR(["Error creating PETSc build directory ${petsc_build_dir}"]))
		])

		# see if it's already been downloaded and/or built
		basedir=$(pwd)
		cd ${petsc_build_dir}
		AC_MSG_CHECKING([PETSc download])
		AS_IF([test -f "configure"],
			[AC_MSG_RESULT([yes])],
			[
				AC_MSG_RESULT([no, cloning])
				git clone -b ${default_petsc_branch} ${petsc_url} .
			]
		)

		# Version
		AS_IF([test "x${petsc_version}" != "x"],
			[
				AC_MSG_NOTICE([Checking out branch ${petsc_version}])
				git checkout ${petsc_version}
				AS_IF([test $? != 0],
					[AC_MSG_ERROR([Error switching to PETSc version ${petsc_version}])]
				)
			]
		)

		AC_MSG_CHECKING([for existing ${petsc_local_arch_name_real} architecture])
		AC_CHECK_FILE( ${petsc_build_dir}/${petsc_local_arch_name_real}/lib/petsc/conf/petscvariables,
			[AC_MSG_RESULT([yes])],
			[
				AC_MSG_RESULT([no, building])

				# build the real version with indicated options
				./configure PETSC_ARCH=${petsc_local_arch_name_real} PETSC_DIR=${petsc_build_dir} ${petsc_extra_config}
				
				AS_IF([test $? != 0],
					[AC_MSG_ERROR([Error running PETSc configure command: ./configure PETSC_ARCH=${petsc_local_arch_name_real} PETSC_DIR=${petsc_build_dir} ${petsc_extra_config}])]
				)
				
				make PETSC_ARCH=${petsc_local_arch_name_real} PETSC_DIR=${petsc_build_dir} all
				AS_IF([test $? != 0],
					[AC_MSG_ERROR([Error running make command for real PETSc: make PETSC_ARCH=${petsc_local_arch_name_real} PETSC_DIR=${petsc_build_dir} all])]
				)

				# see if it made the include file
				AC_MSG_CHECKING([successful ${petsc_local_arch_name_real} PETSc build])

				# clear the cached variable
				unset `$as_echo "ac_cv_file_${petsc_build_dir}/${petsc_local_arch_name_real}/lib/petsc/conf/petscvariables" | $as_tr_sh`

				AC_CHECK_FILE( ${petsc_build_dir}/${petsc_local_arch_name_real}/lib/petsc/conf/petscvariables,,
					[
						AC_MSG_ERROR([Error building real PETSc, ${petsc_build_dir}/${petsc_local_arch_name_real}/lib/petsc/conf/petscvariables not found])
					]
				)
			]
		)
		
		AC_MSG_CHECKING([for existing ${petsc_local_arch_name_complex} architecture])
		AC_CHECK_FILE( ${petsc_build_dir}/${petsc_local_arch_name_complex}/lib/petsc/conf/petscvariables,
			[AC_MSG_RESULT([yes])],
			[
				AC_MSG_RESULT([no, building])

				# build the real version with indicated options
				./configure PETSC_ARCH=${petsc_local_arch_name_complex} PETSC_DIR=${petsc_build_dir} ${petsc_extra_config} --with-scalar-type=complex
				
				AS_IF([test $? != 0],
					[AC_MSG_ERROR([Error running PETSc configure command: ./configure PETSC_ARCH=${petsc_local_arch_name_complex} PETSC_DIR=${petsc_build_dir} ${petsc_extra_config} --with-scalar-type=complex])]
				)
				
				make PETSC_ARCH=${petsc_local_arch_name_complex} PETSC_DIR=${petsc_build_dir} all
				AS_IF([test $? != 0],
					[AC_MSG_ERROR([Error running make command for complex PETSc: make PETSC_ARCH=${petsc_local_arch_name_complex} PETSC_DIR=${petsc_build_dir} all])]
				)

				# see if it made the include file
				AC_MSG_CHECKING([successful ${petsc_local_arch_name_complex} PETSc build])

				# clear the cached variable
				unset `$as_echo "ac_cv_file_${petsc_build_dir}/${petsc_local_arch_name_complex}/lib/petsc/conf/petscvariables" | $as_tr_sh`

				AC_CHECK_FILE( ${petsc_build_dir}/${petsc_local_arch_name_complex}/lib/petsc/conf/petscvariables,,
					[
						AC_MSG_ERROR([Error building complex PETSc, ${petsc_build_dir}/${petsc_local_arch_name_complex}/lib/petsc/conf/petscvariables not found])
					]
				)
			]
		)

		AC_LANG_POP
		export PETSC_DIR=${petsc_build_dir}
		export PETSC_ARCH_REAL=${petsc_local_arch_name_real}
		export PETSC_ARCH_COMPLEX=${petsc_local_arch_name_complex}
		cd $basedir
	],[AC_MSG_RESULT([no])
])

# check validity of PETSc installation
# first, environmental variables
AS_IF([test "x$PETSC_DIR" = "x"],[
	AC_MSG_FAILURE([No value provided for PETSC_DIR!])
])
AS_IF([test "x$PETSC_ARCH_REAL" = "x"],[
	AC_MSG_FAILURE([No value provided for PETSC_ARCH_REAL!])
])
AS_IF([test "x$PETSC_ARCH_COMPLEX" = "x"],[
	AC_MSG_FAILURE([No value provided for PETSC_ARCH_COMPLEX!])
])

# check directories
AC_CHECK_FILE( $PETSC_DIR,[
	AC_MSG_NOTICE([PETSc directory $PETSC_DIR OK.])
],[
	AC_MSG_ERROR([PETSc directory $PETSC_DIR not found.  Check your PETSC_DIR variable.])
])
AS_IF([test -d $PETSC_DIR],,[
	AC_MSG_ERROR([Specified PETSc directory ${PETSC_DIR} is not a directory!])
])

# Now see if the various include files are present
PETSC_INCLUDE_FILE_GENERIC="${PETSC_DIR}/lib/petsc/conf/variables"
AC_CHECK_FILE( $PETSC_INCLUDE_FILE_GENERIC,,[
	AC_MSG_ERROR([PETSc generic include file $PETSC_INCLUDE_FILE_GENERIC not found.  Your PETSc build may be incomplete."])
])
AC_CHECK_FILE( ${PETSC_DIR}/${PETSC_ARCH_REAL}/lib/petsc/conf/petscvariables,,[
	AC_MSG_ERROR([PETSc real architecture-specific include file ${PETSC_DIR}/${PETSC_ARCH_REAL}/lib/petsc/conf/petscvariables not found.  Check your PETSC_ARCH_REAL variable and your PETSc installation.])
])
AC_CHECK_FILE( ${PETSC_DIR}/${PETSC_ARCH_COMPLEX}/lib/petsc/conf/petscvariables,,
[
	AC_MSG_ERROR([PETSc complex architecture-specific include file ${PETSC_DIR}/${PETSC_ARCH_COMPLEX}/lib/petsc/conf/petscvariables not found.  Check your PETSC_ARCH_COMPLEX variable and your PETSc installation.])
])

AC_MSG_CHECKING([automatic local SLEPc installation])
AS_IF([test "x${with_localslepc}" = "xyes"],[
	AC_MSG_RESULT([yes])
	basedir=$(pwd)

	downloader=""
	AC_CHECK_PROGS(downloader,[wget curl])
	AC_MSG_CHECKING([command-line downloader])
	AS_IF([test "x${downloader}" = "x"],[
		AC_MSG_RESULT([no])
		AC_MSG_ERROR("Neither wget nor curl detected")
	],[
		AC_MSG_RESULT([$downloader])
	])

	# See if default SLEPc directory is overridden
	AS_IF([test "x${local_slepc_dir}" != "x"],
		[slepc_build_dir=${local_slepc_dir}],
		[slepc_build_dir="$(pwd)/extern/slepc"]
	)

	# See if default SLEPc version is overridden
	AS_IF([test "x$local_slepc_version" = "x"],
		[slepc_version=$default_slepc_version],
		[slepc_version=$local_slepc_version]
	)

	slepc_tarfile="slepc-${slepc_version}.tar.gz"
	AC_MSG_NOTICE([slepc_build_dir = ${slepc_build_dir}])
	AC_MSG_NOTICE([slepc_version = ${slepc_version}])
	AC_MSG_NOTICE([slepc_tarfile = ${slepc_tarfile}])

	# check for extra configuration options.
	slepc_extra_config=""

	# make the install dir
	AC_MSG_CHECKING([SLEPc build directory])
	AS_IF([test -d "$slepc_build_dir"],[
		AC_MSG_RESULT([yes])
	],[
		AC_MSG_RESULT([no, creating])
		mkdir -p "${slepc_build_dir}"
		AS_IF([test -d "$slepc_build_dir"],,[
			AC_MSG_ERROR(["Error creating SLEPc build directory ${slepc_build_dir}"])
		])
	])

	# see if it's already been downloaded and/or built
	cd ${slepc_build_dir}
	slepc_build_dir="${slepc_build_dir}/slepc-${slepc_version}"
	AC_MSG_CHECKING([SLEPc download])
	AS_IF([test -f "${slepc_build_dir}/configure"],
		[
			AC_MSG_RESULT([yes])
		],[
			AC_MSG_RESULT([no, installing])
			AC_CHECK_FILE([$slepc_tarfile],,[
				slepc_url="${slepc_url_base}/${slepc_tarfile}"
				downloader_flags=""
				AS_IF([test "x$downloader" = "xcurl"],[
					downloader_flags="--output ${slepc_tarfile}"
				])
				AS_IF([test "x$with_insecuressl" = "xyes"],[
					AS_IF([test "x$downloader" = "xwget"],[
						downloader_flags="${downloader_flags} --no-check-certificate"
					],[
						downloader_flags="${downloader_flags} --insecure"
					])
				])
				echo "Downloading using command: ${downloader} ${slepc_url} ${downloader_flags}"
				$downloader ${slepc_url} ${downloader_flags}
				AS_IF([test $? != 0],[
					AC_MSG_ERROR([Error downloading ${slepc_url}, check internet connection and version number])
				])
			])
			
			tar xvzf ${slepc_tarfile}
			AS_IF([test $? != 0],[
				AC_MSG_ERROR([Error decompressing slepc-${slepc_version}.tar.gz, check download])
			])
		]
	)
	export SLEPC_DIR=${slepc_build_dir}
	cd $SLEPC_DIR

	AC_CHECK_FILE( ${SLEPC_DIR}/${PETSC_ARCH_REAL}/lib/slepc/conf/slepcvariables,,[
		AC_MSG_RESULT([no, building])

		# build the real version with indicated options
		export PETSC_ARCH=${PETSC_ARCH_REAL}
		./configure
		
		AS_IF([test $? != 0],
			[AC_MSG_ERROR([Error running SLEPc configure command: ./configure])]
		)
		
		make PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} PETSC_ARCH=${PETSC_ARCH}
		AS_IF([test $? != 0],
			[AC_MSG_ERROR([Error running make command for real SLEPc: make PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} PETSC_ARCH=${PETSC_ARCH}])]
		)

		# see if it made the include file
		AC_MSG_CHECKING([successful ${PETSC_ARCH} SLEPc build])

		# clear the cached variable
		unset `$as_echo "ac_cv_file_${SLEPC_DIR}/${PETSC_ARCH_REAL}/lib/slepc/conf/slepcvariables" | $as_tr_sh`

		AC_CHECK_FILE( ${SLEPC_DIR}/${PETSC_ARCH_REAL}/lib/slepc/conf/slepcvariables,,[
			AC_MSG_ERROR([Error building real SLEPc, ${SLEPC_DIR}/${PETSC_ARCH_REAL}/lib/slepc/conf/slepcvariables not found])
		])
	])

	AC_MSG_CHECKING([for existing ${PETSC_ARCH_COMPLEX} SLEPc architecture])
	AC_CHECK_FILE( ${slepc_build_dir}/${PETSC_ARCH_COMPLEX}/lib/slepc/conf/slepcvariables,,[
		AC_MSG_RESULT([no, building])

		# build the complex version with indicated options
		export PETSC_ARCH=${PETSC_ARCH_COMPLEX}
		./configure
		
		AS_IF([test $? != 0],
			[AC_MSG_ERROR([Error running SLEPc configure command: ./configure])]
		)
		
		make PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR}
		AS_IF([test $? != 0],
			[AC_MSG_ERROR([Error running make command for complex SLEPc: make PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR}])]
		)

		# see if it made the include file
		AC_MSG_CHECKING([successful ${PETSC_ARCH} SLEPc build])

		# clear the cached variable
		unset `$as_echo "ac_cv_file_${slepc_build_dir}/${PETSC_ARCH_COMPLEX}/lib/slepc/conf/slepcvariables" | $as_tr_sh`

		AC_CHECK_FILE( ${slepc_build_dir}/${PETSC_ARCH_COMPLEX}/lib/slepc/conf/slepcvariables,,[
			AC_MSG_ERROR([Error building complex SLEPc, ${slepc_build_dir}/${PETSC_ARCH_COMPLEX}/lib/slepc/conf/slepcvariables not found])
		])
	])

	cd $basedir
],[
	AC_MSG_RESULT([no])
])


# Checks for environmental variables
AS_IF([test "x$SLEPC_DIR" = "x"],[
	AC_MSG_ERROR([No value provided for SLEPC_DIR!])
])

# Check the SLEPc installation
# First, check that the directories exist
AC_CHECK_FILE( $SLEPC_DIR,[
	AC_MSG_NOTICE([SLEPc directory $SLEPC_DIR OK])
],[
	AC_MSG_ERROR([SLEPc directory $SLEPC_DIR not found.  Check your SLEPC_DIR variable.])
])

AS_IF([test ! -d $SLEPC_DIR],[
	AC_MSG_ERROR([Specified PETSc directory ${PETSC_DIR} is not a directory!])
])

# Now see if the generic include files are present
SLEPC_INCLUDE_FILE_GENERIC="${SLEPC_DIR}/lib/slepc/conf/slepc_common"
AC_CHECK_FILE( $SLEPC_INCLUDE_FILE_GENERIC,[
	AC_MSG_NOTICE([SLEPc generic include file OK])
],[
	AC_MSG_ERROR([SLEPc generic include file $SLEPC_INCLUDE_FILE_GENERIC not found.  Your SLEPc build may be incomplete.])
])

# Finally check the architecture-specific files
AC_CHECK_FILE( ${SLEPC_DIR}/${PETSC_ARCH_REAL}/lib/slepc/conf/slepcvariables,[
	AC_MSG_NOTICE([SLEPc real architecture-specific include file OK])
],[
	AC_MSG_ERROR([SLEPc real architecture-specific include file ${SLEPC_DIR}/${PETSC_ARCH_REAL}/lib/slepc/conf/slepcvariables not found.  Check your PETSC_ARCH_REAL variable and your SLEPc installation.])
])
AC_CHECK_FILE( ${SLEPC_DIR}/${PETSC_ARCH_COMPLEX}/lib/slepc/conf/slepcvariables,[
	AC_MSG_NOTICE([SLEPc complex architecture-specific include file OK])
],[
	AC_MSG_ERROR([SLEPc complex architecture-specific include file ${SLEPC_DIR}/${PETSC_ARCH_COMPLEX}/lib/slepc/conf/slepcvariables not found.  Check your PETSC_ARCH_COMPLEX variable and your SLEPc installation.])
])


# Set output variables to propagate into Makefiles
#AC_SUBST([INCLUDEFLAGS],"-I. -I../common -I../units -I../attenuation -I../atmosphere -I../libmodes -I../libpe -I../libbroadband -I/usr/local/include -I/usr/include")
AC_SUBST([INCLUDEFLAGS],"-I../include -I/usr/local/include -I/usr/include")
AC_SUBST([LIBS],"-lgsl -lgslcblas -lm -lfftw3")
AC_SUBST([WARNINGFLAGS],${warning_compiler_flags})
AC_SUBST([CXXFLAGS], "-fpic -c -Wall $CXXFLAGS")
LDFLAGS="${LDFLAGS} -L/usr/lib"
AC_SUBST([STATICLIBS],"../../lib/libatmosphere.a ../../lib/libunits.a ../../lib/libcommon.a")
AC_SUBST([MODELIBS],"../../lib/libmodes.a")
AC_SUBST([PELIBS],"../../lib/libpe.a")
AC_SUBST([BROADBANDLIBS],"../../lib/libbroadband.a")
AC_SUBST([LDFLAGS],$LDFLAGS)
AC_SUBST([PETSC_DIR],$PETSC_DIR)
AC_SUBST([SLEPC_DIR],$SLEPC_DIR)
AC_SUBST([PETSC_INCLUDE_FILE_GENERIC],$PETSC_INCLUDE_FILE_GENERIC)
AC_SUBST([SLEPC_INCLUDE_FILE_GENERIC],$SLEPC_INCLUDE_FILE_GENERIC)

AC_SUBST([PETSC_ARCH_REAL],$PETSC_ARCH_REAL)
AC_SUBST([PETSC_ARCH_COMPLEX],$PETSC_ARCH_COMPLEX)


AC_CONFIG_FILES([
		Makefile
		src/common/Makefile
		src/common/test/Makefile
		src/units/Makefile
		src/units/test/Makefile
		src/atmosphere/Makefile
		src/libmodes/Makefile
		src/libpe/Makefile
		src/libbroadband/Makefile
		src/modess/Makefile
		src/modbb/Makefile
		src/epade_pe/Makefile
		src/wmod/Makefile
		src/gfpe/Makefile
])


AC_OUTPUT

AC_MSG_NOTICE([
	Configuration complete.  Run 'make' to finish building $PACKAGE_TARNAME.
])
AC_MSG_NOTICE

cat <<EOF

Configuration complete.  Run 'make' to finish building $PACKAGE_TARNAME.

Add the following lines to your .bashrc, .bash_profile, or .profile login file:

export PETSC_DIR=${PETSC_DIR}
export SLEPC_DIR=${SLEPC_DIR}
export PETSC_ARCH_REAL=${PETSC_ARCH_REAL}
export PETSC_ARCH_COMPLEX=${PETSC_ARCH_COMPLEX}

EOF
