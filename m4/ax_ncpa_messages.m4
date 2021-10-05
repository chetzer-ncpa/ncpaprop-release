
# AX_NCPA_MESSAGES_GXX_NOT_FOUND
# Outputs a customized error messages with instructions for 
# installing g++ depending on the detected OS.
#
AC_DEFUN([AX_NCPA_MESSAGES_GXX_NOT_FOUND],[
	AC_REQUIRE([AC_CANONICAL_HOST])
	case $host_os in
		linux* )
			cat <<EOF

Configure failed to detect g++, the GNU C++ compiler.  This is required to compile $PACKAGE_NAME.  Please install the 'g++' package or compile from source.  Depending on your package manager, this may be accomplished with 
   sudo apt-get install g++
or
   sudo yum install gcc-c++

EOF
			AC_MSG_ERROR([missing g++ compiler])
			;;

		darwin* )
			cat <<EOF

Configure failed to detect g++, the GNU C++ compiler.  This is required to compile $PACKAGE_NAME.  To obtain Apple's base version, install the XCode package from the App Store, then from within XCode install the Command Line Developer Tools.  This will provide both gcc and g++.

EOF
			AC_MSG_ERROR([missing g++ compiler])
			;;

		*)
			AC_MSG_ERROR([operating system not recognized.])
			;;
	esac
])



# AX_NCPA_MESSAGES_LFFTW3_NOT_FOUND
# Outputs a customized error message with instructions for installing
# fftw3, depending on the detected OS.
#
AC_DEFUN([AX_NCPA_MESSAGES_LFFTW3_NOT_FOUND],[
        AC_REQUIRE([AC_CANONICAL_HOST])
        case $host_os in
                linux* )
                        cat <<EOF

Configure failed to detect libfftw3, the FFTW version 3 library.  This is required to compile $PACKAGE_NAME.  Please install the fftw3 dev package or compile from source.  Depending on your package manager, this may be accomplished with
   sudo apt-get install libfftw3-dev
or
   sudo yum install fftw-devel

EOF
                        AC_MSG_ERROR([missing libfftw3])
                        ;;

                darwin* )
                        cat <<EOF

Configure failed to detect libfftw3, the FFTW version 3 library.  This is required to compile $PACKAGE_NAME.  Please install using your package manager of choice or from source.

EOF
                        AC_MSG_ERROR([missing libfftw3])
                        ;;

                *)
                        AC_MSG_ERROR([operating system not supported.])
                        ;;
        esac
])

# AX_NCPA_MESSAGES_LGSL_NOT_FOUND
# Outputs a customized error message with instructions for installing
# GSL, depending on the detected OS.
#
AC_DEFUN([AX_NCPA_MESSAGES_LGSL_NOT_FOUND],[
        AC_REQUIRE([AC_CANONICAL_HOST])
        case $host_os in
                linux* )
                        cat <<EOF

Configure failed to detect libgsl, the GSL library.  This is required to compile $PACKAGE_NAME.  Please install the GSL dev package or compile from source.  Depending on your package manager, this may be accomplished with
   sudo apt-get install libgsl0-dev
or
   sudo yum install gsl-devel

EOF
                        AC_MSG_ERROR([missing libgsl])
                        ;;

                darwin* )
                        cat <<EOF

Configure failed to detect libgsl, the GSL library.  This is required to compile $PACKAGE_NAME.  Please install using your package manager of choice.

EOF
                        AC_MSG_ERROR([missing libgsl])
                        ;;

                *)
                        AC_MSG_ERROR([operating system not supported.])
                        ;;
        esac
])

