dnl ======================================================================
dnl autoconf macros originally developed for the BBH project
dnl
dnl Copyright 1995--2001 Matthew William Choptuik
dnl The University of British Columbia
dnl               and                    
dnl The University of Texas at Austin
dnl
dnl Please AVOID making modified versions of this file.  If ammendments 
dnl are necessary, send e-mail to 
dnl
dnl     choptuik@physics.ubc.ca
dnl 
dnl immediately.  
dnl
dnl ======================================================================
dnl NOTES:
dnl
dnl (1) April 6 2000
dnl 
dnl     setenv  BBH_CHECK_DEFAULTS NONE
dnl 
dnl     will (should) inhibit use of default lib/include search locations
dnl
dnl (2) July 18 2000
dnl 
dnl     Following Scott Hawley's lead, will define 'CCF77LIBS' to be the 
dnl     standard set of libraries needed to link generic Fortran 77 code
dnl     with C compiler.
dnl
dnl (3) January 27 2001
dnl
dnl     Removed '-mips3' settings for SGIs 
dnl
dnl (4) April 22 2001,   
dnl  
dnl     Fixed bug whereby 'configure' invoked without explicit prefix 
dnl     setting results in prefix=NONE (ac_default_prefix=/usr/local)
dnl     (Thanks to Reinoud Slagter for unearthing this one.)
dnl
dnl (4) May 6 2001,   
dnl  
dnl     Added '-lm' to CCF77LIBS
dnl
dnl (5) December 1 2001
dnl     
dnl     Removed 'glut GLU Xmu Xi' from LIBGL settings, and added LIBGLUT
dnl
dnl (5) December 7 2001
dnl     
dnl     DAGH_NO_MPI environment var now sets -DACE_NO_MPI
dnl
dnl (6) July 31 2003    
dnl     
dnl     Discovered that there is a '-fno-second-underscore' option to g77, 
dnl     changing defaults so that all environments use single-underscore
dnl     convention, can be overridden by setting DOUBLE_UNDERSCORE, e.g.
dnl   
dnl     setenv DOUBLE_UNDERSCORE
dnl ======================================================================

dnl ----------------------------------------------------------------------
dnl BBH_SYS_GETSYSTEM attempts to determine what flavour of UNIX
dnl is being used
dnl ----------------------------------------------------------------------

dnl ----------------------------------------------------------------------
dnl Template for new macros
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_TEMPLATE,[
])

dnl ----------------------------------------------------------------------
dnl BBH_SYS_GETSYSTEM attempts to determine what flavour of UNIX
dnl is being used
dnl 
dnl Sets Output Variable BBH_SYSTEM and possibly others, dependent on 
dnl the pecularities of the given OS.
dnl
dnl Modified Nov 3 1999, to include definitions of CPPFLAGS, which
dnl DAGH/GRACE now wants.
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_SYS_GETSYSTEM,[
   if test "X$BBH_SYSTEM" = "X"; then
dnl   Ensure 'prefix' is sanely set
      if test "X${prefix}" = "XNONE"; then
         prefix=$ac_default_prefix
      fi
      BBH_GET_NPROC
      BBH_SUBSYSTEM="NONE"
		if test "$prefix" != "NONE"; then
			 BBH_PUSHUNIQUE(INCLUDE_PATHS,$prefix/include)
			 BBH_PUSHUNIQUE(LIB_PATHS,$prefix/lib)
		fi
      RPCGEN_CPPFLAGS=""
      AC_MSG_CHECKING(for Unix flavour)
      AC_CACHE_VAL(ac_cv_bbh_system,
      [
         changequote(,)dnl
         __UNAME=`uname -a | tr '[a-z]' '[A-Z]'`
         changequote([, ])dnl
         for s in $__UNAME
         do
            case $s in
            IRIX64)                       
dnl ----------------------------------------------------------------------
dnl              Tried using this flag setting on SGI's to rationalize
dnl              signal behaviour.  Needs further investigation.
dnl
dnl              CFLAGS="$CFLAGS -D_BSD_SIGNALS"
dnl ----------------------------------------------------------------------
dnl              DAGH: Masking CXX warnings 
dnl                    1233 
dnl                    1021  (type qualifiers are meaningless in this 
dnl                           declaration)
dnl                    3262  (vbl. declared and never referenced)
dnl              Disabled this as well since wasting too much time 
dnl              Just set CXXFLAGS to include appropriate -woff string, e.g.
dnl 
dnl              setenv CXXFLAGS "-woff 1233,1021,3262"
dnl              setenv CXXFLAGS "$CXXFLAGS -woff 1233,1021,3262"
dnl ----------------------------------------------------------------------
dnl              CXXFLAGS="$CXXFLAGS -ptused -DSGI -woff 1233,1021,3262"
dnl ----------------------------------------------------------------------
                 CXXFLAGS="$CXXFLAGS -ptused -DSGI"
                 CFLAGS="$CFLAGS -woff 1174,1178,1552"
                 CPPFLAGS="$CPPFLAGS -DSGI"
                 CXXREPOSITORY="ii_files"
dnl ----------------------------------------------------------------------
dnl              Linking to Fortran 90 code via C++ will generally 
dnl              require '-lftn -lftn90'
dnl ----------------------------------------------------------------------
                 CXXF90LIBS='-lftn -lftn90'
dnl ----------------------------------------------------------------------
dnl              For linkage to Fortran 77 code via C
dnl ----------------------------------------------------------------------
                 if test "X$CCF77LIBS" = "X"; then
                    CCF77LIBS='-lftn -lm'
                 fi
dnl ----------------------------------------------------------------------
dnl              For GL/OpenGL
dnl ----------------------------------------------------------------------
                 if test "X$LIBGL" = "X"; then
                    LIBGL="-lGLU -lGL -lXext -lX11 -lm"
                    LIBGLUT="-lglut -lGLU -lGL -lXmu -lXi -lXext -lX11 -lm"
                 fi
dnl ----------------------------------------------------------------------
dnl              Can't get load warnings turned off satisfactorily
dnl               LDFLAGS="$LDFLAGS -Wl,\"-woff 15,84\""
dnl              LDFLAGS="$LDFLAGS -Wl,\"-wall\""
dnl              echo "LDFLAGS: $LDFLAGS"
dnl ----------------------------------------------------------------------
                 case "X$IRIXMODE" in
                 X32)
                     BBH_SYSTEM=IRIX32
                     CFLAGS="$CFLAGS -DSGI -32";
                     CXXFLAGS="CXXFLAGS -DSGI -32";
                     LDFLAGS="$LDFLAGS -32";
                 ;;
                 Xn32|XN32)
                     BBH_SYSTEM=IRIXN32
                     CFLAGS="$CFLAGS -DSGI -n32"
                     CXXFLAGS="$CXXFLAGS -DSGI -n32";
                     LDFLAGS="$LDFLAGS -n32";
                 ;;
                 *)
                     BBH_SYSTEM=IRIX64;        
                     CFLAGS="$CFLAGS -DSGI -64";
                     CXXFLAGS="$CXXFLAGS -DSGI -64";
                     LDFLAGS="$LDFLAGS -64";
                 ;;
                 esac

                 case $BBH_NPROC in
                 1)  ;;
                 *)  
                     CFLAGS="$CFLAGS -WK,-listoption=O -pfa -mp";
                     CXXFLAGS="$CXXFLAGS -pfa -mp";
                     F77FLAGS="$F77FLAGS -pfa -mp";
                     F77FLAGSNOOPT="$F77FLAGSNOOPT -pfa -mp";
                     LDFLAGS="$LDFLAGS -pfa -mp";
                 ;;
                 esac;
                 F77FLAGS="$F77FLAGS -NC149999"
                 CC_TRANSFORM="cat"
dnl                                       F77_TRANSFORM="touch";
                                                                  break;;
            IRIX*)                        BBH_SYSTEM=IRIX;        
dnl                                       F77_TRANSFORM="touch";
                 CFLAGS="$CFLAGS -DSGI"
                 CPPFLAGS="$CPPFLAGS -DSGI"
                 CC_TRANSFORM="cat"
dnl ----------------------------------------------------------------------
dnl              See comment above
dnl              CXXFLAGS="$CXXFLAGS -ptused -DSGI -woff 1233";
dnl ----------------------------------------------------------------------
                 CXXFLAGS="$CXXFLAGS -ptused -DSGI"
                 CXXREPOSITORY="ii_files"
                 CXXF90LIBS='-lftn -lftn90'
                 case "X$IRIXMODE" in
                 Xn32|XN32)
                     BBH_SYSTEM=IRIXN32
                     CFLAGS="$CFLAGS -n32"
                     CXXFLAGS="$CXXFLAGS -n32"
                     LDFLAGS="$LDFLAGS -n32";
                 ;;
                 X64)
                     BBH_SYSTEM=IRIX64
                     CFLAGS="$CFLAGS -64"
                     CXXFLAGS="$CXXFLAGS -64";
                     LDFLAGS="$LDFLAGS -64";
                 ;;
                 esac
                                                                  break;;
            SUNOS)                        
                 BBH_SYSTEM=SUNOS;       
dnl              F77_TRANSFORM="touch";
                 LDFLAGS="$LDFLAGS -L/usr/ccs/lib ";               
                 CXXFLAGS="$CXXFLAGS -pto"
                 CPPFLAGS="$CPPFLAGS -xCC -DCPPFLAGS_UNKNOWN"
                 CXXREPOSITORY="Templates.DB"
                 CXXF90LIBS='-lsunmath -lF77'
                 CC_TRANSFORM="cat"
dnl ----------------------------------------------------------------------
dnl              For linkage to Fortran 77 code via C
dnl ----------------------------------------------------------------------
                 if test "X$CCF77LIBS" = "X"; then
                    CCF77LIBS='-lF77 -lM77 -lsunmath -lm'
                 fi
                                                                  break;;
            AIX)                          
                 BBH_SYSTEM=AIX;        
dnl              F77_TRANSFORM="touch";
                 CXXF90LIBS='-lxlf -lxlf90'
                 CC_TRANSFORM="cat"
dnl ----------------------------------------------------------------------
dnl              Check to see whether this is Cornell SP-2
dnl ----------------------------------------------------------------------
dnl              for ss in $__UNAME
dnl              do
dnl                 case $ss in
dnl                 SP*) BBH_SUBSYSTEM="SP2";;
dnl                 *);;
dnl                 esac
dnl              done
                 BBH_SUBSYSTEM=""
                 case $BBH_SUBSYSTEM in
                 SP2) 
                    BBH_SUBSYSTEM="SP2"
                    CFLAGS="$CFLAGS -DSPX"
                    CPPFLAGS="$CPPFLAGS -DCPPFLAGS_UNKNOWN"
                    CXXFLAGS="$CXXFLAGS -DSPX"
                    CC=mpcc
                    CXX=mpCC
                 ;;
                 *)
                    CFLAGS="$CFLAGS -DRS6000 -D_ALL_SOURCE"
                    CPPFLAGS="$CPPFLAGS -DCPPFLAGS_UNKNOWN"
                    CXXFLAGS="$CXXFLAGS -DRS6000"
                    if test "X$CC" = "X"; then
                       CC=xlc
                    fi
                    if test "X$CXX" = "X"; then
                       CXX=xlC
                    fi
dnl ----------------------------------------------------------------------
dnl  This works on the TICAM RS/6000's.  The path may be wrong on others.
dnl ----------------------------------------------------------------------
dnl                         LDFLAGS="$LDFLAGS -bI:/usr/lpp/xlf/lib/lowsys.exp";

dnl ----------------------------------------------------------------------
dnl                 These flags are DAGH-specific; if needed (at this 
dnl                 time they are apparently *not needed) they should
dnl                 probably be set somewhere else.
dnl ----------------------------------------------------------------------
dnl                 CFLAGS="$CFLAGS -DMPICH -D_ALL_SOURCE"
dnl                 CFLAGS="$CFLAGS -DFORTRANNOUNDERSCORE -DMPI_rs6000"
dnl                 CXXFLAGS="$CXXFLAGS -DMPICH -D_ALL_SOURCE"
dnl                 CXXFLAGS="$CXXFLAGS -DFORTRANNOUNDERSCORE -DMPI_rs6000"
dnl ----------------------------------------------------------------------
                 ;;
                 esac
                 CXXREPOSITORY="tempinc"
                                                                  break;;
            A/UX)                         
                 BBH_SYSTEM=A/UX;       
dnl              F77_TRANSFORM="touch";
                 CC_TRANSFORM="cat"
                                                                  break;;
            OSF1)                         
                 BBH_SYSTEM=OSF1;        
                 F77_TRANSFORM="touch";
                 CC_TRANSFORM="decomment-c-c++"
                 CPPFLAGS="$CPPFLAGS -std1 -DDEC_ALPHA";
                 CXXFLAGS="$CXXFLAGS -x cxx"
                 BBH_DEFS="$BBH_DEFS -std1 -DDEC_ALPHA";
                 if test "X$CCF77LIBS" = "X"; then
                    CCF77LIBS='-lfor -lFutil -lots -lUfor -lm';
                 fi
                                                                  break;;
dnl ----------------------------------------------------------------------
dnl         CRAYS 
dnl ----------------------------------------------------------------------
            UNICOS|CRAY)                  
                 BBH_SYSTEM=UNICOS;      
dnl              F77_TRANSFORM="../bin/f77transcray";
                 BBH_RNPL_FLAGS="-csufE";
                 CPPFLAGS="$CPPFLAGS -D_CRAY"
                 CC_TRANSFORM="cat"
dnl ----------------------------------------------------------------------
dnl              Get subsystem (C90, J90, T3E, T90, SV1)
dnl ----------------------------------------------------------------------
                 for ss in $__UNAME; do
                    case $ss in 
                    C90*)   BBH_SUBSYSTEM="C90";
                            CXXFLAGS="$CXXFLAGS -DCRAYC90"
                            CPPFLAGS="$CPPFLAGS -DCRAYNOIEEE"
                                                    break;;
                    J90*)   BBH_SUBSYSTEM="J90";
                            CXXFLAGS="$CXXFLAGS -DCRAYJ90"
                            CPPFLAGS="$CPPFLAGS -DCRAYNOIEEE"
                                                    break;;
                    T3E*)   BBH_SUBSYSTEM="T3E";
dnl ----------------------------------------------------------------------
dnl                         CXXFLAGS="$CXXFLAGS -Tcray-t3e -DCRAYT3D";
dnl ----------------------------------------------------------------------
                            CXXFLAGS="$CXXFLAGS -DCRAYT3D";
                            CXXFLAGS="$CXXFLAGS -h instantiate=local";
                            CPPFLAGS="$CPPFLAGS -DCPPFLAGS_UNKNOWN"
                                                    break;;
                    T90*|T94*)   BBH_SUBSYSTEM="T90";
dnl ----------------------------------------------------------------------
dnl                         CXXFLAGS="$CXXFLAGS -Tcray-t3e -DCRAYT3D";
dnl ----------------------------------------------------------------------
                            CXXFLAGS="$CXXFLAGS -DCRAYT90";
                            CXXFLAGS="$CXXFLAGS -h instantiate=local";
                            CPPFLAGS="$CPPFLAGS -DCPPFLAGS_UNKNOWN"
                                                    break;;
                    SV1*)   BBH_SUBSYSTEM="SV1";
dnl ----------------------------------------------------------------------
dnl                         CXXFLAGS="$CXXFLAGS -Tcray-t3e -DCRAYT3D";
dnl ----------------------------------------------------------------------
                            CXXFLAGS="$CXXFLAGS -DCRAYT90";
                            CXXFLAGS="$CXXFLAGS -h instantiate=local";
                            CPPFLAGS="$CPPFLAGS -DCPPFLAGS_UNKNOWN"
                                                    break;;
                    *)                                   ;;
                    esac
                 done
dnl ----------------------------------------------------------------------
dnl              Detect mario
dnl ----------------------------------------------------------------------
                 if test `hostname` = "mario"; then
                    BBH_SYSTEM_MARIO="yes";
                 fi
                                                                  break;;
dnl ----------------------------------------------------------------------
dnl         LINUX
dnl ----------------------------------------------------------------------
            LINUX)                        
dnl ----------------------------------------------------------------------
dnl         Use Fortran compiler to fix BBH_SYSTEM
dnl ----------------------------------------------------------------------
                 case "X$F77" in
                 Xpgf77)  
                    BBH_SYSTEM=LINUX_PG
                    F77_TRANSFORM="touch";
						  CC_TRANSFORM="cat"
                    BBH_DEFS="$BBH_DEFS -DLINUX_PG"
                    CFLAGS="$CFLAGS"
dnl                 CPPFLAGS="$CPPFLAGS -DLINUX -D__ -DLINUX_PG -DPORT_GROUP -DWant_c_files"
dnl                 CXXFLAGS="$CXXFLAGS -DLINUX -D__ -DLINUX_PG -DPORT_GROUP -DWant_c_files"
                     CPPFLAGS="$CPPFLAGS -DLINUX -DLINUX_PG -DPORT_GROUP -DWant_c_files"
                     CXXFLAGS="$CXXFLAGS -DLINUX -DLINUX_PG -DPORT_GROUP -DWant_c_files"
                    if test "X$CCF77LIBS" = "X"; then
                       CCF77LIBS="-lpgftnrtl -lpgf90rtl -lm"
                    fi
                 ;;
                 Xefc)  
                    BBH_SYSTEM=LINUX_IA64
                    F77_TRANSFORM="touch";
						  CC_TRANSFORM="cat"
                    BBH_DEFS="$BBH_DEFS"
                    CFLAGS="$CFLAGS"
                    CPPFLAGS="$CPPFLAGS -DLINUX -DPORT_GROUP -DWant_c_files"
                    CXXFLAGS="$CXXFLAGS -DLINUX -DPORT_GROUP -DWant_c_files"
                    if test "X$CCF77LIBS" = "X"; then
                       CCF77LIBS="-lCEPCF90 -lF90 -lintrins -lm"
                    fi
                 ;;
                 Xifc)  
                    BBH_SYSTEM=LINUX_IA32
                    F77_TRANSFORM="touch";
						  CC_TRANSFORM="cat"
                    BBH_DEFS="$BBH_DEFS"
                    CFLAGS="$CFLAGS"
                    CPPFLAGS="$CPPFLAGS -DLINUX -DPORT_GROUP -DWant_c_files"
                    CXXFLAGS="$CXXFLAGS -DLINUX -DPORT_GROUP -DWant_c_files"
                    if test "X$CCF77LIBS" = "X"; then
                       CCF77LIBS="-lCEPCF90 -lF90 -lintrins -lm"
                    fi
                 ;;
                 *) 
                    BBH_SYSTEM=LINUX
                    BBH_SETDEF(CC,gcc)
                    BBH_SETDEF(CXX,gcc)
                    BBH_SETDEF(CFLAGS,-O3)
                    BBH_SETDEF(CPPFLAGS,"")
                    BBH_SETDEF(F77,g77)
                    if test "X$DOUBLE_UNDERSCORE" = "X"; then
                       BBH_SETDEF(F77FLAGS,"-O3 -fno-second-underscore")
                       BBH_SETDEF(F77LFLAGS,"-O3 -fno-second-underscore")
                       BBH_SETDEF(F90FLAGS,"-O3 -fno-second-underscore")
                    else
                       BBH_SETDEF(F77FLAGS,-O3)
                       BBH_SETDEF(F77LFLAGS,-O3)
                       BBH_SETDEF(F90FLAGS,-O3)
                    fi
                    BBH_SETDEF(F90,g77)
                    F77_TRANSFORM="touch";
						  CC_TRANSFORM="cat"
                    BBH_DEFS="$BBH_DEFS -DLINUX"
                    CFLAGS="$CFLAGS"
dnl                 CPPFLAGS="$CPPFLAGS -DLINUX -D__ -DPORT_GROUP -DWant_c_files"
dnl                 CXXFLAGS="$CXXFLAGS -DLINUX -D__ -DPORT_GROUP -DWant_c_files"
                    CPPFLAGS="$CPPFLAGS -DLINUX -DPORT_GROUP -DWant_c_files"
                    CXXFLAGS="$CXXFLAGS -DLINUX -DPORT_GROUP -DWant_c_files"
                    if test "X$CCF77LIBS" = "X"; then
                       CCF77LIBS='-lgfortran -lm'
                    fi
                 ;;
                 esac;
                 RPCGEN_CPPFLAGS="-U__STDC__"
dnl ----------------------------------------------------------------------
dnl              For GL/OpenGL
dnl ----------------------------------------------------------------------
                 LDFLAGS="$LDFLAGS -L/usr/X11R6/lib"
                 if test "X$LIBGL" = "X"; then
                    LIBGL="-lGLU -lGL -lXext -lX11 -lm"
                    LIBGLUT="-lglut -lGLU -lGL -lXmu -lXi -lXext -lX11 -lm"
                 fi
                                                                  break;;
dnl ----------------------------------------------------------------------
dnl         DARWIN
dnl ----------------------------------------------------------------------
            DARWIN)                        
                 BBH_SYSTEM=DARWIN
                 BBH_SETDEF(CC,gcc)
                 BBH_SETDEF(CXX,gcc)
                 BBH_SETDEF(CFLAGS,-O3)
                 BBH_SETDEF(CPPFLAGS,"")
                 BBH_SETDEF(F77,g77)
                 BBH_SETDEF(F77FLAGS,-O3)
                 BBH_SETDEF(F77LFLAGS,-O3)
                 BBH_SETDEF(F90,g77)
                 BBH_SETDEF(F90FLAGS,-O3)
                 F77_TRANSFORM="touch";
                 CC_TRANSFORM="cat"
                 CFLAGS="$CFLAGS"
                 CPPFLAGS="$CPPFLAGS -DDARWIN -D__ -DWant_c_files"
                 CXXFLAGS="$CXXFLAGS -DDARWIN -D__ -DWant_c_files"
                 if test "X$CCF77LIBS" = "X"; then
                    CCF77LIBS='-lgfortran -lm'
                 fi
                 RPCGEN_CPPFLAGS="-U__STDC__"
dnl ----------------------------------------------------------------------
dnl              For GL/OpenGL and g77
dnl ----------------------------------------------------------------------
                 LDFLAGS="$LDFLAGS -L/usr/X11R6/lib -L/sw/lib"
                 if test "X$LIBGL" = "X"; then
                    LIBGL="-lGLU -lGL -lXext -lX11 -lm"
                    LIBGLUT="-lglut -lGLU -lGL -lXmu -lXi -lXext -lX11 -lm"
                 fi
                                                                  break;;
            HP-UX)                        
                 BBH_SYSTEM=HP-UX;       
dnl              F77_TRANSFORM="touch";
                 CC_TRANSFORM="cat";
                 BBH_DEFS="$BBH_DEFS -DHP9000";
                 CFLAGS="$CFLAGS -Aa -D_HPUX_SOURCE -w";        
                 CPPFLAGS="$CPPFLAGS -DCPPFLAGS_UNKNOWN"
                                                                  break;;
            *)                            
                 BBH_SYSTEM=UNKNOWN           ;;
            esac
         done
         ac_cv_bbh_system="$BBH_SYSTEM"
dnl      ac_cv_bbh_f77_transform="$F77_TRANSFORM";
      ]
      )

      case $BBH_NPROC in
      1) AC_MSG_RESULT($ac_cv_bbh_system);;
      *) AC_MSG_RESULT("$BBH_NPROC-processor $ac_cv_bbh_system");;
      esac

      case $BBH_SUBSYSTEM in
      NONE) ;;
      *) AC_MSG_CHECKING("for Unix subflavour"); AC_MSG_RESULT($BBH_SUBSYSTEM);;
      esac

      BBH_SYSTEM=$ac_cv_bbh_system
dnl   F77_TRANSFORM=$ac_cv_bbh_f77_transform

dnl ----------------------------------------------------------------------
dnl   Export some values
dnl ----------------------------------------------------------------------
      AC_SUBST(BBH_SYSTEM)
      AC_SUBST(BBH_RNPL_FLAGS)
      AC_SUBST(CC_TRANSFORM)
dnl   AC_SUBST(F77_TRANSFORM)
dnl   AC_SUBST(CPPFLAGS)
      AC_SUBST(CXXREPOSITORY)
      AC_SUBST(CXXFLAGS)
      AC_SUBST(CXXF90LIBS)
      AC_SUBST(CCF77LIBS)
      AC_SUBST(LIBGL)
      AC_SUBST(LIBGLUT)
      AC_SUBST(RPCGEN_CPPFLAGS)

      case $BBH_SYSTEM in
      SUNOS)
         AC_CHECK_LIB(socket, socket,
           LIBS="$LIBS -lsocket"
           BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -lsocket"
           BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -lsocket"
           )
         AC_CHECK_LIB(nsl, clnt_create,
           LIBS="$LIBS -lnsl"
           BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -lnsl"
           BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -lnsl"
           )
                                                                  break;;
dnl ----------------------------------------------------------------------
dnl Used to set "-lbsd" on AIX in attempt to get ^C handling working
dnl more  consistently, but that breaks 'pow' and probably other things.
dnl ----------------------------------------------------------------------
      esac
   fi
])

dnl ----------------------------------------------------------------------
dnl BBH_GET_NPROC attempts to determine the number of processors in
dnl the system.
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_GET_NPROC,[
   if test "X$BBH_NPROC" = "X"; then
      AC_CACHE_VAL(ac_cv_bbh_nproc,
      [
         changequote(,)dnl
         __UNAME=`uname -a | tr '[a-z]' '[A-Z]'`
         for s in $__UNAME
         do
            case $s in
            IRIX64) ac_cv_bbh_nproc=`hinv | grep 'IP.*Processor' | sed 's? .*??' | grep '^[0-9]*$'`;
                                                                     break;;
            *)      ac_cv_bbh_nproc=1;
                                                                     break;;
            esac
         done
         changequote([, ])dnl
      ])
      BBH_NPROC=$ac_cv_bbh_nproc
      AC_SUBST(BBH_NPROC)
   fi
])


dnl ----------------------------------------------------------------------
dnl Locate lib, include paths for DAGH
dnl
dnl "Standard" places people install DAGH can be hard-coded here
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_CHECK_DAGH,[
dnl ----------------------------------------------------------------------
dnl   Look for DAGH header ...
dnl ----------------------------------------------------------------------
   BBH_CHECK_FATAL="no"

   AC_MSG_CHECKING(for DAGH header)
   if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
      BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS /usr/local/include"
      BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS $HOME/dagh/include $HOME/include"
   fi
   BBH_FINDFILE(BBH_CHECK_INC_PATHS,"DAGH.h",BBH_INC_PATH)
   if test $BBH_INC_PATH != "no"; then
      BBH_PUSHUNIQUE(BBH_RNPLAPP_INCPATHS,$BBH_INC_PATH)
      AC_MSG_RESULT($BBH_INC_PATH/DAGH.h)
   else
      AC_MSG_RESULT(not found ... will have to abort configuration)
      BBH_CHECK_FATAL="yes"
   fi
 
   AC_MSG_CHECKING(for DAGH library)
   if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /usr/local/lib"
      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS $HOME/dagh/lib $HOME/lib"
   fi
   BBH_FINDLIB(BBH_CHECK_LIB_PATHS,"dagh",BBH_LIB_PATH)
   if test $BBH_LIB_PATH != "no"; then
      BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)
      BBH_RNPLAPP_FLIBS="-ldagh $BBH_RNPLAPP_FLIBS"
      BBH_RNPLAPP_CLIBS="-ldagh $BBH_RNPLAPP_CLIBS"
      BBH_RNPLBLD_CLIBS="-ldagh $BBH_RNPLBLD_CLIBS"
      AC_MSG_RESULT($BBH_LIB_PATH/libdagh.a)
   else
      AC_MSG_RESULT(not found ... will have to abort configuration)
      BBH_CHECK_FATAL="yes"
   fi
   if test $BBH_CHECK_FATAL = "yes"; then
      BBH_MSG_DAGH_EXAMPLES_NEED_DAGH
   fi

])

AC_DEFUN(BBH_MSG_DAGH_EXAMPLES_NEED_DAGH,[
cat<<.

Installation of the DAGH examples directory requires that DAGH 
be intalled on your system and that this script be able to locate 
the header file 'DAGH.h' and the library file 'libdagh.a'.
One or both of these files could not be found in the "usual places"
If you know the directories in which the 
headers and library are installed on this system, set the 
environment variables

LIB_PATHS
INCLUDE_PATHS

to those directory and re-configure.

It is probably wisest 
to remove the configuration cache:

/bin/rm config.cache

before reconfiguring.
.
AC_MSG_ERROR(Exiting)
])

AC_DEFUN(BBH_MSG_NEED_MPI,[
cat<<.

Complete installation of DAGH and related software requires the 
MPI header file 'mpi.h' and the MPI library 'libmpi.a'; at least one 
of which could not be found in any of the usual locations.  
If you know the directories in which the headers and library are installed 
on this system, set the environment variables

LIB_PATHS
INCLUDE_PATHS

to include those directories and re-configure.

Alternatively, you can 

setenv DAGH_NO_MPI on

and re-configure to build without MPI.  It is probably wisest 
to remove the configuration cache:

/bin/rm config.cache

before reconfiguring.
.
AC_MSG_ERROR(Exiting)
])
dnl ----------------------------------------------------------------------
dnl Locate lib, include paths for auxiliary software needed by DAGH 
dnl and not currently needed by RNPL (including librnpl.a).  
dnl BBH_CHECK_RNPL_LIBS should always be invoked before this macro
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_CHECK_DAGH_WITH_RNPL_LIBS,[
   BBH_CHECK_FATAL="no"
   if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
      BBH_CHECK_LIB_PATHS="$LIB_PATHS $BBH_CHECK_LIB_PATHS /usr/local/lib"
      BBH_CHECK_INC_PATHS="$INCLUDE_PATHS $BBH_CHECK_INC_PATHS /usr/local/include"
   else
      BBH_CHECK_LIB_PATHS="$LIB_PATHS $BBH_CHECK_LIB_PATHS"
      BBH_CHECK_INC_PATHS="$INCLUDE_PATHS $BBH_CHECK_INC_PATHS"
   fi
dnl ----------------------------------------------------------------------
dnl   Look for the MPI headers and library ... Fatal exit can be 
dnl   overwritten by setting the environment variable DAGH_NO_MPI to some 
dnl   non-trivial value.
dnl ----------------------------------------------------------------------
   if test "X$DAGH_NO_MPI" = "X"; then
dnl ----------------------------------------------------------------------
dnl   Building with MPI
dnl ----------------------------------------------------------------------
      case "$BBH_SYSTEM/$BBH_SUBSYSTEM" in
dnl ----------------------------------------------------------------------
dnl   Don't need MPI on SP2
dnl ----------------------------------------------------------------------
      AIX/SP2) 
      ;;
      *) 
         AC_MSG_CHECKING(for mpi headers)
dnl         BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS /usr/local/mpi/include"
         if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
            BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS /pds/MPI/mpich/include"
         fi
         BBH_FINDFILE(BBH_CHECK_INC_PATHS,"mpi.h",BBH_INC_PATH)
         if test $BBH_INC_PATH != "no"; then
            AC_MSG_RESULT($BBH_INC_PATH/mpi.h)
            BBH_PUSHUNIQUE(BBH_RNPLAPP_INCPATHS,$BBH_INC_PATH)
dnl ----------------------------------------------------------------------
dnl         Guess whether we should be compiling with -DMPICH
dnl ----------------------------------------------------------------------
            MPICH=`echo $BBH_INC_PATH | grep mpich`
            AC_MSG_CHECKING(if we're using mpich)
            if test "X$MPICH" != "X"; then
               AC_MSG_RESULT(yes ... setting -DMPICH for C++)
               CXXFLAGS="$CXXFLAGS -DMPICH"
            else
               AC_MSG_RESULT(no)
            fi
         else
                AC_CHECK_HEADER(mpi.h,
                  HAVE_MPIH=1
                  )
                  if test "X$HAVE_MPIH" = "X"; then
                     if test "X$DAGH_NO_MPI" = "X"; then
                         AC_MSG_RESULT(not found)
                         AC_MSG_WARN(Can't find mpi headers ... see message below)
                         BBH_CHECK_FATAL="yes"
                     else
                         AC_MSG_RESULT(not found ... building without MPI)
                     fi
                  else
                     MPICH=`echo $BBH_INC_PATH | grep mpich`
                     AC_MSG_CHECKING(if we're using mpich)
                     if test "X$MPICH" != "X"; then
                         AC_MSG_RESULT(yes ... setting -DMPICH for C++)
                         CXXFLAGS="$CXXFLAGS -DMPICH"
                     else
                         AC_MSG_RESULT(no)
                     fi
                  fi
         fi

         AC_MSG_CHECKING(for mpi library)
dnl      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /usr/local/mpi/lib/IRIX64/ch_p4"
dnl      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /usr/local/mpi/lib/IRIXN32/ch_p4"
         BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /pds/MPI/mpich/lib/rs6000/ch_p4"
         BBH_FINDLIB(BBH_CHECK_LIB_PATHS,"mpi",BBH_LIB_PATH)
         if test $BBH_LIB_PATH != "no"; then
            AC_MSG_RESULT($BBH_LIB_PATH/libmpi.a)
            BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)
dnl ----------------------------------------------------------------------
dnl         MPI requires some OS-dependent additional libaries
dnl ----------------------------------------------------------------------
            case $BBH_SYSTEM in
            IRIX|IRIX64|IRIX32|IRIXN32) BBH_MPI_XTRA_LIBS=" "               ;; 
            SUNOS)                      BBH_MPI_XTRA_LIBS=" "               ;;
            AIX)                        BBH_MPI_XTRA_LIBS="-lbsd"           ;;
            A/UX)                       BBH_MPI_XTRA_LIBS=" "               ;; 
            OSF1)                       BBH_MPI_XTRA_LIBS=" "               ;;
            UNICOS)                     BBH_MPI_XTRA_LIBS=" "               ;;
            *)                          BBH_MPI_XTRA_LIBS=" "               ;;
            esac
            BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -lmpi $BBH_MPI_XTRA_LIBS"
            BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -lmpi $BBH_MPI_XTRA_LIBS"
            BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -lmpi $BBH_MPI_XTRA_LIBS"
            BBH_DEFS="$BBH_DEFS -DHAVE_LIBMPI=1"
         else
                AC_CHECK_LIB(mpi, MPI_Init,
                HAVE_LIBMPI=1
                )
                if test "X$HAVE_LIBMPI" = "X"; then
                     AC_MSG_RESULT(not found)
                     BBH_CHECK_FATAL="yes"
                else
            case $BBH_SYSTEM in
            IRIX|IRIX64|IRIX32|IRIXN32) BBH_MPI_XTRA_LIBS=" "               ;; 
            SUNOS)                      BBH_MPI_XTRA_LIBS=" "               ;;
            AIX)                        BBH_MPI_XTRA_LIBS="-lbsd"           ;;
            A/UX)                       BBH_MPI_XTRA_LIBS=" "               ;; 
            OSF1)                       BBH_MPI_XTRA_LIBS=" "               ;;
            UNICOS)                     BBH_MPI_XTRA_LIBS=" "               ;;
            *)                          BBH_MPI_XTRA_LIBS=" "               ;;
            esac
            BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -lmpi $BBH_MPI_XTRA_LIBS"
            BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -lmpi $BBH_MPI_XTRA_LIBS"
            BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -lmpi $BBH_MPI_XTRA_LIBS"
            BBH_DEFS="$BBH_DEFS -DHAVE_LIBMPI=1"
                fi
         fi

         if test $BBH_CHECK_FATAL = "yes"; then
            BBH_MSG_NEED_MPI
         fi
      ;;
      esac
   else
dnl ----------------------------------------------------------------------
dnl   Building without MPI
dnl ----------------------------------------------------------------------
      AC_MSG_WARN("Building without MPI")
      BBH_DEFS="$BBH_DEFS -DDAGH_NO_MPI -DACE_NO_MPI"
   fi

dnl ----------------------------------------------------------------------
dnl   Look for location of RNPL (bbhutil) headers and library (rnpl).
dnl ----------------------------------------------------------------------
   if test "X$NORNPL" = "X"; then
      AC_MSG_CHECKING(for bbhutil headers)
      if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
         BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS /usr/local/include /usr/local/rnpl/include"
      fi
      BBH_FINDFILE(BBH_CHECK_INC_PATHS,"bbhutil.h",BBH_INC_PATH)
      if test $BBH_INC_PATH != "no"; then
         CPPFLAGS="$CPPFLAGS -DIO_RNPLIO"
         BBH_PUSHUNIQUE(BBH_RNPLAPP_INCPATHS,$BBH_INC_PATH)
         AC_MSG_RESULT($BBH_INC_PATH/bbhutil.h)
      else
         AC_MSG_RESULT(not found ... some test apps may not build)
      fi
 
      AC_MSG_CHECKING(for rnpl library)
      BBH_FINDLIB(BBH_CHECK_LIB_PATHS,"rnpl",BBH_LIB_PATH)
      if test $BBH_LIB_PATH != "no"; then
         BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)
dnl ----------------------------------------------------------------------
dnl      Don't have to add '-lrnpl' since is done by BBH_CHECK_RNPL_LIBS 
dnl ----------------------------------------------------------------------
         AC_MSG_RESULT($BBH_LIB_PATH/librnpl.a)
      else
         AC_MSG_RESULT(not found ... some test apps may not build)
      fi
   else
      AC_MSG_WARN("Building without RNPL ... some test apps may not build")
   fi
])

AC_DEFUN(BBH_MSG_NEED_MPI,[
cat<<.

Complete installation on DAGH requires the MPI header file 'mpi.h' and 
the MPI library 'libmpi.a'; at least one of which could not be found in 
any of the usual locations.  If you know the directories in which the 
headers and library are installed on this system, set the 
environment variables

LIB_PATHS
INCLUDE_PATHS

to those directory and re-configure.

Alternatively, you can 

setenv DAGH_NO_MPI on

and re-configure to build DAGH without MPI.  It is probably wisest 
to remove the configuration cache:

/bin/rm config.cache

before reconfiguring.
.
AC_MSG_ERROR(Exiting)
])


dnl ----------------------------------------------------------------------
dnl Locate lib, include paths for auxiliary software needed by RNPL
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_CHECK_RNPL_LIBS,[
   BBH_CHECK_FATAL="no"
   BBH_RNPLAPP_LIBPATHS="$BBH_RNPLAPP_LIBPATHS $prefix/lib"
   BBH_RNPLAPP_INCPATHS="$BBH_RNPLAPP_INCPATHS $prefix/include"

   BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -lrnpl $BBH_MISC_FLIBS"
   BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -lrnpl "
   BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -lrnpl $LIBS $BBH_MISC_FLIBS"

   if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
      BBH_CHECK_INC_PATHS="$INCLUDE_PATHS /usr/include"
      BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS $HOME/include $HOME/hdf/include"
      BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS /usr/local/include /local/include"
      BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS /usr/local/hdf/include"
      BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS /usr/local/lib/hdf/include"
      BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS $prefix/include"
      BBH_CHECK_INC_PATHS="$BBH_CHECK_INC_PATHS $prefix/hdf/include"

      BBH_CHECK_LIB_PATHS="/usr/lib64 $LIB_PATHS"
      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS $HOME/lib $HOME/hdf/lib"
      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /usr/local/lib /local/lib" 
      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /usr/local/hdf/lib"
      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /usr/local/lib/hdf/lib"
      BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS $prefix/lib $prefix/hdf/lib"
   else
      BBH_CHECK_LIB_PATHS="$LIB_PATHS $BBH_CHECK_LIB_PATHS"
      BBH_CHECK_INC_PATHS="$INCLUDE_PATHS $BBH_CHECK_INC_PATHS"
   fi

dnl ----------------------------------------------------------------------
dnl Handle HDF ... now only execute if environment variable HDF is set.
dnl ----------------------------------------------------------------------
   if test "X$HDF" != "X"; then
      AC_MSG_CHECKING(for netcdf headers)
      BBH_FINDFILE(BBH_CHECK_INC_PATHS,"mfhdf.h",BBH_INC_PATH)
      if test $BBH_INC_PATH = "no"; then
         AC_MSG_RESULT(not found)
         AC_MSG_WARN(Can't find netcdf headers ... see message below)
         BBH_CHECK_FATAL="yes"
      else 
         AC_MSG_RESULT($BBH_INC_PATH/mfhdf.h)
         BBH_PUSHUNIQUE(BBH_RNPLAPP_INCPATHS,$BBH_INC_PATH)
      fi
      AC_MSG_CHECKING(for hdf headers) 
      BBH_FINDFILE(BBH_CHECK_INC_PATHS,"df.h",BBH_INC_PATH)
      if test $BBH_INC_PATH = "no"; then
         AC_MSG_RESULT(not found)
         AC_MSG_WARN(Can't find hdf headers ... see message below)
         BBH_CHECK_FATAL="yes"
      else 
         BBH_PUSHUNIQUE(BBH_RNPLAPP_INCPATHS,$BBH_INC_PATH)
         AC_MSG_RESULT($BBH_INC_PATH/df.h)
      fi
   dnl ----------------------------------------------------------------------
   dnl Name of 'netcdf' library changes with HDF 4.0 (sigh)
   dnl ----------------------------------------------------------------------
      BBH_APPEND_JPEG_Z="no"
      AC_MSG_CHECKING(for mfhdf library)
      BBH_FINDLIB(BBH_CHECK_LIB_PATHS,"mfhdf",BBH_LIB_PATH)
      if test $BBH_LIB_PATH = "no"; then
         AC_MSG_RESULT(not found ... will look for netcdf)
         BBH_FINDLIB(BBH_CHECK_LIB_PATHS,"netcdf",BBH_LIB_PATH)
         AC_MSG_CHECKING(for netcdf library)
         if test $BBH_LIB_PATH = "no"; then
            AC_MSG_RESULT(not found)
            AC_MSG_WARN(Can't find netcdf library ... see message below)
            BBH_CHECK_FATAL="yes"
         else
            AC_MSG_RESULT($BBH_LIB_PATH/libnetcdf.a)
            BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)

            BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -lnetcdf"
            BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -lnetcdf"
            BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -lnetcdf"
    
            BBH_HDFAPP_FLIBS="$BBH_HDFAPP_FLIBS -lnetcdf"
            BBH_HDFAPP_CLIBS="$BBH_HDFAPP_FLIBS -lnetcdf"

            BBH_DEFS="$BBH_DEFS -DHAVE_LIBNETCDF=1"
         fi 
      else 
         AC_MSG_RESULT($BBH_LIB_PATH/libmfhdf.a)
         BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)
         BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -lmfhdf"
         BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -lmfhdf"
         BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -lmfhdf"

         BBH_HDFAPP_FLIBS="$BBH_HDFAPP_FLIBS -lmfhdf"
         BBH_HDFAPP_CLIBS="$BBH_HDFAPP_FLIBS -lmfhdf"

         BBH_DEFS="$BBH_DEFS -DHAVE_LIBNETCDF=1"
         BBH_APPEND_JPEG_Z="yes"
      fi
      BBH_FINDLIB(BBH_CHECK_LIB_PATHS,"df",BBH_LIB_PATH)
      AC_MSG_CHECKING(for hdf library)
      if test $BBH_LIB_PATH = "no"; then
         AC_MSG_RESULT(not found)
         AC_MSG_WARN(Can't find hdf library ... see message below)
         BBH_CHECK_FATAL="yes"
      else 
         AC_MSG_RESULT($BBH_LIB_PATH/libdf.a)
         BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)
   dnl ----------------------------------------------------------------------
   dnl   Kludge for mario: ensure that local copy of 'libdf.a' gets linked
   dnl ----------------------------------------------------------------------
         if test "X$BBH_SYSTEM_MARIO" = "Xyes"; then
            BBH_MARIO_LIBDF="-l/usr/users/8/choptuik/lib/libdf.a"
            BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS $BBH_MARIO_LIBDF"
            BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS $BBH_MARIO_LIBDF"
            BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS $BBH_MARIO_LIBDF"

            BBH_HDFAPP_FLIBS="$BBH_HDFAPP_FLIBS $BBH_MARIO_LIBDF"
            BBH_HDFAPP_CLIBS="$BBH_HDFAPP_FLIBS $BBH_MARIO_LIBDF"
         else
            BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -ldf"
            BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -ldf"
            BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -ldf"

            BBH_HDFAPP_FLIBS="$BBH_HDFAPP_FLIBS -ldf"
            BBH_HDFAPP_CLIBS="$BBH_HDFAPP_FLIBS -ldf"
         fi

         BBH_DEFS="$BBH_DEFS -DHAVE_LIBDF=1"
      fi
      if test $BBH_APPEND_JPEG_Z = "yes"; then
   dnl ----------------------------------------------------------------------
   dnl   Assume that 'libjpeg.a' and 'libz.a' have been installed in same
   dnl   location as 'libdf.a' and 'libmfhdf.a'
   dnl ----------------------------------------------------------------------
         BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -ljpeg -lz"
         BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -ljpeg -lz"
         BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -ljpeg -lz"

         BBH_HDFAPP_FLIBS="$BBH_HDFAPP_FLIBS -ljpeg -lz"
         BBH_HDFAPP_CLIBS="$BBH_HDFAPP_FLIBS -ljpeg -lz"
      fi
   dnl ----------------------------------------------------------------------
   dnl   Define extra executables (for etc/Makefile.in) 
   dnl ----------------------------------------------------------------------
      BBH_RNPL_HDF_UTILITIES="hdfinfo hdftovs"
      BBH_DEFS="$BBH_DEFS -DHAVE_HDF=1"
   fi
dnl ----------------------------------------------------------------------
dnl   'vs' library.  There are now at least four versions of this 
dnl   library.  Clearly time has come to rationalize but will search
dnl   for vsso, then vs, but allow user to override with VS environment
dnl   vbl.
dnl ----------------------------------------------------------------------
   VS=${VS-"vsso vs"}
   for LIB in $VS; do
      AC_MSG_CHECKING(for $LIB library)
      BBH_FINDLIB(BBH_CHECK_LIB_PATHS,$LIB,BBH_LIB_PATH)
      if test $BBH_LIB_PATH != "no"; then
         AC_MSG_RESULT($BBH_LIB_PATH/lib$LIB.a)
         BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)
         BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -l$LIB"
         BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -l$LIB"
         BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -l$LIB"
         BBH_DEFS="$BBH_DEFS -DHAVE_LIBVS=1"
      else 
         AC_MSG_RESULT(not found)
      fi
   done
dnl ----------------------------------------------------------------------
dnl   'sv' library.
dnl ----------------------------------------------------------------------
   LIB="sv"
   AC_MSG_CHECKING(for $LIB library)
   BBH_FINDLIB(BBH_CHECK_LIB_PATHS,$LIB,BBH_LIB_PATH)
   if test $BBH_LIB_PATH != "no"; then
         AC_MSG_RESULT($BBH_LIB_PATH/lib$LIB.a)
         BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)
         BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -l$LIB"
         BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -l$LIB"
         BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -l$LIB"
         BBH_DEFS="$BBH_DEFS -DHAVE_LIBSV=1"
   else
      AC_MSG_RESULT(not found)
   fi
dnl ----------------------------------------------------------------------
dnl   'flex' library.
dnl ----------------------------------------------------------------------
   case BBH_SYSTEM in
   NEED_FL) 
      LIB="fl"
      AC_MSG_CHECKING(for flex library)
      BBH_FINDLIB(BBH_CHECK_LIB_PATHS,$LIB,BBH_LIB_PATH)
      if test $BBH_LIB_PATH != "no"; then
            AC_MSG_RESULT($BBH_LIB_PATH/lib$LIB.a)
            BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$BBH_LIB_PATH)
            BBH_RNPLAPP_FLIBS="$BBH_RNPLAPP_FLIBS -l$LIB"
            BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -l$LIB"
            BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -l$LIB"
            BBH_DEFS="$BBH_DEFS -DHAVE_LIBFLEX=1"
      else
         AC_MSG_RESULT(not found)
      fi
   ;;
   *) ;;
   esac
dnl ----------------------------------------------------------------------
dnl   Miscellaneous C libraries ...
dnl ----------------------------------------------------------------------
   BBH_RNPLAPP_CLIBS="$BBH_RNPLAPP_CLIBS -lm"
   BBH_RNPLBLD_CLIBS="$BBH_RNPLBLD_CLIBS -lm"
dnl ----------------------------------------------------------------------
dnl   Did we find everything?
dnl ----------------------------------------------------------------------
   if test $BBH_CHECK_FATAL = "yes"; then
      BBH_MSG_NEED_HDF
   fi
])

AC_DEFUN(BBH_MSG_NEED_HDF,[
cat<<.

Installation of RNPL requires HDF and NetCDF software which could 
not be found in any of the usual locations.  If this software 
exists on your system, set the following environment variables

INCLUDE_PATHS
LIB_PATHS

to appropriate values and re-configure.

.
AC_MSG_ERROR(Exiting)
])

dnl ----------------------------------------------------------------------
dnl BBH_RNPL_SET_VARS defines various environment variables used 
dnl in RNPL application Makefiles
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_RNPL_SET_VARS,[
   for el in $BBH_RNPLAPP_LIBPATHS; do 
      BBH_PUSHUNIQUE(BBH_RNPLAPP_LPATHS,"-L$el")
   done
   for el in $BBH_RNPLAPP_INCPATHS; do 
      BBH_PUSHUNIQUE(BBH_RNPLAPP_CINC,"-I$el")
      BBH_PUSHUNIQUE(BBH_RNPLAPP_FINC,"${F77INCSTEM}${el}")
   done

   AC_SUBST(BBH_HDFAPP_FLIBS)
   AC_SUBST(BBH_HDFAPP_CLIBS)

   AC_SUBST(BBH_MISC_FLFLAGS)
   AC_SUBST(BBH_MISC_FLIBS)

   AC_SUBST(BBH_RNPLAPP_CINC)
   AC_SUBST(BBH_RNPLAPP_FINC)
   AC_SUBST(BBH_RNPLAPP_DEFS)
   AC_SUBST(BBH_RNPLAPP_LPATHS)
   AC_SUBST(BBH_RNPLAPP_FLIBS)
   AC_SUBST(BBH_RNPLAPP_CLIBS)
   AC_SUBST(BBH_RNPLBLD_CLIBS)
   AC_SUBST(BBH_RNPL_FLAGS)
   RNPL=$prefix/bin/rnpl
   AC_SUBST(RNPL)

   AC_SUBST(BBH_RNPL_HDF_UTILITIES)
])

dnl ----------------------------------------------------------------------
dnl BBH_F77_CONFIGURE invokes various macros to set up f77 environment
dnl for rnpl.  Supply an argument to abort if f77transcray not found
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_F77_CONFIGURE,[
   BBH_PROG_F77
   if test $F77 != "None"; then
      BBH_F77_TRANSLATE($1)
      BBH_F77_T3E_FIX
      BBH_F77_T90_FIX
      BBH_PROG_F77FLAGS
      BBH_PROG_F77INCSTEM
      BBH_F77_HAS_SYSTEM
      BBH_F77_HAS_MKDIR
      BBH_F77_HAS_CHDIR
   else
      AC_MSG_WARN(No Fortran 77 compiler found on this system)
   fi
])

dnl ----------------------------------------------------------------------
dnl Supply an argument to abort if f77transcray not found
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_F77_TRANSLATE,[
   if test "X$BBH_SYSTEM" = "X"; then
      BBH_SYS_GETSYSTEM
   fi
   if test "X$F77_TRANSFORM" = "X"; then
      case $BBH_SYSTEM in
      UNICOS)
         AC_MSG_CHECKING(for f77transcray);
         BBH_CHECK_PROG_NOMESS(F77_TRANSFORM,f77transcray,f77transcray,"None")
         if test $F77_TRANSFORM = "None"; then
            if test "X$1" != "X"; then
               AC_MSG_RESULT(not found ... aborting configuration)
               exit
            else
               AC_MSG_RESULT(not found)
            fi
         else
            AC_MSG_RESULT(f77transcray)
         fi
      ;;
      *)    F77_TRANSFORM="touch";;
      esac;
   fi
   AC_SUBST(F77_TRANSFORM)
])

dnl ----------------------------------------------------------------------
dnl Supply an argument to abort if f77transcray not found
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_F77_TRANSLATE_OLD,[
   if test "X$BBH_SYSTEM" = "X"; then
      BBH_SYS_GETSYSTEM
   fi
   case $BBH_SYSTEM in
   UNICOS)
      AC_MSG_CHECKING(for f77transcray);
      BBH_CHECK_PROG_NOMESS(F77_TRANSFORM,f77transcray,f77transcray,"None")
      if test $F77_TRANSFORM = "None"; then
         if test "X$1" != "X"; then
            AC_MSG_RESULT(not found ... aborting configuration)
            exit
         else
            AC_MSG_RESULT(not found)
         fi
      else
         AC_MSG_RESULT(f77transcray)
      fi
   ;;
   *)    F77_TRANSFORM="touch";;
   esac;
   AC_SUBST(F77_TRANSFORM)
])

dnl ----------------------------------------------------------------------
dnl BBH_F90_CONFIGURE invokes various macros to set up f90 environment
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_F90_CONFIGURE,[
   BBH_PROG_F90
   if test $F90 != "None"; then
dnl   BBH_F77_T3E_FIX
      BBH_PROG_F90FLAGS
dnl   BBH_PROG_F77INCSTEM
dnl   BBH_F77_HAS_SYSTEM
   else
      AC_MSG_WARN(No Fortran 90 compiler found on this system)
   fi
])

dnl ----------------------------------------------------------------------
dnl Checks for libt3e_utill.a on T3Es (makes up for getenv(), getarg()
dnl deficiency.
dnl
dnl Interpolate
dnl
dnl    @BBH_MISC_FLFLAGS@
dnl    @BBH_MISC_FLIBS@
dnl
dnl in Makefile.in
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_F77_T3E_FIX,[
   if test "+$BBH_SUBSYSTEM" = "+T3E"; then

      L_INCLUDE_PATHS="$INCLUDE_PATHS $prefix/include /usr/local/include /local/include $HOME/include $HOMEMWC/include"
      L_LIB_PATHS="$LIB_PATHS $prefix/lib /usr/local/lib /local/lib $HOME/lib $HOMEMWC/lib"

      AC_MSG_CHECKING(for t3e_util library)
      BBH_FINDLIB(L_LIB_PATHS,"t3e_util",L_LIB)
      if test $L_LIB = "no"; then
         AC_MSG_RESULT(not found)
         L_FATAL="yes"
      else
         AC_MSG_RESULT($L_LIB/libt3e_util.a)
         BBH_PUSHUNIQUE(BBH_MISC_FLFLAGS,-L$L_LIB)
         BBH_PUSHUNIQUE(BBH_MISC_FLIBS,-lt3e_util)
      fi

   fi
   AC_SUBST(BBH_MISC_FLFLAGS)
   AC_SUBST(BBH_MISC_FLIBS)
])


dnl ----------------------------------------------------------------------
dnl Checks for libt3e_utill.a on T90s (makes up for getenv(), getarg()
dnl deficiency.
dnl
dnl Interpolate
dnl
dnl    @BBH_MISC_FLFLAGS@
dnl    @BBH_MISC_FLIBS@
dnl
dnl in Makefile.in
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_F77_T90_FIX,[
   if test "+$BBH_SUBSYSTEM" = "+T90"; then

      L_INCLUDE_PATHS="$INCLUDE_PATHS $prefix/include /usr/local/include /local/include $HOME/include $HOMEMWC/include"
      L_LIB_PATHS="$LIB_PATHS $prefix/lib /usr/local/lib /local/lib $HOME/lib $HOMEMWC/lib"

      AC_MSG_CHECKING(for t3e_util library (on t90))
      BBH_FINDLIB(L_LIB_PATHS,"t3e_util",L_LIB)
      if test $L_LIB = "no"; then
         AC_MSG_RESULT(not found)
         L_FATAL="yes"
      else
         AC_MSG_RESULT($L_LIB/libt3e_util.a)
         BBH_PUSHUNIQUE(BBH_MISC_FLFLAGS,-L$L_LIB)
         BBH_PUSHUNIQUE(BBH_RNPLAPP_LIBPATHS,$L_LIB)
         BBH_PUSHUNIQUE(BBH_MISC_FLIBS,-lt3e_util)
      fi

   fi
   AC_SUBST(BBH_MISC_FLFLAGS)
   AC_SUBST(BBH_MISC_FLIBS)
])

dnl ----------------------------------------------------------------------
dnl BBH_PROG_F77 attempts to determine what the Fortran 77 compiler 
dnl is called 
dnl
dnl Sets Output Variable: F77
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_PROG_F77,[
   BBH_SYS_GETSYSTEM
   AC_MSG_CHECKING(for Fortran 77 compiler)
   AC_CACHE_VAL(ac_cv_bbh_prog_f77,[
      if test "X$F77" = "X"; then
         case $BBH_SYSTEM in
         IRIX|IRIX64|IRIX32|IRIXN32)   
                               F77='f77';                                ;;
         SUNOS)                F77='f77';                                ;; 
         AIX)                  
            case $BBH_SUBSYSTEM in
            SP2) F77='mpxlf';;
            *)   F77='xlf';;
            esac                                                         ;;
         A/UX)                 F77='f77';                                ;;
         OSF1)                 F77='f90';                                ;;
         UNICOS)               
            case $BBH_SUBSYSTEM in
            T3E*) F77='f90';;
            T90)  F77='f90'; echo 'BBH_PROG_F77_setting F77 to f90';;
            SV1)  F77='f90';;
            *)    F77='cf77';;
            esac;
                                                                         ;;
         *)                    F77='f77';                                ;;
         esac
         BBH_CHECK_PROG_NOMESS(ac_cv_bbh_prog_f77,$F77,$F77,"None")
      else
         ac_cv_bbh_prog_f77=$F77
      fi
      if test $ac_cv_bbh_prog_f77 != "None"; then
         cat <<END > _foo.f
         write(*,*) 'hello world!'
         stop
         end
END
         if $F77 _foo.f -o _foo > /dev/null 2>&1; then
            ac_cv_bbh_prog_f77=$F77
            rm _foo.f > /dev/null 2>&1 
            rm _foo.o > /dev/null 2>&1 
            rm _foo > /dev/null 2>&1 
         else
            ac_cv_bbh_prog_f77="None"
            rm _foo.f > /dev/null 2>&1
            rm _foo.o > /dev/null 2>&1 
            rm _foo > /dev/null 2>&1
         fi
      fi
   ])
   AC_MSG_RESULT($ac_cv_bbh_prog_f77)
   F77=$ac_cv_bbh_prog_f77
   AC_SUBST(F77)
])

dnl ----------------------------------------------------------------------
dnl BBH_PROG_F90 attempts to determine what the Fortran 77 compiler 
dnl is called 
dnl
dnl Sets Output Variable: F90
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_PROG_F90,[
   BBH_SYS_GETSYSTEM
   AC_MSG_CHECKING(for Fortran 90 compiler)
   AC_CACHE_VAL(ac_cv_bbh_prog_f90,[
      if test "X$F90" = "X"; then
         case $BBH_SYSTEM in
         IRIX|IRIX64|IRIX32|IRIXN32)   
                               F90='f90';                                ;;
         SUNOS)                F90='f90';                                ;; 
         AIX)                  
            case $BBH_SUBSYSTEM in
            SP2) F90='mpxlf';;
            *)   F90='xlf';;
            esac                                                         ;;
         A/UX)                 F90='f90';                                ;;
         OSF1)                 F90='f90';                                ;;
         UNICOS)               F90='f90';                                ;;
         *)                    F90='f90';                                ;;
         esac
         BBH_CHECK_PROG_NOMESS(ac_cv_bbh_prog_f90,$F90,$F90,"None")
      else
         ac_cv_bbh_prog_f90=$F90
      fi
      if test $ac_cv_bbh_prog_f90 != "None"; then
         cat <<END > _foo.F
         write(*,*) 'hello world!'
         stop
         end
END
         if $F90 _foo.F -o _foo > /dev/null 2>&1; then
            ac_cv_bbh_prog_f90=$F90
            rm _foo.F > /dev/null 2>&1 
            rm _foo.o > /dev/null 2>&1 
            rm _foo > /dev/null 2>&1 
         else
            ac_cv_bbh_prog_f90="None"
            rm _foo.f > /dev/null 2>&1
            rm _foo.o > /dev/null 2>&1 
            rm _foo > /dev/null 2>&1
         fi
      fi
   ])
   AC_MSG_RESULT($ac_cv_bbh_prog_f90)
   F90=$ac_cv_bbh_prog_f90
   AC_SUBST(F90)
])


dnl ----------------------------------------------------------------------
dnl BBH_CHECK_PROG_NOMESS: Version of AC_CHECK_PROG which disables 
dnl messages ...
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_CHECK_PROG_NOMESS,
[# Extract the first word of "$2", so it can be a program name with args.
set dummy $2; ac_word=[$]2
AC_CACHE_VAL(ac_cv_prog_$1,
[if test -n "[$]$1"; then
  ac_cv_prog_$1="[$]$1" # Let the user override the test.
else
  IFS="${IFS=   }"; ac_save_ifs="$IFS"; IFS="${IFS}:"
  for ac_dir in $PATH; do
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/$ac_word; then
      ac_cv_prog_$1="$3"
      break
    fi
  done
  IFS="$ac_save_ifs"
dnl If no 4th arg is given, leave the cache variable unset,
dnl so AC_CHECK_PROGS will keep looking.
ifelse([$4], , , [  test -z "[$]ac_cv_prog_$1" && ac_cv_prog_$1="$4"
])dnl
fi])dnl
$1="$ac_cv_prog_$1"
AC_SUBST($1)dnl
])

dnl ----------------------------------------------------------------------
dnl BBH_PROG_F77FLAGS attempts to determine what  Fortran 77 flags 
dnl should be used
dnl
dnl Sets Output Variables: F77FLAGS, F77R8FLAG, F77FLAGSNOOPT
dnl
dnl setenv F77FLAGS to override default flags 
dnl setenv F77OPT to override default optimization level
dnl
dnl F77R8FLAG is system-specific flag for making real*8 default for 
dnl real computations.
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_PROG_F77FLAGS,[
   BBH_SYS_GETSYSTEM
   BBH_GET_NPROC
   AC_MSG_CHECKING(for Fortran 77 flags (with optimization))
   AC_CACHE_VAL(ac_cv_bbh_prog_f77flags,[
      if test "X$F77FLAGS" = "X"; then
         case $BBH_SYSTEM in
         LINUX)
            F77FLAGSNOOPT=" ";
            F77OPT=${F77OPT-"-O"};
                                                                    ;;
         LINUX_PG)
            F77FLAGSNOOPT=" ";
            F77OPT=${F77OPT-"-O"};
                                                                    ;;
         DARWIN)
            F77FLAGSNOOPT=" ";
            F77OPT=${F77OPT-"-O"};
                                                                    ;;
         IRIX)       
            F77FLAGSNOOPT=" ";
            F77OPT=${F77OPT-"-O"};
                                                                    ;;
         IRIX64)     
            case "X$BBH_NPROC" in
            X1) F77FLAGSNOOPT="-64 -NC200";;
            *)  F77FLAGSNOOPT="-64 -WK,-listoption=O -pfa -mp -NC200";;
            esac;
            F77OPT=${F77OPT-"-O3"};
            test "X${BBH_CHECK_DEFAULTS}" != "XNONE" && LDFLAGS="$LDFLAGS -L/usr/local/lib";
            
                                                                    ;;
         IRIXN32)    
            F77FLAGSNOOPT="-n32";                         
            F77OPT=${F77OPT-"-O"};
            test "X${BBH_CHECK_DEFAULTS}" != "XNONE" && LDFLAGS="$LDFLAGS -L/usr/localn32/lib";
                                                                    ;;
         IRIX32)     
            F77FLAGSNOOPT="-32";                    
            F77OPT=${F77OPT-"-O"};
            test "X${BBH_CHECK_DEFAULTS}" != "XNONE" && LDFLAGS="$LDFLAGS -L/usr/local/lib32";
                                                                    ;;
         SUNOS)      
            F77FLAGSNOOPT=" ";
            F77OPT=${F77OPT-"-O"};
                                                                    ;; 
         AIX)        
            F77FLAGSNOOPT="-qextname";                    
            F77OPT=${F77OPT-"-O"};
                                                                    ;;
         A/UX)       
            F77FLAGSNOOPT=" ";
            F77OPT=${F77OPT-"-O"};
                                                                    ;;
         OSF1)       
dnl         F77FLAGSNOOPT="-fpe2 -check nopower -check nounderflow";
dnl         F77OPT=${F77OPT-"-O4"};
            F77FLAGSNOOPT="";
            F77OPT=${F77OPT-"-fast"};
                                                                    ;;
         UNICOS)     
            case $BBH_SUBSYSTEM in
            T3E*) F77FLAGSNOOPT=" ";
                  F77OPT=${F77OPT-"-O3"};
                                                                 ;;
            T90*) F77FLAGSNOOPT=" ";
                  F77OPT=${F77OPT-"-O vector3 -O task3 -O inline3"};
                                                                 ;;
            SV1*) F77FLAGSNOOPT=" ";
                  F77OPT=${F77OPT-"-O vector3 -O task3 -O inline4"};
                                                                 ;;
            *)    F77FLAGSNOOPT=" ";       
                  F77OPT=${F77OPT-"-O3"};
                                                                 ;;
            esac; 
                                                                    ;;
         *)          
            F77FLAGSNOOPT=" ";                            
                                                                    ;;
         esac
         ac_cv_bbh_prog_f77flags="$F77FLAGSNOOPT $F77OPT"
      else
         ac_cv_bbh_prog_f77flags="$F77FLAGS"
      fi
   ])
   AC_MSG_RESULT($ac_cv_bbh_prog_f77flags)
   AC_MSG_CHECKING(for Fortran 77 flags (without optimization))
   AC_MSG_RESULT($F77FLAGSNOOPT)
   F77FLAGS=$ac_cv_bbh_prog_f77flags
   AC_SUBST(F77FLAGS)
   AC_SUBST(F77FLAGSNOPT)

   AC_MSG_CHECKING(for Fortran 77 real*8 default flag)
   AC_CACHE_VAL(ac_cv_bbh_prog_f77r8flag,[
      if test "X$F77R8FLAG" = "X"; then
         case $BBH_SYSTEM in
         IRIX*)       
            ac_cv_bbh_prog_f77r8flag=${F77R8FLAG-"-r8"};                 
                                                                    ;;
         UNICOS)     
            ac_cv_bbh_prog_f77r8flag=${F77R8FLAG-"-dp"};            
                                                                    ;;
         *)          
            ac_cv_bbh_prog_f77r8flag=' ';                           
                                                                    ;;
         esac
      else
         ac_cv_bbh_prog_f77r8flag="$F77R8FLAG"
      fi
   ])
   AC_MSG_RESULT($ac_cv_bbh_prog_f77r8flag)
   F77R8FLAG=$ac_cv_bbh_prog_f77r8flag
   AC_SUBST(F77R8FLAG)
])

dnl ----------------------------------------------------------------------
dnl BBH_PROG_F90FLAGS attempts to determine what  Fortran 90 flags 
dnl should be used
dnl
dnl Sets Output Variables: F90FLAGS, F90R8FLAG, F90FLAGSNOOPT
dnl
dnl setenv F90FLAGS to override default flags 
dnl setenv F90OPT to override default optimization level
dnl
dnl F90R8FLAG is system-specific flag for making real*8 default for 
dnl real computations.
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_PROG_F90FLAGS,[
   BBH_SYS_GETSYSTEM
   BBH_GET_NPROC
   AC_MSG_CHECKING(for Fortran 90 flags (with optimization))
   AC_CACHE_VAL(ac_cv_bbh_prog_f90flags,[
      if test "X$F90FLAGS" = "X"; then
         case $BBH_SYSTEM in
         IRIX)       
            F90FLAGSNOOPT="-freeform";
            F90OPT=${F90OPT-"-O"};
                                                                    ;;
         IRIX64)     
dnl Default here designed to handle ADM code and like
            _OPT="-OPT:const_copy_limit=35000 -OPT:fold_arith_limit=30000"
            _OPT="$_OPT -OPT:global_limit=35000 -OPT:fprop_limit=1221"
            case "X$BBH_NPROC" in
            X1) F90FLAGSNOOPT="-64 -freeform $_OPT";;
            *)  F90FLAGSNOOPT="-64 -freeform $_OPT -WK,-listoption=O -pfa -mp";;
            esac;
            F90OPT=${F90OPT-"-O3"};
            test "X${BBH_CHECK_DEFAULTS}" != "XNONE" && LDFLAGS="$LDFLAGS -L/usr/local/lib";
                                                                    ;;
         IRIXN32)    
            _OPT="-OPT:const_copy_limit=35000 -OPT:fold_arith_limit=30000";
            _OPT="$_OPT -OPT:global_limit=35000 -OPT:fprop_limit=1221";
            F90FLAGSNOOPT="-n32 -freeform $_OPT";
            F90OPT=${F90OPT-"-O3"};
            test "X${BBH_CHECK_DEFAULTS}" != "XNONE" && LDFLAGS="$LDFLAGS -L/usr/localn32/lib";
                                                                    ;;
         IRIX32)     
            F90FLAGSNOOPT="-32";                            
            test "X${BBH_CHECK_DEFAULTS}" != "XNONE" && LDFLAGS="$LDFLAGS -L/usr/local/lib32";
            F90OPT=${F90OPT-"-O3"};
                                                                    ;;
         SUNOS)      
            F90OPT=${F90OPT-"-O"};
            F90FLAGSNOOPT=" ";       
                                                                    ;; 
         AIX)        
            F90FLAGSNOOPT="-qextname -qfree=f90";
            F90OPT=${F90OPT-"-O"};
                                                                    ;;
         A/UX)       
            F90FLAGSNOOPT=" ";        
                                                                    ;;
         OSF1)       
            F90FLAGSNOOPT="-fpe2 -check nopower -check nounderflow";
            F90OPT=${F90OPT-"-O4 -fpe2 -check nopower -check nounderflow"};
                                                                    ;;
         UNICOS)     
            case $BBH_SUBSYSTEM in
            T3E*) F90FLAGSNOOPT="-f free";
                  F90OPT=${F90OPT-"-O3"};
                                                                 ;;
            T90*) F90FLAGSNOOPT="-f free";         
                  F90OPT=${F90OPT-"-O vector3 -O task3 -O inline3"};
                                                                 ;;
            SV1*) F90FLAGSNOOPT="-f free";         
                  F90OPT=${F90OPT-"-O vector3 -O task3 -O inline3"};
                                                                 ;;
            *)    F90FLAGSNOOPT=" ";       
                                                                 ;;
            esac; 
                                                                    ;;
         *)          
            F90FLAGSNOOPT=" ";                                   ;;
         esac
         ac_cv_bbh_prog_f90flags="$F90FLAGSNOOPT $F90OPT"
      else
         ac_cv_bbh_prog_f90flags="$F90FLAGS";
      fi
   ])
   AC_MSG_RESULT($ac_cv_bbh_prog_f90flags)
   AC_MSG_CHECKING(for Fortran 90 flags (without optimization))
   AC_MSG_RESULT($F90FLAGSNOOPT)
   F90FLAGS=$ac_cv_bbh_prog_f90flags
   AC_SUBST(F90FLAGS)
   AC_SUBST(F90FLAGSNOOPT)

   AC_MSG_CHECKING(for Fortran 90 real*8 default flag)
   AC_CACHE_VAL(ac_cv_bbh_prog_f90r8flag,[
      if test "X$F90R8FLAG" = "X"; then
         case $BBH_SYSTEM in
         IRIX*)       
            ac_cv_bbh_prog_f90r8flag=${F90R8FLAG-"-r8"};                 
                                                                    ;;
         UNICOS)     
            ac_cv_bbh_prog_f90r8flag=${F90R8FLAG-"-dp"};            
                                                                    ;;
         *)          
            ac_cv_bbh_prog_f90r8flag=' ';                           
                                                                    ;;
         esac
      else
         ac_cv_bbh_prog_f90r8flag="$F90R8FLAG"
      fi
   ])
   AC_MSG_RESULT($ac_cv_bbh_prog_f90r8flag)
   F90R8FLAG=$ac_cv_bbh_prog_f90r8flag
   AC_SUBST(F90R8FLAG)
])

dnl ----------------------------------------------------------------------
dnl BBH_PROG_F77INCSTEM sets stem for 'include' files 
dnl
dnl Sets Output Variable: F77INCSTEM
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_PROG_F77INCSTEM,[
   BBH_SYS_GETSYSTEM
   AC_MSG_CHECKING(for Fortran 77 include stem)
   AC_CACHE_VAL(ac_cv_bbh_prog_f77incstem,[
      case $BBH_SYSTEM in
      IRIX64)     ac_cv_bbh_prog_f77incstem='-YI,'                       ;;
      IRIXN32)    ac_cv_bbh_prog_f77incstem='-YI,'                       ;;
      IRIX32)     ac_cv_bbh_prog_f77incstem='-Wf,-I'                     ;;
      IRIX)       ac_cv_bbh_prog_f77incstem='-Wf,-I'                     ;;
      SUNOS)      ac_cv_bbh_prog_f77incstem='-I';                        ;; 
      AIX)        ac_cv_bbh_prog_f77incstem='-I';                        ;;
      A/UX)       ac_cv_bbh_prog_f77incstem='-I';                        ;;
      OSF1)       ac_cv_bbh_prog_f77incstem='-I';                        ;;
      UNICOS)     ac_cv_bbh_prog_f77incstem='-I';                        ;;
      LINUX*)     ac_cv_bbh_prog_f77incstem='-I';                        ;;
      DARWIN*)    ac_cv_bbh_prog_f77incstem='-I';                        ;;
      *)          ac_cv_bbh_prog_f77incstem=' ';                         ;;
      esac
   ])
   AC_MSG_RESULT($ac_cv_bbh_prog_f77incstem)
   F77INCSTEM=$ac_cv_bbh_prog_f77incstem
   AC_SUBST(F77INCSTEM)
])

dnl ----------------------------------------------------------------------
dnl BBH_LIB_HDF attempts to locate the HDF and NETCDF include files 
dnl and libraries 
dnl Sets Output Variables: BBH_HDF_INCLUDEDIR, BBH_HDF_LIBDIR
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_LIB_HDF,[
])

dnl ----------------------------------------------------------------------
dnl BBH_F77_HAS_SYSTEM checks to see whether Fortran
dnl accepts 'call system()'
dnl Sets Output Variables: F77_HAS_SYSTEM (iff 'call system()' works)
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_F77_HAS_SYSTEM,[
   AC_MSG_CHECKING(whether Fortran 77 supports 'system()')
   AC_CACHE_VAL(ac_cv_bbh_f77_has_system,[
      if test $F77 != None; then
         cat <<END > _foo.f
            call system('touch _foo.f')
            stop
            end
END
         if $F77 $F77FLAGS $LDFLAGS _foo.f -o _foo > _output 2>&1; then
dnl         This hack necessary since UNICOS only generates warning for 
dnl         unsatisfied external references
            if grep -i Unsatisfied _output > /dev/null; then
               ac_cv_bbh_f77_has_system="no"
            else
               ac_cv_bbh_f77_has_system="yes"
            fi
         else
            ac_cv_bbh_f77_has_system="no"
         fi
         rm _output > /dev/null 2>&1
      fi
      rm _foo.f > /dev/null 2>&1 
      rm _foo.o > /dev/null 2>&1 
      rm _foo > /dev/null 2>&1 
   ])
   AC_MSG_RESULT($ac_cv_bbh_f77_has_system)
   if test $ac_cv_bbh_f77_has_system = "yes"; then
      BBH_DEFS="$BBH_DEFS -DF77_HAS_SYSTEM"
   fi
])

dnl ----------------------------------------------------------------------
dnl BBH_F77_HAS_MKDIR checks to see whether Fortran
dnl accepts 'rc=mkdir()'
dnl Sets Output Variables: F77_HAS_MKDIR (iff 'rc=mkdir()' works)
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_F77_HAS_MKDIR,[
   AC_MSG_CHECKING(whether Fortran 77 supports 'mkdir()')
   AC_CACHE_VAL(ac_cv_bbh_f77_has_mkdir,[
      if test $F77 != None; then
         cat <<END > _foo.f
                integer rc,mkdir
            rc=mkdir('foodir',511)
            stop
            end
END
         if $F77 $F77FLAGS $LDFLAGS _foo.f -o _foo > _output 2>&1; then
dnl         This hack necessary since UNICOS only generates warning for 
dnl         unsatisfied external references
            if grep -i Unsatisfied _output > /dev/null; then
               ac_cv_bbh_f77_has_mkdir="no"
            else
               ac_cv_bbh_f77_has_mkdir="yes"
            fi
         else
            ac_cv_bbh_f77_has_mkdir="no"
         fi
         rm _output > /dev/null 2>&1
      fi
      rm _foo.f > /dev/null 2>&1 
      rm _foo.o > /dev/null 2>&1 
      rm _foo > /dev/null 2>&1 
   ])
   AC_MSG_RESULT($ac_cv_bbh_f77_has_mkdir)
   if test $ac_cv_bbh_f77_has_mkdir = "yes"; then
      BBH_DEFS="$BBH_DEFS -DF77_HAS_MKDIR"
   fi
])

dnl ----------------------------------------------------------------------
dnl BBH_F77_HAS_CHDIR checks to see whether Fortran
dnl accepts 'rc=chdir()'
dnl Sets Output Variables: F77_HAS_CHDIR (iff 'rc=chdir()' works)
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_F77_HAS_CHDIR,[
   AC_MSG_CHECKING(whether Fortran 77 supports 'chdir()')
   AC_CACHE_VAL(ac_cv_bbh_f77_has_chdir,[
      if test $F77 != None; then
         cat <<END > _foo.f
                integer rc,chdir
            rc=chdir('foodir')
            stop
            end
END
         if $F77 $F77FLAGS $LDFLAGS _foo.f -o _foo > _output 2>&1; then
dnl         This hack necessary since UNICOS only generates warning for 
dnl         unsatisfied external references
            if grep -i Unsatisfied _output > /dev/null; then
               ac_cv_bbh_f77_has_chdir="no"
            else
               ac_cv_bbh_f77_has_chdir="yes"
            fi
         else
            ac_cv_bbh_f77_has_chdir="no"
         fi
         rm _output > /dev/null 2>&1
      fi
      rm _foo.f > /dev/null 2>&1 
      rm _foo.o > /dev/null 2>&1 
      rm _foo > /dev/null 2>&1 
   ])
   AC_MSG_RESULT($ac_cv_bbh_f77_has_chdir)
   if test $ac_cv_bbh_f77_has_chdir = "yes"; then
      BBH_DEFS="$BBH_DEFS -DF77_HAS_CHDIR"
   fi
])

AC_DEFUN(BBH_F77_CHECK_LIB,[
   l_lib=$1; l_routine=$2;
   AC_MSG_CHECKING([for library -l$1])
   AC_CACHE_VAL(ac_cv_bbh_f77_lib_$l_lib,[
      if test $F77 != None; then
         cat <<END > _foo.f
            call $l_routine()
            stop
            end
END
         if $F77 $F77FLAGS $LDFLAGS _foo.f -l$l_lib -o _foo > _output 2>&1; then
dnl         This hack necessary since UNICOS only generates warning for 
dnl         unsatisfied external references
            if grep -i Unsatisfied _output > /dev/null; then
               eval "ac_cv_bbh_f77_lib_$l_lib=no"
               AC_MSG_RESULT("no");
            else
               eval "ac_cv_bbh_f77_lib_$l_lib=yes"
               LIBS="$LIBS -l$l_lib"
               AC_MSG_RESULT("yes");
            fi
         else
               eval "ac_cv_bbh_f77_lib_$l_lib=no"
               AC_MSG_RESULT("no");
         fi
         rm _output > /dev/null 2>&1
      fi
      rm _foo.f > /dev/null 2>&1 
      rm _foo.o > /dev/null 2>&1 
      rm _foo > /dev/null 2>&1 
   ])
])

dnl ----------------------------------------------------------------------
dnl BBH_RNPL_DEFAULT_PATH sets CPP variable BBH_RNPL_DEFAULT_PATH
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_RNPL_DEFAULT_PATH,[
   BBH_DEFS="$BBH_DEFS -DBBH_RNPL_DEFAULT_PATH='\"$prefix/lib/rnpl\"'"
])

dnl ----------------------------------------------------------------------
dnl Locate lib, include paths for auxiliary software needed by 
dnl Klasky's interpolator
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_CHECK_INTERP_LIBS,[
   BBH_CHECK_FATAL="no"
   BBH_INTERP_LIBS="-lbbhinterp"

   case $BBH_SYSTEM in
   UNICOS) 
      AC_MSG_WARN(Unicos system: Assuming 'linpack', 'blas' in 'libsci.a');
   ;;
   *)
      if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
         BBH_CHECK_LIB_PATHS="/usr/lib64 $LIB_PATHS"
         BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS $HOME/lib $HOME/linpack/lib"
         BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /usr/lib" 
         BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS /usr/local/lib /local/lib" 
         BBH_CHECK_LIB_PATHS="$BBH_CHECK_LIB_PATHS $prefix/lib $prefix/linpack/lib"
      fi

      AC_MSG_CHECKING(for linpack library)
      BBH_FINDLIB(BBH_CHECK_LIB_PATHS,"linpack",BBH_LIB_PATH)
      if test $BBH_LIB_PATH = "no"; then
         AC_MSG_RESULT(not found)
         AC_MSG_WARN(Can't find linpack library ... see message below)
         BBH_CHECK_FATAL="yes"
      else 
         AC_MSG_RESULT($BBH_LIB_PATH/liblinpack.a)
         BBH_PUSHUNIQUE(BBH_INTERP_LIBPATHS,$BBH_LIB_PATH)
         BBH_INTERP_LIBS="$BBH_INTERP_LIBS -llinpack"
      fi

      AC_MSG_CHECKING(for blas library)
      BBH_FINDSO(BBH_CHECK_LIB_PATHS,"blas",BBH_LIB_PATH)
      if test $BBH_LIB_PATH = "no"; then 
         BBH_FINDLIB(BBH_CHECK_LIB_PATHS,"blas",BBH_LIB_PATH)
         if test $BBH_LIB_PATH = "no"; then
            AC_MSG_RESULT(not found)
            AC_MSG_WARN(Can't find blas library ... see message below)
            BBH_CHECK_FATAL="yes"
         else 
            AC_MSG_RESULT($BBH_LIB_PATH/libblas.a)
            BBH_PUSHUNIQUE(BBH_INTERP_LIBPATHS,$BBH_LIB_PATH)
            BBH_INTERP_LIBS="$BBH_INTERP_LIBS -lblas"
         fi
      else
         AC_MSG_RESULT($BBH_LIB_PATH/libblas.so)
         BBH_PUSHUNIQUE(BBH_INTERP_LIBPATHS,$BBH_LIB_PATH)
         BBH_INTERP_LIBS="$BBH_INTERP_LIBS -lblas"
      fi
   ;;
   esac

   if test $BBH_CHECK_FATAL = "yes"; then
      BBH_MSG_NEED_LINPACK_BLAS
   fi
])

dnl ----------------------------------------------------------------------
dnl BBH_INTERP_SET_VARS defines various environment variables used
dnl by interpolation software
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_INTERP_SET_VARS,[
   for el in $BBH_INTERP_LIBPATHS; do
      BBH_PUSHUNIQUE(BBH_INTERP_LPATHS,"-L$el")
   done
   case $BBH_SYSTEM in
   AIX)        F77FLAGS='-O';                           ;;
   esac
   AC_SUBST(F77FLAGS)
   AC_SUBST(BBH_INTERP_LPATHS)
   AC_SUBST(BBH_INTERP_LIBS)
])

AC_DEFUN(BBH_MSG_NEED_LINPACK_BLAS,[
cat<<.

Installation of 'bbhinterp' requires LINPACK (liblinpack.a)
and BLAS (libblas.a).  These could not be found in any of the 
usual locations.  If you know that this software exists on your 
system, set the following environment variable

LIB_PATHS

to an appropriate value and re-configure.

.
AC_MSG_ERROR(Exiting)
])

dnl ----------------------------------------------------------------------
dnl Locate lib, include paths for BBH utilities (libbbhutil.a) and 
dnl interpolate variables 
dnl
dnl Include:
dnl 
dnl   BBH_BBHUTIL_SETUP
dnl
dnl in configure.in, and use
dnl
dnl   @BBH_UTIL_LDFLAGS@
dnl   @BBH_UTIL_CCFLAGS@
dnl   @BBH_UTIL_LIBS@
dnl
dnl in Makefile.in files.
dnl
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_BBHUTIL_SETUP,[
   L_FATAL="no"

   BBH_UTIL_CCFLAGS=""
   BBH_UTIL_LDFLAGS=""
   BBH_UTIL_LIBS=""

   if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
      L_INCLUDE_PATHS="$INCLUDE_PATHS /usr/local/include /local/include"
      L_LIB_PATHS="$LIB_PATHS /usr/local/lib /local/lib"
   else
      L_INCLUDE_PATHS="$INCLUDE_PATHS"
      L_LIB_PATHS="$LIB_PATHS"
   fi

   AC_MSG_CHECKING(for bbhutil headers)
   BBH_FINDFILE(L_INCLUDE_PATHS,"bbhutil.h",L_INC)
   if test $L_INC = "no"; then
      AC_MSG_RESULT(not found)
      L_FATAL="yes"
   else
      AC_MSG_RESULT($L_INC/bbhutil.h)
      BBH_PUSHUNIQUE(BBH_UTIL_CCFLAGS,-I$L_INC)
   fi
   AC_SUBST(BBH_UTIL_CCFLAGS)

   AC_MSG_CHECKING(for bbhutil library)
   BBH_FINDLIB(L_LIB_PATHS,"bbhutil",L_LIB)
   if test $L_LIB = "no"; then
      AC_MSG_RESULT(not found)
      L_FATAL="yes"
   else
      AC_MSG_RESULT($L_LIB/libbhutil.a)
      BBH_PUSHUNIQUE(BBH_UTIL_LDFLAGS,-L$L_LIB)
      BBH_PUSHUNIQUE(BBH_UTIL_LIBS,-lbbhutil)
   fi
   AC_SUBST(BBH_UTIL_LDFLAGS)
   AC_SUBST(BBH_UTIL_LIBS)

])

dnl ----------------------------------------------------------------------
dnl Locate lib, include paths for BBH io utilities (libbbhutil.a) and 
dnl interpolate variables 
dnl
dnl Include:
dnl 
dnl   BBH_BBIO_SETUP
dnl
dnl in configure.in, and use
dnl
dnl   @BBH_IO_LDFLAGS@
dnl   @BBH_IO_CCFLAGS@
dnl   @BBH_IO_LIBS@
dnl
dnl in Makefile.in files.
dnl
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_BBHIO_SETUP,[
   L_FATAL="no"

   BBH_IO_CCFLAGS=""
   BBH_IO_LDFLAGS=""
   BBH_IO_LIBS=""

   if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
      L_INCLUDE_PATHS="$INCLUDE_PATHS /usr/local/include /local/include"
      L_LIB_PATHS="$LIB_PATHS /usr/local/lib /local/lib"
   else
      L_INCLUDE_PATHS="$INCLUDE_PATHS"
      L_LIB_PATHS="$LIB_PATHS"
   fi

   AC_MSG_CHECKING(for bbhio headers)
   BBH_FINDFILE(L_INCLUDE_PATHS,"bbhio.h",L_INC)
   if test $L_INC = "no"; then
      AC_MSG_RESULT(not found)
      L_FATAL="yes"
   else
      AC_MSG_RESULT($L_INC/bbhio.h)
      BBH_PUSHUNIQUE(BBH_IO_CCFLAGS,-I$L_INC)
   fi
   AC_SUBST(BBH_IO_CCFLAGS)

   AC_MSG_CHECKING(for bbhio library)
   BBH_FINDLIB(L_LIB_PATHS,"bbhio",L_LIB)
   if test $L_LIB = "no"; then
      AC_MSG_RESULT(not found)
      L_FATAL="yes"
   else
      AC_MSG_RESULT($L_LIB/libbhio.a)
      BBH_PUSHUNIQUE(BBH_IO_LDFLAGS,-L$L_LIB)
      BBH_PUSHUNIQUE(BBH_IO_LIBS,-lbbhio)
   fi
   AC_SUBST(BBH_IO_LDFLAGS)
   AC_SUBST(BBH_IO_LIBS)

])

dnl ----------------------------------------------------------------------
dnl Locate Choptuik's f77 support libs 
dnl
dnl   BBH_MWC_F77_SETUP
dnl
dnl in configure.in, and use
dnl
dnl   @BBH_MWC_F77_LDFLAGS@
dnl   @BBH_MWC_F77_LIBS@
dnl   @BBH_MWC_F77_FLIBS@
dnl
dnl in Makefile.in files.
dnl
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_MWC_F77_SETUP,[
   L_FATAL="no"

   BBH_MWC_F77_LIBS=""

   if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
      L_LIB_PATHS="$LIB_PATHS /usr/local/lib $HOME/lib $HOMEMWC/lib"
   else
      L_LIB_PATHS="$LIB_PATHS $HOME/lib $HOMEMWC/lib"
   fi

   AC_MSG_CHECKING("for libvs.a")
   BBH_FINDLIB(L_LIB_PATHS,vs,L_LIB)
   if test $L_LIB = "no"; then
      AC_MSG_RESULT(not found)
   else
      AC_MSG_RESULT($L_LIB/libvs.a)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_LIBS,-lvs)
   fi
   AC_MSG_CHECKING("for libfvs.a")
   BBH_FINDLIB(L_LIB_PATHS,fvs,L_LIB)
   if test $L_LIB = "no"; then
      AC_MSG_RESULT(not found)
   else
      AC_MSG_RESULT($L_LIB/libfvs.a)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_FLIBS,-lfvs)
   fi
   AC_MSG_CHECKING("for libsvs.a")
   BBH_FINDLIB(L_LIB_PATHS,svs,L_LIB)
   if test $L_LIB = "no"; then
      AC_MSG_RESULT(not found)
   else
      AC_MSG_RESULT($L_LIB/libsvs.a)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_SLIBS,-lsvs)
      AC_MSG_CHECKING("for libbbhutil.a")
      BBH_FINDLIB(L_LIB_PATHS,bbhutil,L_LIB)
      if test $L_LIB = "no"; then
         AC_MSG_RESULT(not found)
      else
         AC_MSG_RESULT($L_LIB/libbbhutil.a)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_SLIBS,-lbbhutil)
      fi
   fi
   AC_MSG_CHECKING("for libxvs.a")
   BBH_FINDLIB(L_LIB_PATHS,xvs,L_LIB)
   if test $L_LIB = "no"; then
      AC_MSG_RESULT(not found)
   else
      AC_MSG_RESULT($L_LIB/libxvs.a)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_XLIBS,-lxvs)
      AC_MSG_CHECKING("for libbbhutil.a")
      BBH_FINDLIB(L_LIB_PATHS,bbhutil,L_LIB)
      if test $L_LIB = "no"; then
         AC_MSG_RESULT(not found)
      else
         AC_MSG_RESULT($L_LIB/libbbhutil.a)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_XLIBS,-lbbhutil)
      fi
   fi
   AC_MSG_CHECKING("for libjvs.a")
   BBH_FINDLIB(L_LIB_PATHS,jvs,L_LIB)
   if test $L_LIB = "no"; then
      AC_MSG_RESULT(not found)
   else
      AC_MSG_RESULT($L_LIB/libjvs.a)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_JLIBS,-ljvs)
   fi
   AC_MSG_CHECKING("for libsv.a")
   BBH_FINDLIB(L_LIB_PATHS,sv,L_LIB)
   if test $L_LIB = "no"; then
      AC_MSG_RESULT(not found)
   else
      AC_MSG_RESULT($L_LIB/libsv.a)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_JLIBS,-lsv)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_XLIBS,-lsv)
      BBH_APPEND_UNIQUE(BBH_MWC_F77_SLIBS,-lsv)
   fi
   for l in f77sig vutil utilmath utilio; do
      AC_MSG_CHECKING("for lib$l.a")
      BBH_FINDLIB(L_LIB_PATHS,$l,L_LIB)
      if test $L_LIB = "no"; then
         AC_MSG_RESULT(not found)
      else
         AC_MSG_RESULT($L_LIB/lib$l.a)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_LDFLAGS,-L$L_LIB)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_LIBS,-l$l)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_FLIBS,-l$l)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_JLIBS,-l$l)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_XLIBS,-l$l)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_XLIBS,-l$l)
         BBH_APPEND_UNIQUE(BBH_MWC_F77_SLIBS,-l$l)
      fi
   done
   AC_SUBST(BBH_MWC_F77_LDFLAGS)
   AC_SUBST(BBH_MWC_F77_LIBS)
   AC_SUBST(BBH_MWC_F77_FLIBS)
   AC_SUBST(BBH_MWC_F77_JLIBS)
   AC_SUBST(BBH_MWC_F77_XLIBS)
   AC_SUBST(BBH_MWC_F77_SLIBS)
])


dnl ----------------------------------------------------------------------
dnl Some utility macros
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_PUSHUNIQUE,[
   _pu_push="yes"
   for _pu_el in ${$1}; do
      if test $_pu_el = $2; then
         _pu_push="no"
      fi
   done
   if test $_pu_push = "yes"; then
      $1="${$1} "" $2"
   fi
])

AC_DEFUN(BBH_PREPEND_UNIQUE,[
   _pu_pre="yes"
   for _pu_el in ${$1}; do
      if test $_pu_el = $2; then
         _pu_pre="no"
      fi
   done
   if test $_pu_pre = "yes"; then
      $1="${$1} "" $2"
   fi
])

AC_DEFUN(BBH_APPEND_UNIQUE,[
   _pu_app="yes"
   for _pu_el in ${$1}; do
      if test $_pu_el = $2; then
         _pu_app="no"
      fi
   done
   if test $_pu_app = "yes"; then
      $1="${$1} "" $2"
   fi
])

AC_DEFUN(BBH_FINDLIB,[
   $3="no"
   findlib_done=""
   for _stem in ${$1}; do
      if test -f $_stem/lib$2.a; then
         $3=$_stem
         findlib_done="yes"
         break
      fi
   done
   if test "X$findlib_done" = "X"; then
      for _stem in ${$1}; do
         if test -f $_stem/lib$2.so; then
            $3=$_stem
            break
         fi
      done
   fi
])

AC_DEFUN(BBH_FINDSO,[
   $3="no"
   for _stem in ${$1}; do
      if test -f $_stem/lib$2.so; then
         $3=$_stem
         break
      fi
   done
])

AC_DEFUN(BBH_FINDFILE,[
   $3="no"
   for _stem in ${$1}; do
      if test -f $_stem/$2; then
         $3=$_stem
			BBH_FILE=$_stem/$2
         break
      fi
   done
])

AC_DEFUN(BBH_SETDEF,[
   if test "X${$1}" = "X"; then
      $1=$2
   fi
])

dnl ----------------------------------------------------------------------
dnl Usage: 
dnl    BBH_CHECK_LIBS(libs,paths,LDPATHS,LIBS,ABORT?)
dnl 
dnl Example:
dnl    BBH_CHECK_LIBS("vutil utilmath","$HOME/lib",MYLDPATHS,MYLIBS,yes)
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_CHECK_LIBS,[
   BBH_CHECK_FATAL="no"
   BBH_BAD_ARGS="no"
   if test "X$3" = "X"; then
      echo "[BBH_CHECK_LIBS]: Error: Third argument is missing"
      BBH_BAD_ARGS="yes"
   fi
   if test "X$4" = "X"; then
      echo "[BBH_CHECK_LIBS]: Error: Fourth argument is missing"
      BBH_BAD_ARGS="yes"
   fi
   if test $BBH_BAD_ARGS = "no"; then
      if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
         _PATHS="$LIB_PATHS "$2" /usr/lib64 /usr/lib64/mips4 $_PATHS /usr/local/lib $HOME/lib"
      else
         _PATHS="$LIB_PATHS "$2" $_PATHS"
      fi
      for lib in "$1"; do
         AC_MSG_CHECKING(for $lib library)
         BBH_FINDLIB(_PATHS,$lib,_PATH)
         if test "$_PATH" != "no"; then
            AC_MSG_RESULT($_PATH/lib$lib.a)
dnl Don't add -L/usr/lib to load flags
            if test $_PATH != "/usr/lib"; then
               BBH_PUSHUNIQUE($3,-L$_PATH)
            fi
            BBH_APPEND_UNIQUE($4,-l$lib)
         else
            if test "X$5" != "X"; then
               AC_MSG_RESULT(not found ... will have to abort configuration)
cat<<END
setenv LIB_PATHS if you know the location of the library
END
               BBH_CHECK_FATAL="yes"
            else
               AC_MSG_RESULT(not found ... will ignore)
            fi
         fi
      done
      AC_SUBST($3)
      AC_SUBST($4)
   fi
   if test $BBH_CHECK_FATAL = "yes"; then
      exit
   fi
])

dnl ----------------------------------------------------------------------
dnl Usage: 
dnl    BBH_CHECK_HEADERS(includes,paths,INCPATHS,ABORT?)
dnl 
dnl Example:
dnl    BBH_CHECK_HEADERS("mpi.h","/usr/local/mpi/include",MYINCPATHS,yes)
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_CHECK_HEADERS,[
   BBH_CHECK_FATAL="no"
   BBH_BAD_ARGS="no"
   if test "X$3" = "X"; then
      echo "[BBH_CHECK_HEADERS]: Error: Third argument is missing"
      BBH_BAD_ARGS="yes"
   fi
   if test $BBH_BAD_ARGS = "no"; then
      if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
         _PATHS="$INCLUDE_PATHS "$2" $_PATHS /usr/local/include $HOME/include"
      else
         _PATHS="$INCLUDE_PATHS "$2" $_PATHS"
      fi
      for include in "$1"; do
         AC_MSG_CHECKING(for $include)
         BBH_FINDFILE(_PATHS,$include,_PATH)
         if test "$_PATH" != "no"; then
            AC_MSG_RESULT($_PATH/$include)
            if test "$_PATH" != "/usr/include"; then
               BBH_PUSHUNIQUE($3,-I$_PATH)
            fi
         else
            if test "X$4" != "X"; then
               AC_MSG_RESULT(not found ... will have to abort configuration)
cat<<END
setenv INCLUDE_PATHS if you know the location of the headers
END
               BBH_CHECK_FATAL="yes"
            else
               AC_MSG_RESULT(not found ... will ignore)
            fi
         fi
      done
      AC_SUBST($3)
   fi
   if test $BBH_CHECK_FATAL = "yes"; then
      exit
   fi
])

dnl ----------------------------------------------------------------------
dnl Usage: 
dnl    BBH_LOCATE_HEADER(header,paths,var,ABORT?)
dnl 
dnl Example:
dnl    BBH_LOCATE_HEADER("mpi.h","/usr/local/mpi/include",MPIH,yes)
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_LOCATE_HEADER,[
   BBH_CHECK_FATAL="no"
   BBH_BAD_ARGS="no"
   if test "X$3" = "X"; then
      echo "[BBH_LOCATE_HEADERS]: Error: Third argument is missing"
      BBH_BAD_ARGS="yes"
   fi
   if test $BBH_BAD_ARGS = "no"; then
		include=$1
      if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
         _PATHS="$INCLUDE_PATHS "$2" $_PATHS /usr/local/include $HOME/include"
      else
         _PATHS="$INCLUDE_PATHS "$2" $_PATHS"
      fi
		AC_MSG_CHECKING(for location of $include)
		BBH_FINDFILE(_PATHS,$include,_PATH)
		if test "$_PATH" != "no"; then
			AC_MSG_RESULT($_PATH/$include)
			$3=$BBH_FILE
			AC_SUBST($3)
		else
			if test "X$4" != "X"; then
				AC_MSG_RESULT(not found ... will have to abort configuration)
cat<<END
setenv INCLUDE_PATHS if you know the location of the header file
END
				BBH_CHECK_FATAL="yes"
			else
				AC_MSG_RESULT(not found ... will ignore)
			fi
		fi
   fi
   if test $BBH_CHECK_FATAL = "yes"; then
      exit
   fi
])

dnl ----------------------------------------------------------------------
dnl Usage: 
dnl    BBH_CHECK_F77_HEADERS(includes,paths,INCPATHS,ABORT?)
dnl 
dnl Example:
dnl    BBH_CHECK_F77_HEADERS("mpi.h","/usr/local/mpi/include",MYINCPATHS,yes)
dnl
dnl Be sure to invoke BBH_F77_CONFIGURE before this macro.
dnl ----------------------------------------------------------------------

AC_DEFUN(BBH_CHECK_F77_HEADERS,[
   BBH_CHECK_FATAL="no"
   BBH_BAD_ARGS="no"
   if test "X$3" = "X"; then
      echo "[BBH_CHECK_F77_HEADERS]: Error: Third argument is missing"
      BBH_BAD_ARGS="yes"
   fi
   if test $BBH_BAD_ARGS = "no"; then
      if test "X${BBH_CHECK_DEFAULTS}" != "XNONE"; then
         _PATHS="$INCLUDE_PATHS "$2" $_PATHS /usr/local/include $HOME/include"
      else
         _PATHS="$INCLUDE_PATHS "$2" $_PATHS"
      fi
      for include in "$1"; do
         AC_MSG_CHECKING(for $include)
         BBH_FINDFILE(_PATHS,$include,_PATH)
         if test "$_PATH" != "no"; then
            AC_MSG_RESULT($_PATH/$include)
            BBH_PUSHUNIQUE($3,"$F77INCSTEM""$_PATH")
         else
            if test "X$4" != "X"; then
               AC_MSG_RESULT(not found ... will have to abort configuration)
cat<<END
setenv INCLUDE_PATHS if you know the location of the headers
END
               BBH_CHECK_FATAL="yes"
            else
               AC_MSG_RESULT(not found ... will ignore)
            fi
         fi
      done
      AC_SUBST($3)
   fi
   if test $BBH_CHECK_FATAL = "yes"; then
      exit
   fi
])

dnl ----------------------------------------------------------------------
dnl Installs appropriate routine for producing f77 external name.
dnl BBH_SYS_GETSYSTEM must be called before this macro.
dnl
dnl August 19, 2003: Removed LINUX_PG from LINUX/DARWIN equivalent.
dnl
dnl September 22, 2003: Made LINUX/DARWIN equivalent to SUN.
dnl ----------------------------------------------------------------------
AC_DEFUN(BBH_INSTALL_GENF77EXTERN,[
   AC_MSG_CHECKING(for f77 external name routine )
   if test "X$BBH_SYSTEM" != "X"; then
      case $BBH_SYSTEM in
      UNICOS) 
         if test -f src/genf77extern_unicos.c; then
            cp src/genf77extern_unicos.c src/genf77extern.c
            AC_MSG_RESULT(src/genf77extern_unicos.c -> src/genf77extern.c)
         else
            AC_MSG_RESULT(Warning: src/genf77extern_unicos.c not found)
         fi
      ;;
dnl   LINUX|DARWIN)
dnl      if test -f src/genf77extern_linux.c; then
dnl         cp src/genf77extern_linux.c src/genf77extern.c
dnl         AC_MSG_RESULT(src/genf77extern_linux.c -> src/genf77extern.c)
dnl      else
dnl         AC_MSG_RESULT(Warning: src/genf77extern_linux.c not found)
dnl      fi
dnl   ;;
      *)
         if test -f src/genf77extern_sun.c; then
            cp src/genf77extern_sun.c src/genf77extern.c
            AC_MSG_RESULT(src/genf77extern_sun.c -> src/genf77extern.c)
         else
            AC_MSG_RESULT(Warning: src/genf77extern_sun.c not found)
         fi
      ;;
      esac
   else
      AC_MSG_RESULT(BBH_SYSTEM not defined ... no routine installed)
   fi
])

