#!/usr/bin/perl

use strict;

my($i, $plot_path, $parallel, %paths, $cc, $CC, $make_daemon, $lapack_installed, $lapack_version, $blas_version, $plot_version, $cgi_installed, $boost_installed, @out, $myloc, $j);

$lapack_installed=0;
$boost_installed=0;
$cgi_installed=0;
$plot_path ="";
$parallel=0;
$make_daemon=0;
$plot_version="";

if(glob("/usr/lib/x86_64-linux-gnu/libcgicc*")) {
    $cgi_installed=1;
}

if(glob("/usr/lib/libcgicc*")) {
    $cgi_installed=1;
}

if(glob("/usr/lib/x86_64-linux-gnu/libboost_filesystem*")) {
    $boost_installed=1;
}

if(glob("/usr/lib/libboost_filesystem*")) {
    $boost_installed=1;
}

if (@ARGV > 0) {
	for($i=0; $i<@ARGV; $i++) {
		if ($ARGV[$i] =~ /^-p/i) {
			$plot_path =$ARGV[$i];
			$plot_path = substr($plot_path, 3, length($plot_path)-3);
            print "Build will include visualization tools: $plot_path\n";
            $myloc = $plot_path . "libplot*";
            #print "Checking for plot library: $myloc\n";
            @out=`ls $myloc`;
            for($j=0; $j<@out; $j++) {
                if ($out[$j] =~ /(libplot\.so\.\d)/) {
                    $plot_version=$1;
                }
            }
            
		}
		if ($ARGV[$i] =~ /omp/i) {
			$parallel=1;
            		print "Build will use OpenMP\n";
		}
		if ($ARGV[$i] =~ /browser/i) {
			$make_daemon=1;
			print "Will build executables for POInT browser\n";
		}
	}
}

if (($plot_path eq "") && ($make_daemon==1)) {
    print "ERROR: No plot library given to build plotting daemon\n";
    exit();
}

if (($cgi_installed == 0) && ($make_daemon==1)) {
    print "ERROR: No cgi library installed to build plotting daemon\n";
    exit();
}

if (($boost_installed == 0) && ($make_daemon==1)) {
    print "ERROR: No boost libraries installed to build plotting daemon\n";
    exit();
}

#if ( -e "Makefile") {`rm Makefile`;}
if (-e "libf2c/Makefile") {`rm libf2c/Makefile`}
if (-e "ranlib/Makefile") {`rm ranlib/Makefile`}
if (-e "lapack/Makefile") {`rm lapack/Makefile`}
if (-e "lapack/make.inc") {`rm lapack/make.inc`}
if (-e "src/Makefile") {`rm src/Makefile`}

$lapack_version="";
$blas_version="";

if((glob("/usr/lib/liblapack*")) && (glob("/usr/lib/libblas*"))) {
    $lapack_installed=1;
    @out=`ls /usr/lib/liblapack*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(liblapack\.so\.\d)/) {
            $lapack_version=$1;
        }
    }
    
    @out=`ls /usr/lib/libblas*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(libblas\.so\.\d)/) {
            $blas_version=$1;
        }
    }
}

if((glob("/usr/local/lib/liblapack*")) && (glob("/usr/local/lib/libblas*"))) {
    $lapack_installed=1;
    
    @out=`ls /usr/local/lib/liblapack*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(liblapack\.so\.\d)/) {
            $lapack_version=$1;
        }
    }
    
    @out=`ls /usr/local/lib/libblas*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(libblas\.so\.\d)/) {
            $blas_version=$1;
        }
    }
    
}

if((glob("/lib/liblapack*")) && (glob("/lib/libblas*"))) {
    $lapack_installed=1;
    
    @out=`ls /lib/liblapack*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(liblapack\.so\.\d)/) {
            $lapack_version=$1;
        }
    }
    
    @out=`ls /lib/libblas*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(libblas\.so\.\d)/) {
            $blas_version=$1;
        }
    }
}

if((glob("/usr/lib/x86_64-linux-gnu/liblapack*")) && (glob("/usr/lib/x86_64-linux-gnu/libblas*"))) {
    $lapack_installed=1;
    
    @out=`ls /usr/lib/x86_64-linux-gnu/liblapack*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(liblapack\.so\.\d)/) {
            $lapack_version=$1;
        }
    }
    
    @out=`ls /usr/lib/x86_64-linux-gnu/libblas*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(libblas\.so\.\d)/) {
            $blas_version=$1;
        }
    }
}

$CC=`which g++`;
$cc=`which gcc`;

$CC =~ s/\n//;

if ($CC eq "") {
    $CC=`which c++`;
    if ($CC eq "") {
        $CC=`which xlC`;
        if ($CC eq "") {
            $CC=`which Clang`;
                if($CC eq "") {
                    $CC = `which icpc`;
                
                    if ($CC eq "") {
                        print "ERROR: No c++ compiler found\n";
                        exit();
                    }
                }
        }
	}
}

print "Found $CC for c++\n";

$cc=`which gcc`;

$cc =~ s/\n//;

if ($cc eq "") {
    $cc=`which cc`;
    if ($cc eq "") {
        $cc=`which xlc`;
        if ($cc eq "") {
                $cc=`which clang`;
                        if($cc eq "") {
                                $cc = `which icc`;

                                if ($cc eq "") {
                                        print "ERROR: No c compiler found\n";
                                        exit();
                                }
                        }
        }

    }
}

print "Found $cc for c\n";

if($lapack_installed == 0) {
    open(WRITEMAKE, ">libf2c/Makefile") or die;
    print WRITEMAKE "# Unix makefile: see README.\n# For C++, first \"make hadd\".\n# If your compiler does not recognize ANSI C, add\n";
    print WRITEMAKE "#    -DKR_headers\n# to the CFLAGS = line below.\n# On Sun and other BSD systems that do not provide an ANSI sprintf, add\n#    -DUSE_STRLEN\n";
    print WRITEMAKE "# to the CFLAGS = line below.\n# On Linux systems, add\n#    -DNON_UNIX_STDIO\n# to the CLFAGS = line below.\n";
    print WRITEMAKE ".SUFFIXES: .c .o\n";
    print WRITEMAKE "CC = $cc\n";
    print WRITEMAKE "SHELL = /bin/sh\n";
    print WRITEMAKE "CFLAGS = -O\n";

    print WRITEMAKE "# compile, then strip unnecessary symbols\n.c.o:\n";
    print WRITEMAKE "\t\$(CC) -c -DSkip_f2c_Undefs \$(CFLAGS) \$*.c\n";
    print WRITEMAKE "\tld -r -x -o \$*.xxx \$*.o\n\tmv \$*.xxx \$*.o\n";
    print WRITEMAKE "## Under Solaris (and other systems that do not understand ld -x),\n## omit -x in the ld line above.\n## If your system does not have the ld command, comment out\n";
    print WRITEMAKE "## or remove both the ld and mv lines above.\n";

    print WRITEMAKE "MISC =    f77vers.o i77vers.o main.o s_rnge.o abort_.o exit_.o getarg_.o iargc_.o\\\ngetenv_.o signal_.o s_stop.o s_paus.o system_.o cabs.o\\\nderf_.o derfc_.o erf_.o erfc_.o sig_die.o uninit.o\n";
    print WRITEMAKE "POW =    pow_ci.o pow_dd.o pow_di.o pow_hh.o pow_ii.o pow_ri.o pow_zi.o pow_zz.o\n";
    print WRITEMAKE "CX =    c_abs.o c_cos.o c_div.o c_exp.o c_log.o c_sin.o c_sqrt.o\n";
    print WRITEMAKE "DCX =    z_abs.o z_cos.o z_div.o z_exp.o z_log.o z_sin.o z_sqrt.o\n";
    print WRITEMAKE "REAL =    r_abs.o r_acos.o r_asin.o r_atan.o r_atn2.o r_cnjg.o r_cos.o\\\nr_cosh.o r_dim.o r_exp.o r_imag.o r_int.o\\\nr_lg10.o r_log.o r_mod.o r_nint.o r_sign.o\\\n";
    print WRITEMAKE "r_sin.o r_sinh.o r_sqrt.o r_tan.o r_tanh.o\n";
    print WRITEMAKE "DBL =    d_abs.o d_acos.o d_asin.o d_atan.o d_atn2.o\\\nd_cnjg.o d_cos.o d_cosh.o d_dim.o d_exp.o\\\nd_imag.o d_int.o d_lg10.o d_log.o d_mod.o\\\nd_nint.o d_prod.o d_sign.o d_sin.o d_sinh.o\\\nd_sqrt.o d_tan.o d_tanh.o\n";
    print WRITEMAKE "INT =    i_abs.o i_dim.o i_dnnt.o i_indx.o i_len.o i_mod.o i_nint.o i_sign.o lbitbits.o lbitshft.o\n";
    print WRITEMAKE "HALF =    h_abs.o h_dim.o h_dnnt.o h_indx.o h_len.o h_mod.o h_nint.o h_sign.o\n";
    print WRITEMAKE "CMP =    l_ge.o l_gt.o l_le.o l_lt.o hl_ge.o hl_gt.o hl_le.o hl_lt.o\n";
    print WRITEMAKE "EFL =    ef1asc_.o ef1cmc_.o\n";
    print WRITEMAKE "CHAR =    f77_aloc.o s_cat.o s_cmp.o s_copy.o\n";
    print WRITEMAKE "I77 =    backspac.o close.o dfe.o dolio.o due.o endfile.o err.o\\\nfmt.o fmtlib.o ftell_.o iio.o ilnw.o inquire.o lread.o lwrite.o\\\nopen.o rdfmt.o rewind.o rsfe.o rsli.o rsne.o sfe.o sue.o\\\ntypesize.o uio.o util.o wref.o wrtfmt.o wsfe.o wsle.o wsne.o xwsne.o\n";
    print WRITEMAKE "QINT =    pow_qq.o qbitbits.o qbitshft.o\n";
    print WRITEMAKE "TIME =    dtime_.o etime_.o\n";

    print WRITEMAKE "# If you get an error compiling dtime_.c or etime_.c, try adding\n# -DUSE_CLOCK to the CFLAGS assignment above; if that does not work,\n# omit $(TIME) from the dependency list for libf2c.a below.\n";

    print WRITEMAKE "# For INTEGER*8 support (which requires system-dependent adjustments to\n# f2c.h), add $(QINT) to the libf2c.a dependency list below...\n";

    print WRITEMAKE "all: f2c.h signal1.h libf2c.a\n";

    print WRITEMAKE "libf2c.a: \$(MISC) \$(POW) \$(CX) \$(DCX) \$(REAL) \$(DBL) \$(INT) \\\n\$(HALF) \$(CMP) \$(EFL) \$(CHAR) \$(I77) \$(TIME)\n";
    print WRITEMAKE "\tar r libf2c.a \$?\n";
    print WRITEMAKE "\tranlib libf2c.a\n";

    print WRITEMAKE "### If your system lacks ranlib, you don't need it; see README.\n";

    print WRITEMAKE "f77vers.o: f77vers.c\n";
    print WRITEMAKE "\t\$(CC) -c f77vers.c\n";

    print WRITEMAKE "i77vers.o: i77vers.c\n\t\$(CC) -c i77vers.c\n";

    print WRITEMAKE "# To get an \"f2c.h\" for use with \"f2c -C++\", first \"make hadd\"\n";
    print WRITEMAKE "hadd: f2c.h0 f2ch.add\n\tcat f2c.h0 f2ch.add >f2c.h\n";

    print WRITEMAKE "# For use with \"f2c\" and \"f2c -A\":\n";
    print WRITEMAKE "f2c.h: f2c.h0\n\tcp f2c.h0 f2c.h\n";

    print WRITEMAKE "# You may need to adjust signal1.h suitably for your system...\n";
    print WRITEMAKE "signal1.h: signal1.h0\n\tcp signal1.h0 signal1.h\n";

    print WRITEMAKE "# If your system lacks onexit() and you are not using an\n# ANSI C compiler, then you should uncomment the following\n# two lines (for compiling main.o):\n#main.o: main.c\n";
    print WRITEMAKE "#    $(CC) -c -DNO_ONEXIT -DSkip_f2c_Undefs main.c\n# On at least some Sun systems, it is more appropriate to\n# uncomment the following two lines:\n#main.o: main.c\n#    $(CC) -c -Donexit=on_exit -DSkip_f2c_Undefs main.c\n";

    print WRITEMAKE "install: libf2c.a\n\tcp libf2c.a \$(LIBDIR)\n\tranlib \$(LIBDIR)/libf2c.a\n";

    print WRITEMAKE "clean:\n\trm -f libf2c.a *.o arith.h\n";

    print WRITEMAKE "backspac.o:    fio.h\nclose.o:    fio.h\ndfe.o:        fio.h\ndfe.o:        fmt.h\ndue.o:        fio.h\n";
    print WRITEMAKE "endfile.o:    fio.h rawio.h\nerr.o:        fio.h rawio.h\nfmt.o:        fio.h\nfmt.o:        fmt.h\n";
    print WRITEMAKE "iio.o:        fio.h\niio.o:        fmt.h\nilnw.o:        fio.h\nilnw.o:        lio.h\ninquire.o:    fio.h\n";
    print WRITEMAKE "lread.o:    fio.h\nlread.o:    fmt.h\nlread.o:    lio.h\nlread.o:    fp.h\nlwrite.o:    fio.h\nlwrite.o:    fmt.h\n";
    print WRITEMAKE "lwrite.o:    lio.h\nopen.o:        fio.h rawio.h\nrdfmt.o:    fio.h\nrdfmt.o:    fmt.h\nrdfmt.o:    fp.h\nrewind.o:    fio.h\n";
    print WRITEMAKE "rsfe.o:        fio.h\nrsfe.o:        fmt.h\nrsli.o:        fio.h\nrsli.o:        lio.h\nrsne.o:        fio.h\nrsne.o:        lio.h\n";
    print WRITEMAKE "sfe.o:        fio.h\nsue.o:        fio.h\nuio.o:        fio.h\nuninit.o:    arith.h\nutil.o:        fio.h\nwref.o:        fio.h\n";
    print WRITEMAKE "wref.o:        fmt.h\nwref.o:        fp.h\nwrtfmt.o:    fio.h\nwrtfmt.o:    fmt.h\nwsfe.o:        fio.h\nwsfe.o:        fmt.h\n";
    print WRITEMAKE "wsle.o:        fio.h\nwsle.o:        fmt.h\nwsle.o:        lio.h\nwsne.o:        fio.h\nwsne.o:        lio.h\nxwsne.o:    fio.h\n";
    print WRITEMAKE "xwsne.o:    lio.h\nxwsne.o:    fmt.h\n\n";

    print WRITEMAKE "arith.h: arithchk.c\n\t\$(CC) \$(CFLAGS) -DNO_FPINIT arithchk.c || \$(CC) -DNO_LONG_LONG \$(CFLAGS) -DNO_FPINIT arithchk.c\n";
    print WRITEMAKE "\t./a.out >arith.h\n";
    print WRITEMAKE "\trm -f a.out arithchk.o\n\n";

    print WRITEMAKE "check:\n\txsum Notice README abort_.c arithchk.c backspac.c c_abs.c c_cos.c \\\nc_div.c c_exp.c c_log.c c_sin.c c_sqrt.c cabs.c close.c comptry.bat \\\nd_abs.c d_acos.c d_asin.c d_atan.c d_atn2.c d_cnjg.c d_cos.c d_cosh.c \\\nd_dim.c d_exp.c d_imag.c d_int.c d_lg10.c d_log.c d_mod.c \\\nd_nint.c d_prod.c d_sign.c d_sin.c d_sinh.c d_sqrt.c d_tan.c \\\nd_tanh.c derf_.c derfc_.c dfe.c dolio.c dtime_.c due.c ef1asc_.c \\\n";
    print WRITEMAKE "ef1cmc_.c endfile.c erf_.c erfc_.c err.c etime_.c exit_.c f2c.h0 \\\nf2ch.add f77_aloc.c f77vers.c fio.h fmt.c fmt.h fmtlib.c \\\nfp.h ftell_.c \\\ngetarg_.c getenv_.c h_abs.c h_dim.c h_dnnt.c h_indx.c h_len.c \\\nh_mod.c h_nint.c h_sign.c hl_ge.c hl_gt.c hl_le.c hl_lt.c \\\ni77vers.c i_abs.c i_dim.c i_dnnt.c i_indx.c i_len.c i_mod.c \\\ni_nint.c i_sign.c iargc_.c iio.c ilnw.c inquire.c l_ge.c l_gt.c \\\nl_le.c l_lt.c lbitbits.c lbitshft.c libf2c.lbc libf2c.sy lio.h \\\nlread.c lwrite.c main.c makefile.sy makefile.u makefile.vc \\\nmakefile.wat math.hvc mkfile.plan9 open.c pow_ci.c pow_dd.c \\\npow_di.c pow_hh.c pow_ii.c pow_qq.c pow_ri.c pow_zi.c pow_zz.c \\\nqbitbits.c qbitshft.c r_abs.c r_acos.c r_asin.c r_atan.c r_atn2.c \\\nr_cnjg.c r_cos.c r_cosh.c r_dim.c r_exp.c r_imag.c r_int.c r_lg10.c \\\n";
    print WRITEMAKE "r_log.c r_mod.c r_nint.c r_sign.c r_sin.c r_sinh.c r_sqrt.c \\\nr_tan.c r_tanh.c rawio.h rdfmt.c rewind.c rsfe.c rsli.c rsne.c \\\ns_cat.c s_cmp.c s_copy.c s_paus.c s_rnge.c s_stop.c scomptry.bat \\\nsfe.c sig_die.c signal1.h0 signal_.c sue.c system_.c typesize.c \\\nuio.c uninit.c util.c wref.c wrtfmt.c wsfe.c wsle.c wsne.c xwsne.c \\\nz_abs.c z_cos.c z_div.c z_exp.c z_log.c z_sin.c z_sqrt.c >xsum1.out\n";
    print WRITEMAKE "\tcmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out\n";
    close(WRITEMAKE);


    open(WRITEMAKE, ">lapack/make.inc") or die;
    print WRITEMAKE "####################################################################\n#  LAPACK make include file.                                       #\n#  LAPACK, Version 3.0                                             #\n#  June 30, 1999                                                  #\n####################################################################\n#\nSHELL = /bin/sh\n";
    print WRITEMAKE "cc = $cc\n";
    print WRITEMAKE "CC = $CC\n";
    print WRITEMAKE "LINUX_BUILD = -Xlinker -defsym -Xlinker MAIN__=main\n";

    print WRITEMAKE "DRVOPTS  = \$(OPTS)\nNOOPT    = -u -f\nLOADER   = \$(CC)\n";
    print WRITEMAKE "ARCH     = ar\nARCHFLAGS= cr\nRANLIB   = ranlib\n";
    print WRITEMAKE "LAPACKLIB    = liblapack.a\n";
    close(WRITEMAKE);

    open(WRITEMAKE, ">lapack/Makefile") or die;
    print WRITEMAKE "LIN_SRC_DIR = ./\ninclude \$(LIN_SRC_DIR)make.inc\n\n\n";

    print WRITEMAKE "ALLAUX =  \$(LIN_SRC_DIR)xerbla.o \\\n\$(LIN_SRC_DIR)lsame.o \\\n\$(LIN_SRC_DIR)idamax.o \\\n\$(LIN_SRC_DIR)dswap.o \\\n\$(LIN_SRC_DIR)dcopy.o \\\n\$(LIN_SRC_DIR)dgemm.o \\\n\$(LIN_SRC_DIR)dtrmm.o \\\n\$(LIN_SRC_DIR)dgemv.o \\\n\$(LIN_SRC_DIR)dtrmv.o \\\n\$(LIN_SRC_DIR)dger.o \\\n\$(LIN_SRC_DIR)ddot.o \\\n\$(LIN_SRC_DIR)daxpy.o \\\n\$(LIN_SRC_DIR)ilaenv.o \\\n\$(LIN_SRC_DIR)ieeeck.o \\\n\$(LIN_SRC_DIR)dbdsqr.o \\\n\$(LIN_SRC_DIR)dorglq.o \\\n\$(LIN_SRC_DIR)dgelqf.o \\\n\$(LIN_SRC_DIR)dgebrd.o \\\n\$(LIN_SRC_DIR)dorgbr.o \\\n\$(LIN_SRC_DIR)dgeqrf.o \\\n\$(LIN_SRC_DIR)dormbr.o \\\n\$(LIN_SRC_DIR)dlasr.o \\\n\$(LIN_SRC_DIR)dlas2.o \\\n\$(LIN_SRC_DIR)dlasv2.o \\\n\$(LIN_SRC_DIR)dgelq2.o \\\n\$(LIN_SRC_DIR)dlabrd.o \\\n\$(LIN_SRC_DIR)dlasq1.o \\\n\$(LIN_SRC_DIR)dorgl2.o \\\n\$(LIN_SRC_DIR)dgebd2.o \\\n\$(LIN_SRC_DIR)dlabrd.o \\\n\$(LIN_SRC_DIR)dbdsdc.o \\\n\$(LIN_SRC_DIR)dgeqr2.o \\\n\$(LIN_SRC_DIR)dormqr.o \\\n\$(LIN_SRC_DIR)dormlq.o \\\n\$(LIN_SRC_DIR)dlasrt.o \\\n\$(LIN_SRC_DIR)dlasq2.o \\\n\$(LIN_SRC_DIR)dlasq3.o \\\n\$(LIN_SRC_DIR)dorml2.o \\\n\$(LIN_SRC_DIR)dorm2r.o \\\n\$(LIN_SRC_DIR)dlasd0.o \\\n\$(LIN_SRC_DIR)dlanst.o \\\n\$(LIN_SRC_DIR)dlasdq.o \\\n\$(LIN_SRC_DIR)dlasda.o \\\n\$(LIN_SRC_DIR)dlasq4.o \\\n\$(LIN_SRC_DIR)dlasq5.o \\\n\$(LIN_SRC_DIR)dlasq6.o \\\n\$(LIN_SRC_DIR)dlasdt.o \\\n\$(LIN_SRC_DIR)dlasd1.o \\\n\$(LIN_SRC_DIR)dlasd2.o \\\n\$(LIN_SRC_DIR)dlasd3.o \\\n\$(LIN_SRC_DIR)dlasd6.o \\\n\$(LIN_SRC_DIR)dlamrg.o \\\n\$(LIN_SRC_DIR)dlasd4.o \\\n\$(LIN_SRC_DIR)dlasd7.o \\\n\$(LIN_SRC_DIR)dlasd8.o \\\n\$(LIN_SRC_DIR)dlasd5.o \\\n\$(LIN_SRC_DIR)dlaed6.o \\\n\n";

    print WRITEMAKE "DZLAUX = \\\n\$(LIN_SRC_DIR)dlabad.o  \$(LIN_SRC_DIR)dlacpy.o \\\n\$(LIN_SRC_DIR)dlamch.o  \\\n\$(LIN_SRC_DIR)dscal.o  \\\n\$(LIN_SRC_DIR)dlascl.o  \\\n\$(LIN_SRC_DIR)drot.o  \\\n\$(LIN_SRC_DIR)dnrm2.o  \\\n\$(LIN_SRC_DIR)dhseqr.o \\\n\$(LIN_SRC_DIR)dlapy2.o \\\n\$(LIN_SRC_DIR)dlahrd.o \\\n\$(LIN_SRC_DIR)dlarfx.o \\\n\$(LIN_SRC_DIR)dlahqr.o \\\n\$(LIN_SRC_DIR)dlarf.o \\\n\$(LIN_SRC_DIR)dlarfg.o \\\n\$(LIN_SRC_DIR)dlanhs.o \\\n\$(LIN_SRC_DIR)dlanv2.o \\\n\$(LIN_SRC_DIR)dorghr.o \\\n\$(LIN_SRC_DIR)dtrevc.o \\\n\$(LIN_SRC_DIR)dlaln2.o \\\n\$(LIN_SRC_DIR)dladiv.o \\\n\$(LIN_SRC_DIR)dorgqr.o \\\n\$(LIN_SRC_DIR)dlassq.o \\\n\$(LIN_SRC_DIR)dorg2r.o \\\n\$(LIN_SRC_DIR)dlarft.o \\\n\$(LIN_SRC_DIR)dlaset.o \\\n\$(LIN_SRC_DIR)dgetri.o \\\n\$(LIN_SRC_DIR)dtrtri.o \\\n\$(LIN_SRC_DIR)dtrsm.o \\\n\$(LIN_SRC_DIR)dtrti2.o \\\n\$(LIN_SRC_DIR)dgetrf.o \\\n\$(LIN_SRC_DIR)dgetf2.o \\\n\$(LIN_SRC_DIR)dlaswp.o \\\n\n\n";


    print WRITEMAKE "DLASRC = \\\n\$(LIN_SRC_DIR)dgebal.o \\\n\$(LIN_SRC_DIR)dgeev.o  \\\n\$(LIN_SRC_DIR)dgehrd.o \\\n\$(LIN_SRC_DIR)dlange.o \\\n\$(LIN_SRC_DIR)dlartg.o \\\n\$(LIN_SRC_DIR)dlarfb.o \\\n\$(LIN_SRC_DIR)dgebak.o \\\n\$(LIN_SRC_DIR)dgehd2.o \\\n\$(LIN_SRC_DIR)dgesdd.o \\\n\$(LIN_SRC_DIR)dgesvd.o\n";





    print WRITEMAKE "all: double\n\n";


    print WRITEMAKE "double:  \$(DLASRC) \$(ALLAUX) \$(DZLAUX)\n";
    print WRITEMAKE "\t\$(ARCH) \$(ARCHFLAGS) \$(LAPACKLIB)  \$(DLASRC) \$(ALLAUX) \\\n\$(DZLAUX)\n";
    print WRITEMAKE "\t\$(RANLIB) \$(LAPACKLIB)\n";


    print WRITEMAKE "%.o: %.c\n";
    print WRITEMAKE "\t\$(cc) -O1 -c \$(LIN_SRC_DIR)\$<\n";
    close(WRITEMAKE);
}


open(WRITEMAKE, ">ranlib/Makefile") or die;
print WRITEMAKE "# Unix makefile for random number generator: see README.\n";
print WRITEMAKE ".SUFFIXES: .c .o\n";
print WRITEMAKE "cc = $CC\n";

print WRITEMAKE "SHELL = /bin/sh\n";
print WRITEMAKE "CFLAGS = -O\n";
print WRITEMAKE "O = o\n";
print WRITEMAKE "ARCH     = ar\nARCHFLAGS= cr\n";
print WRITEMAKE "RANDLIB    = libranlib.a\n";
print WRITEMAKE "all: \$(RANDLIB)\n\n";


print WRITEMAKE "\$(RANDLIB):  linpack.\$(O) com.\$(O) ranlib.\$(O)\n";
print WRITEMAKE "\t\$(ARCH) \$(ARCHFLAGS) \$(RANDLIB) linpack.\$(O) com.\$(O) ranlib.\$(O)\n";
print WRITEMAKE "\tranlib \$(RANDLIB)\n";


print WRITEMAKE "%.o: %.c\n";
print WRITEMAKE "\t\$(cc) -O1 -c \$<\n";
close(WRITEMAKE);

open(WRITEMAKE, ">src/Makefile") or die;

print WRITEMAKE "#Makefile for POInT core routines\n";
print WRITEMAKE "#G. Conant, 4/26/19\n\n";
print WRITEMAKE "cc = $cc\n";
print WRITEMAKE "CC = $CC\n";

print WRITEMAKE "O = o\nSRC_DIR = ./\n";
if($lapack_installed == 0) {
    print WRITEMAKE "LIBRARY_DIR = -L../lapack \\\n-L../libf2c \\\n -L../ranlib\\\n";
}
else {
    print WRITEMAKE "LIBRARY_DIR = -L../ranlib\\\n";
}
if ($plot_path eq "") { print WRITEMAKE "-L.\n";}
else {print WRITEMAKE "-L.\\\n-L$plot_path\n";}

if($lapack_installed == 0) {
    if ($parallel == 0) {
        if ($plot_path eq "") {
            print WRITEMAKE "CFLAGS = -g -D___LOCAL_BLAS___  -I./\n";
        }
        else {
            print WRITEMAKE "CFLAGS = -g -D_DO_PLOT_ -D___LOCAL_BLAS___  -I./\n";
        }
    }
    else {
        if ($plot_path eq "") {
            print WRITEMAKE "CFLAGS = -DGCC_COMPILE -D_OPEN_MP_VERSION_ -fopenmp -D___LOCAL_BLAS___\n"
        }
        else {
            print WRITEMAKE "CFLAGS =-D_DO_PLOT_ -DGCC_COMPILE -D_OPEN_MP_VERSION_ -fopenmp -D___LOCAL_BLAS___ -I.\n"
        }
    }
}
else {
    if ($parallel == 0) {
        if ($plot_path eq "") {
            print WRITEMAKE "CFLAGS = -g   -I./\n";
        }
        else {
            print WRITEMAKE "CFLAGS = -g -D_DO_PLOT_   -I./\n";
        }
    }
    else {
        if ($plot_path eq "") {
            print WRITEMAKE "CFLAGS = -DGCC_COMPILE -D_OPEN_MP_VERSION_ -fopenmp \n"
        }
        else {
            print WRITEMAKE "CFLAGS =-D_DO_PLOT_ -DGCC_COMPILE -D_OPEN_MP_VERSION_ -fopenmp  -I.\n"
        }
    }
}
print WRITEMAKE "OPTIM_SPEED = -O3\nOPTIM_SIZE = -O1\nMATH_LIB = -lm\n";
if ($lapack_version eq "") {
    print WRITEMAKE "LAPACK_LIB = -llapack\nBLAS_LIB = -lblas\n";
}
else {
    print WRITEMAKE "LAPACK_LIB = -l:", $lapack_version, "\nBLAS_LIB = -l:", $blas_version, "\n";
}

print WRITEMAKE "F2C_LIB = -lf2c\nDNAFUNCS_LIB= -ldnafuncs\nLIKELIHOOD_LIB = -llikelihood\nRANDLIB = -lranlib\n";
if ($plot_path ne "") {
    if ($plot_version eq "") {
        print WRITEMAKE "PLOT_LIB=-lplot\n";
    }
    else {
        print WRITEMAKE "PLOT_LIB=-l:", $plot_version, "\n";
    }
}

print WRITEMAKE "OPTIONS = \$(CFLAGS) \$(INCLUDE)\n";


if ($make_daemon ==0) {
    print WRITEMAKE "all:  libdnafuncs.a liblikelihood.a POInT POInT_genome_scaffold POInT_ances_order POInT_simulate\n";
    print WRITEMAKE "install:\n\tln -s ../POInT /usr/local/bin\n";
}
else {
    print WRITEMAKE "all:  libdnafuncs.a liblikelihood.a POInT POInT_genome_scaffold POInT_ances_order POInT_simulate POInT_daemon POInT_browser POInT_download\n";
    print WRITEMAKE "install:\n\tln -s ../POInT /usr/local/bin\n";
}
print WRITEMAKE "DNALIBRARY = libdnafuncs.a\nDNALIBRARY_OBJS =        \\\n\tread_seq.\$(O)    \\\n\tgen_dna_funcs.\$(O) \\\n\tscore_matrix.\$(O) \\\n";
print WRITEMAKE "\texchange.\$(O)   \\\n\tgen_code.\$(O)   \\\n\twrite_seq.\$(O)\n";

print WRITEMAKE "\n\n\nlibdnafuncs : \$(DNALIBRARY)\n";

print WRITEMAKE "\$(DNALIBRARY): \$(DNALIBRARY_OBJS)\n";
print WRITEMAKE "\tar  cr \$(DNALIBRARY) \$(DNALIBRARY_OBJS)\n";
print WRITEMAKE "\tranlib \$(DNALIBRARY)\n";

if($lapack_installed == 0) {
    print WRITEMAKE "\n\n\nLIN_SRC_DIR= ../lapack/\n";

    print WRITEMAKE "ALLAUX =  \$(LIN_SRC_DIR)xerbla.o \\\n\$(LIN_SRC_DIR)lsame.o \\\n\$(LIN_SRC_DIR)idamax.o \\\n\$(LIN_SRC_DIR)dswap.o \\\n\$(LIN_SRC_DIR)dcopy.o \\\n\$(LIN_SRC_DIR)dgemm.o \\\n\$(LIN_SRC_DIR)dtrmm.o \\\n\$(LIN_SRC_DIR)dgemv.o \\\n\$(LIN_SRC_DIR)dtrmv.o \\\n\$(LIN_SRC_DIR)dger.o \\\n\$(LIN_SRC_DIR)ddot.o \\\n\$(LIN_SRC_DIR)daxpy.o \\\n\$(LIN_SRC_DIR)ilaenv.o \\\n\$(LIN_SRC_DIR)ieeeck.o \\\n\$(LIN_SRC_DIR)dbdsqr.o \\\n\$(LIN_SRC_DIR)dorglq.o \\\n\$(LIN_SRC_DIR)dgelqf.o \\\n\$(LIN_SRC_DIR)dgebrd.o \\\n\$(LIN_SRC_DIR)dorgbr.o \\\n\$(LIN_SRC_DIR)dgeqrf.o \\\n\$(LIN_SRC_DIR)dormbr.o \\\n\$(LIN_SRC_DIR)dlasr.o \\\n\$(LIN_SRC_DIR)dlas2.o \\\n\$(LIN_SRC_DIR)dlasv2.o \\\n\$(LIN_SRC_DIR)dgelq2.o \\\n\$(LIN_SRC_DIR)dlabrd.o \\\n\$(LIN_SRC_DIR)dlasq1.o \\\n\$(LIN_SRC_DIR)dorgl2.o \\\n\$(LIN_SRC_DIR)dgebd2.o \\\n\$(LIN_SRC_DIR)dlabrd.o \\\n\$(LIN_SRC_DIR)dbdsdc.o \\\n\$(LIN_SRC_DIR)dgeqr2.o \\\n\$(LIN_SRC_DIR)dormqr.o \\\n\$(LIN_SRC_DIR)dormlq.o \\\n\$(LIN_SRC_DIR)dlasrt.o \\\n\$(LIN_SRC_DIR)dlasq2.o \\\n\$(LIN_SRC_DIR)dlasq3.o \\\n\$(LIN_SRC_DIR)dorml2.o \\\n\$(LIN_SRC_DIR)dorm2r.o \\\n\$(LIN_SRC_DIR)dlasd0.o \\\n\$(LIN_SRC_DIR)dlanst.o \\\n\$(LIN_SRC_DIR)dlasdq.o \\\n\$(LIN_SRC_DIR)dlasda.o \\\n\$(LIN_SRC_DIR)dlasq4.o \\\n\$(LIN_SRC_DIR)dlasq5.o \\\n\$(LIN_SRC_DIR)dlasq6.o \\\n\$(LIN_SRC_DIR)dlasdt.o \\\n\$(LIN_SRC_DIR)dlasd1.o \\\n\$(LIN_SRC_DIR)dlasd2.o \\\n\$(LIN_SRC_DIR)dlasd3.o \\\n\$(LIN_SRC_DIR)dlasd6.o \\\n\$(LIN_SRC_DIR)dlamrg.o \\\n\$(LIN_SRC_DIR)dlasd4.o \\\n\$(LIN_SRC_DIR)dlasd7.o \\\n\$(LIN_SRC_DIR)dlasd8.o \\\n\$(LIN_SRC_DIR)dlasd5.o \\\n\$(LIN_SRC_DIR)dlaed6.o \\\n\n";

    print WRITEMAKE "DZLAUX = \\\n\$(LIN_SRC_DIR)dlabad.o  \$(LIN_SRC_DIR)dlacpy.o \\\n\$(LIN_SRC_DIR)dlamch.o  \\\n\$(LIN_SRC_DIR)dscal.o  \\\n\$(LIN_SRC_DIR)dlascl.o  \\\n\$(LIN_SRC_DIR)drot.o  \\\n\$(LIN_SRC_DIR)dnrm2.o  \\\n\$(LIN_SRC_DIR)dhseqr.o \\\n\$(LIN_SRC_DIR)dlapy2.o \\\n\$(LIN_SRC_DIR)dlahrd.o \\\n\$(LIN_SRC_DIR)dlarfx.o \\\n\$(LIN_SRC_DIR)dlahqr.o \\\n\$(LIN_SRC_DIR)dlarf.o \\\n\$(LIN_SRC_DIR)dlarfg.o \\\n\$(LIN_SRC_DIR)dlanhs.o \\\n\$(LIN_SRC_DIR)dlanv2.o \\\n\$(LIN_SRC_DIR)dorghr.o \\\n\$(LIN_SRC_DIR)dtrevc.o \\\n\$(LIN_SRC_DIR)dlaln2.o \\\n\$(LIN_SRC_DIR)dladiv.o \\\n\$(LIN_SRC_DIR)dorgqr.o \\\n\$(LIN_SRC_DIR)dlassq.o \\\n\$(LIN_SRC_DIR)dorg2r.o \\\n\$(LIN_SRC_DIR)dlarft.o \\\n\$(LIN_SRC_DIR)dlaset.o \\\n\$(LIN_SRC_DIR)dgetri.o \\\n\$(LIN_SRC_DIR)dtrtri.o \\\n\$(LIN_SRC_DIR)dtrsm.o \\\n\$(LIN_SRC_DIR)dtrti2.o \\\n\$(LIN_SRC_DIR)dgetrf.o \\\n\$(LIN_SRC_DIR)dgetf2.o \\\n\$(LIN_SRC_DIR)dlaswp.o \\\n\n\n";


    print WRITEMAKE "DLASRC = \\\n\$(LIN_SRC_DIR)dgebal.o \\\n\$(LIN_SRC_DIR)dgeev.o  \\\n\$(LIN_SRC_DIR)dgehrd.o \\\n\$(LIN_SRC_DIR)dlange.o \\\n\$(LIN_SRC_DIR)dlartg.o \\\n\$(LIN_SRC_DIR)dlarfb.o \\\n\$(LIN_SRC_DIR)dgebak.o \\\n\$(LIN_SRC_DIR)dgehd2.o \\\n\$(LIN_SRC_DIR)dgesdd.o \\\n\$(LIN_SRC_DIR)dgesvd.o\n";

    print WRITEMAKE "\n\n\nF2C_DIR = ../libf2c/\n";
    print WRITEMAKE "MISC =  \$(F2C_DIR)f77vers.o \$(F2C_DIR)i77vers.o \$(F2C_DIR)main.o \$(F2C_DIR)s_rnge.o \\\n\$(F2C_DIR)abort_.o \$(F2C_DIR)exit_.o \$(F2C_DIR)getarg_.o \$(F2C_DIR)iargc_.o\\\n\$(F2C_DIR)getenv_.o \$(F2C_DIR)signal_.o \$(F2C_DIR)s_stop.o \$(F2C_DIR)s_paus.o \$(F2C_DIR)system_.o \$(F2C_DIR)cabs.o\\\n\$(F2C_DIR)derf_.o \$(F2C_DIR)derfc_.o \$(F2C_DIR)erf_.o \$(F2C_DIR)erfc_.o \$(F2C_DIR)sig_die.o \$(F2C_DIR)uninit.o\n";

    print WRITEMAKE "POW =   \$(F2C_DIR)pow_ci.o \$(F2C_DIR)pow_dd.o \$(F2C_DIR)pow_di.o \$(F2C_DIR)pow_hh.o\\\n\$(F2C_DIR)pow_ii.o \$(F2C_DIR)pow_ri.o \$(F2C_DIR)pow_zi.o \$(F2C_DIR)pow_zz.o\n";

    print WRITEMAKE "CX =    \$(F2C_DIR)c_abs.o \$(F2C_DIR)c_cos.o \$(F2C_DIR)c_div.o \$(F2C_DIR)c_exp.o \$(F2C_DIR)c_log.o \$(F2C_DIR)c_sin.o \$(F2C_DIR)c_sqrt.o\n";
    print WRITEMAKE "DCX =   \$(F2C_DIR)z_abs.o \$(F2C_DIR)z_cos.o \$(F2C_DIR)z_div.o \$(F2C_DIR)z_exp.o \$(F2C_DIR)z_log.o \$(F2C_DIR)z_sin.o \$(F2C_DIR)z_sqrt.o\n";
    print WRITEMAKE "REAL =  \$(F2C_DIR)r_abs.o \$(F2C_DIR)r_acos.o \$(F2C_DIR)r_asin.o \$(F2C_DIR)r_atan.o \$(F2C_DIR)r_atn2.o \$(F2C_DIR)r_cnjg.o \$(F2C_DIR)r_cos.o\\\n\$(F2C_DIR)r_cosh.o \$(F2C_DIR)r_dim.o \$(F2C_DIR)r_exp.o \$(F2C_DIR)r_imag.o \$(F2C_DIR)r_int.o\\\n\$(F2C_DIR)r_lg10.o \$(F2C_DIR)r_log.o \$(F2C_DIR)r_mod.o \$(F2C_DIR)r_nint.o \$(F2C_DIR)r_sign.o\\\n\$(F2C_DIR)r_sin.o \$(F2C_DIR)r_sinh.o \$(F2C_DIR)r_sqrt.o \$(F2C_DIR)r_tan.o \$(F2C_DIR)r_tanh.o\n";
    print WRITEMAKE "DBL =   \$(F2C_DIR)d_abs.o \$(F2C_DIR)d_acos.o \$(F2C_DIR)d_asin.o \$(F2C_DIR)d_atan.o \$(F2C_DIR)d_atn2.o\\\n\$(F2C_DIR)d_cnjg.o \$(F2C_DIR)d_cos.o \$(F2C_DIR)d_cosh.o \$(F2C_DIR)d_dim.o \$(F2C_DIR)d_exp.o\\\n\$(F2C_DIR)d_imag.o \$(F2C_DIR)d_int.o \$(F2C_DIR)d_lg10.o \$(F2C_DIR)d_log.o \$(F2C_DIR)d_mod.o\\\n\$(F2C_DIR)d_nint.o \$(F2C_DIR)d_prod.o \$(F2C_DIR)d_sign.o \$(F2C_DIR)d_sin.o \$(F2C_DIR)d_sinh.o\\\n\$(F2C_DIR)d_sqrt.o \$(F2C_DIR)d_tan.o \$(F2C_DIR)d_tanh.o\n";
    print WRITEMAKE "INT =   \$(F2C_DIR)i_abs.o \$(F2C_DIR)i_dim.o \$(F2C_DIR)i_dnnt.o \$(F2C_DIR)i_indx.o \\\n\$(F2C_DIR)i_len.o \$(F2C_DIR)i_mod.o \$(F2C_DIR)i_nint.o \$(F2C_DIR)i_sign.o\\\n\$(F2C_DIR)lbitbits.o \$(F2C_DIR)lbitshft.o\n";

    print WRITEMAKE "HALF =  \$(F2C_DIR)h_abs.o \$(F2C_DIR)h_dim.o \$(F2C_DIR)h_dnnt.o \$(F2C_DIR)h_indx.o \$(F2C_DIR)h_len.o \$(F2C_DIR)h_mod.o \$(F2C_DIR)h_nint.o \$(F2C_DIR)h_sign.o\n";
    print WRITEMAKE "CMP =   \$(F2C_DIR)l_ge.o \$(F2C_DIR)l_gt.o \$(F2C_DIR)l_le.o \$(F2C_DIR)l_lt.o \$(F2C_DIR)hl_ge.o \$(F2C_DIR)hl_gt.o \$(F2C_DIR)hl_le.o \$(F2C_DIR)hl_lt.o\n";
    print WRITEMAKE "EFL =   \$(F2C_DIR)ef1asc_.o \$(F2C_DIR)ef1cmc_.o\n";
    print WRITEMAKE "CHAR =  \$(F2C_DIR)f77_aloc.o \$(F2C_DIR)s_cat.o \$(F2C_DIR)s_cmp.o \$(F2C_DIR)s_copy.o\n";
    print WRITEMAKE "I77 =   \$(F2C_DIR)backspac.o \$(F2C_DIR)close.o \$(F2C_DIR)dfe.o \$(F2C_DIR)dolio.o \$(F2C_DIR)due.o \$(F2C_DIR)endfile.o \$(F2C_DIR)err.o\\\n\$(F2C_DIR)fmt.o \$(F2C_DIR)fmtlib.o \$(F2C_DIR)ftell_.o \$(F2C_DIR)iio.o \$(F2C_DIR)ilnw.o \$(F2C_DIR)inquire.o \$(F2C_DIR)lread.o \$(F2C_DIR)lwrite.o\\\n\$(F2C_DIR)open.o \$(F2C_DIR)rdfmt.o \$(F2C_DIR)rewind.o \$(F2C_DIR)rsfe.o \$(F2C_DIR)rsli.o \$(F2C_DIR)rsne.o \$(F2C_DIR)sfe.o \$(F2C_DIR)sue.o\\\n\$(F2C_DIR)typesize.o \$(F2C_DIR)uio.o \$(F2C_DIR)util.o \$(F2C_DIR)wref.o \$(F2C_DIR)wrtfmt.o \$(F2C_DIR)wsfe.o \$(F2C_DIR)wsle.o \$(F2C_DIR)wsne.o \$(F2C_DIR)xwsne.o\n";
    print WRITEMAKE "QINT =  \$(F2C_DIR)pow_qq.o \$(F2C_DIR)qbitbits.o \$(F2C_DIR)qbitshft.o\n";
    print WRITEMAKE "TIME =  \$(F2C_DIR)dtime_.o \$(F2C_DIR)etime_.o\n";
}

print WRITEMAKE "\n\n\nLIKELIHOODLIBRARY = liblikelihood.a\n";


print WRITEMAKE "LIKELIHOODLIBRARY_OBJS =        \\\n\tmaxlike.\$(O)    \\\n\tcodon_like.\$(O) \\\n\tnucleotide_like.\$(O)   \\\n\ttree.\$(O)          \\\n\tpowell.\$(O)       \\\n\tlin_alg.\$(O)    \\\n\tread_tree.\$(O)  \\\n\twrite_tree.\$(O)\\\n\tsearch_trees.\$(O)\\\n\tother_like.\$(O)\n";


print WRITEMAKE "\n\n\nliblikelihood : \$(LIKELIHOODLIBRARY)\n";

print WRITEMAKE "\$(LIKELIHOODLIBRARY): \$(LIKELIHOODLIBRARY_OBJS)\n";
if($lapack_installed == 0) {
    print WRITEMAKE "\tar cr \$(LIKELIHOODLIBRARY) \$(LIKELIHOODLIBRARY_OBJS) \$(DLASRC) \$(ALLAUX) \$(DZLAUX) \$(MISC) \$(POW) \$(CX) \$(DCX) \$(REAL) \$(DBL) \$(INT)  \$(HALF) \$(CMP) \$(EFL) \$(CHAR) \$(I77) \$(TIME)\n";
}
else {
    print WRITEMAKE "\tar cr \$(LIKELIHOODLIBRARY) \$(LIKELIHOODLIBRARY_OBJS)\n";
}
print WRITEMAKE "\tranlib  \$(LIKELIHOODLIBRARY)\n";

print WRITEMAKE "SCAFFOLD_WGX_OBJS = scaffold_WGX.\$(O) random.\$(O)\n";
print WRITEMAKE "\n\nPOInT_genome_scaffold: \$(SCAFFOLD_WGX_OBJS) \$(LIKELIHOOD_LIB) \$(DNAFUNCS_LIB)\n";
print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_genome_scaffold \$(OPTIONS) \$(SCAFFOLD_WGX_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(RANDLIB)\n";

print WRITEMAKE "minimize_track_breaks.\$(O): minimize_track_breaks.cpp track_anneal_recombo.cpp\n";

print WRITEMAKE "MINIMIZE_TRACK_BREAKS_OBJS =    minimize_track_breaks.\$(O) phylo_model_matrix.\$(O) genome_list.\$(O) genome_tripl_list.\$(O) random.\$(O) nrutil.\$(O)\\\n";
print WRITEMAKE "\tgenome_ploidy_like.\$(O)\n";
#if ($plot_path eq "") {print WRITEMAKE "\tgenome_ploidy_like.\$(O)\n";}
#else {print WRITEMAKE "\tdraw_tracking_WGX.\$(O)\\\n\tgenome_ploidy_like.\$(O)\n";}



print WRITEMAKE "\n\nPOInT_ances_order: \$(MINIMIZE_TRACK_BREAKS_OBJS) \$(LIKELIHOOD_LIB) \$(DNAFUNCS_LIB)\n";
if($plot_path eq "") {
    if($lapack_installed == 0) {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_ances_order \$(OPTIONS) \$(MINIMIZE_TRACK_BREAKS_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(RANDLIB)\n";
    }
    else {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_ances_order \$(OPTIONS) \$(MINIMIZE_TRACK_BREAKS_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(LAPACK_LIB) \$(BLAS_LIB) \$(RANDLIB)\n";
    }
}
else {
    if($lapack_installed == 0) {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_ances_order \$(OPTIONS) \$(MINIMIZE_TRACK_BREAKS_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(RANDLIB) \$(PLOT_LIB)\n";
    }
    else {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_ances_order \$(OPTIONS) \$(MINIMIZE_TRACK_BREAKS_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(LAPACK_LIB) \$(BLAS_LIB) \$(RANDLIB) \$(PLOT_LIB)\n";
    }
}


print WRITEMAKE "SIM_WGX_OBJS =    sim_WGX.\$(O) sim_data.\$(O) phylo_model_matrix.\$(O) genome_list.\$(O) genome_tripl_list.\$(O) random.\$(O) nrutil.\$(O)\\\n";
print WRITEMAKE "\tgenome_ploidy_like.\$(O)\n";
#if ($plot_path eq "") {print WRITEMAKE "\tgenome_ploidy_like.\$(O)\n";}
#else {print WRITEMAKE "\tdraw_tracking_WGX.\$(O)\\\n\tgenome_ploidy_like.\$(O)\n";}



print WRITEMAKE "\n\nPOInT_simulate: \$(SIM_WGX_OBJS) \$(LIKELIHOOD_LIB) \$(DNAFUNCS_LIB)\n";
if($plot_path eq "") {
    if($lapack_installed == 0) {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_simulate \$(OPTIONS) \$(SIM_WGX_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(RANDLIB)\n";
    }
    else {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_simulate \$(OPTIONS) \$(SIM_WGX_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(LAPACK_LIB) \$(BLAS_LIB) \$(RANDLIB)\n";
    }
}
else {
    if($lapack_installed == 0) {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_simulate \$(OPTIONS) \$(SIM_WGX_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(RANDLIB) \$(PLOT_LIB)\n";
    }
    else {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT_simulate \$(OPTIONS) \$(SIM_WGX_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(LAPACK_LIB) \$(BLAS_LIB) \$(RANDLIB) \$(PLOT_LIB)\n";
    }
}


print WRITEMAKE "SEARCH_WGX_ASSIGNS_OBJS =       nrutil.\$(O)\\\n\tsearch_WGX_assigns.\$(O)\\\n\tgenome_list.\$(O) \\\n\tgenome_tripl_list.\$(O)\\\n\tphylo_model_matrix.\$(O) \\\n";
if ($plot_path eq "") {print WRITEMAKE "\tgenome_ploidy_like.\$(O)\n";}
else {print WRITEMAKE "\tdraw_tracking_WGX.\$(O)\\\n\tgenome_ploidy_like.\$(O)\n";}

print WRITEMAKE "\n\n\nPOInT: \$(SEARCH_WGX_ASSIGNS_OBJS) \$(LIKELIHOOD_LIB) \$(DNAFUNCS_LIB)\n";
if($plot_path eq "") {
        if($lapack_installed == 0) {
            print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT \$(OPTIONS) \$(SEARCH_WGX_ASSIGNS_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB)\n";
        }
        else {
            print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT \$(OPTIONS) \$(SEARCH_WGX_ASSIGNS_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(LAPACK_LIB) \$(BLAS_LIB)\n";
        }
}
else {
    if($lapack_installed == 0) {
        print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT \$(OPTIONS) \$(SEARCH_WGX_ASSIGNS_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(PLOT_LIB)\n";
    }
    else {
       print WRITEMAKE "\t\$(CC) \$(LIBRARY_DIR) -o ../POInT \$(OPTIONS) \$(SEARCH_WGX_ASSIGNS_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB)  \$(DNAFUNCS_LIB) \$(LAPACK_LIB) \$(BLAS_LIB) \$(PLOT_LIB)\n";
    }
}

if (($plot_path ne "") && ($make_daemon == 1)) {
    print WRITEMAKE "search_WGX_assigns_daemon.\$(O): search_WGX_assigns.cpp\n\t\$(CC)  \$(OPTIONS) -DPOInT_daemon \$(OPTIM_SPEED) -c \$< -o search_WGX_assigns_daemon.\$(O)\n\n";
    print WRITEMAKE "POINT_DAEMON_OBJS =      nrutil.\$(O)  genome_tripl_list.\$(O)  phylo_model_matrix.\$(O) tree_plot.\$(O) ";
    print WRITEMAKE " genome_ploidy_like.\$(O)     genome_list.\$(O)  POInT_IPC.\$(O)  draw_tracking_WGX.\$(O) search_WGX_assigns_daemon.\$(O)\n";
    print WRITEMAKE "POInT_daemon: \$(POINT_DAEMON_OBJS) \$(LIKELIHOOD_LIB)\n";
    if($lapack_installed == 0) {
        print WRITEMAKE "\t\$(CC) -std=c++11 -Xlinker -defsym -Xlinker MAIN__=main \$(LIBRARY_DIR) -DPOInT_daemon   -o ../POInT_daemon \$(OPTIONS) \$(POINT_DAEMON_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB) \$(DNAFUNCS_LIB) \$(PLOT_LIB)\n";
    }
    else {
        print WRITEMAKE "\t\$(CC) -std=c++11 -Xlinker -defsym -Xlinker MAIN__=main \$(LIBRARY_DIR) -DPOInT_daemon   -o ../POInT_daemon \$(OPTIONS) \$(POINT_DAEMON_OBJS) \$(DNAFUNCS_LIB) \$(LIKELIHOOD_LIB) \$(DNAFUNCS_LIB) \$(LAPACK_LIB) \$(BLAS_LIB) \$(PLOT_LIB)\n";
    }
    
    print WRITEMAKE "POInT_browser:\n\t\$(CC) -I /usr/include/cgicc/   POInT_browse.cpp -lcgicc -o ../POInT_browse\n";
    print WRITEMAKE "POInT_download:\n\t\$(CC) -I /usr/include/cgicc/  POInT_download.cpp -lcgicc -o ../POInT_download -lboost_system -DBOOST_NO_CXX11_SCOPED_ENUMS -lboost_filesystem\n";
}

if ($plot_path ne "") {
	print WRITEMAKE "\n\ndraw_tracking_WGX.\$(O): draw_tracking_WGX.cpp\n\t\$(CC) -std=c++11  \$(OPTIONS) -c draw_tracking_WGX.cpp\n\n";
	if ($make_daemon == 1) {
        print WRITEMAKE "\n\nPOInT_IPC.\$(O): POInT_IPC.cpp\n\t\$(CC) -std=c++11  \$(OPTIONS) -c POInT_IPC.cpp\n";
	}
}

print WRITEMAKE "\n\n\n%.o: %.cpp\n";
print WRITEMAKE "\t\$(CC)  \$(OPTIONS) \$(OPTIM_SPEED) -c \$<\n";


print WRITEMAKE "%.o: %.c\n";
print WRITEMAKE "\t\$(CC)  \$(OPTIONS) \$(OPTIM_SPEED) -c \$<\n";
close(WRITEMAKE);

open(WRITEMAKE, ">Makefile") or die;
print WRITEMAKE "#Master Makefile for POInT--uses MAkefiles generated by configure.pl to compile required libraries and source files\n#G. Conant 04/25/19\n\n";


print WRITEMAKE "default: all\n";

if($lapack_installed == 0) {
    print WRITEMAKE "all: ./lapack/liblapack.a ./libf2c/libf2c.a ./ranlib/libranlib.a progs\n";
    print WRITEMAKE "./lapack/liblapack.a:\n";
    print WRITEMAKE "\tcd lapack; make\n\n";
    
    print WRITEMAKE "./libf2c/libf2c.a:\n";
    print WRITEMAKE "\tcd libf2c; make\n\n";
    
    print WRITEMAKE "./ranlib/libranlib.a:\n";
    print WRITEMAKE "\tcd ranlib; make\n\n";
    
    print WRITEMAKE "progs: ./lapack/liblapack.a ./libf2c/libf2c.a ./ranlib/libranlib.a\n";
    print WRITEMAKE "\tcd src;  make\n";
    
    print WRITEMAKE "clean:\n";
    print WRITEMAKE "\trm -f POInT  src/*.o src/*.a ranlib/*.o ranlib/*.a lapack/*.o lapack/*.a lapack/ libf2c/*.o libf2c/*.a \n";
    
}
else {
    print WRITEMAKE "all:  progs ./ranlib/libranlib.a \n";
    print WRITEMAKE "./ranlib/libranlib.a:\n";
    print WRITEMAKE "\tcd ranlib; make\n\n";
    
    print WRITEMAKE "progs: ./ranlib/libranlib.a\n";
    print WRITEMAKE "\tcd src;  make\n";
    print WRITEMAKE "clean:\n";
    print WRITEMAKE "\trm -f POInT s ranlib/*.o ranlib/*.a src/*.o src/*.a\n";
}









close(WRITEMAKE);
