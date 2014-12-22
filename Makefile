# This Makefile is for the FEELnc extension to perl.
#
# It was generated automatically by MakeMaker version
# 6.66 (Revision: 66600) from the contents of
# Makefile.PL. Don't edit this file, edit Makefile.PL instead.
#
#       ANY CHANGES MADE HERE WILL BE LOST!
#
#   MakeMaker ARGV: (q[PREFIX=/home/genouest/umr6061/recomgen/tderrien/bin/perl/FEELnc_v0.1])
#

#   MakeMaker Parameters:

#     AUTHOR => [q[ Fabrice Legeai ; Audrey David ; Thomas Derrien]]
#     BUILD_REQUIRES => {  }
#     CONFIGURE_REQUIRES => {  }
#     NAME => q[FEELnc]
#     PREREQ_PM => { Bio::SeqFeature::Generic=>q[0], Bio::Tools::GFF=>q[0], Parallel::ForkManager=>q[1.07] }
#     TEST_REQUIRES => {  }
#     VERSION_FROM => q[lib/FEELnc.pm]

# --- MakeMaker post_initialize section:


# --- MakeMaker const_config section:

# These definitions are from config.sh (via /local/perl/5.18.2/lib/5.18.2/x86_64-linux/Config.pm).
# They may have been overridden via Makefile.PL or on the command line.
AR = ar
CC = cc
CCCDLFLAGS = -fPIC
CCDLFLAGS = -Wl,-E
DLEXT = so
DLSRC = dl_dlopen.xs
EXE_EXT = 
FULL_AR = /usr/bin/ar
LD = cc
LDDLFLAGS = -shared -O2 -L/usr/local/lib -fstack-protector
LDFLAGS =  -fstack-protector -L/usr/local/lib
LIBC = /lib/libc-2.5.so
LIB_EXT = .a
OBJ_EXT = .o
OSNAME = linux
OSVERS = 2.6.18-194.11.4.el5
RANLIB = :
SITELIBEXP = /local/perl/5.18.2/lib/site_perl/5.18.2
SITEARCHEXP = /local/perl/5.18.2/lib/site_perl/5.18.2/x86_64-linux
SO = so
VENDORARCHEXP = $(VENDORPREFIX)/lib/5.18.2/x86_64-linux
VENDORLIBEXP = $(VENDORPREFIX)/lib


# --- MakeMaker constants section:
AR_STATIC_ARGS = cr
DIRFILESEP = /
DFSEP = $(DIRFILESEP)
NAME = FEELnc
NAME_SYM = FEELnc
VERSION = 0.1
VERSION_MACRO = VERSION
VERSION_SYM = 0_1
DEFINE_VERSION = -D$(VERSION_MACRO)=\"$(VERSION)\"
XS_VERSION = 0.1
XS_VERSION_MACRO = XS_VERSION
XS_DEFINE_VERSION = -D$(XS_VERSION_MACRO)=\"$(XS_VERSION)\"
INST_ARCHLIB = blib/arch
INST_SCRIPT = blib/script
INST_BIN = blib/bin
INST_LIB = blib/lib
INST_MAN1DIR = blib/man1
INST_MAN3DIR = blib/man3
MAN1EXT = 1
MAN3EXT = 3
INSTALLDIRS = site
DESTDIR = 
PREFIX = /home/genouest/umr6061/recomgen/tderrien/bin/perl/FEELnc_v0.1
PERLPREFIX = $(PREFIX)
SITEPREFIX = $(PREFIX)
VENDORPREFIX = $(PREFIX)
INSTALLPRIVLIB = $(PERLPREFIX)/lib/5.18.2
DESTINSTALLPRIVLIB = $(DESTDIR)$(INSTALLPRIVLIB)
INSTALLSITELIB = $(SITEPREFIX)/lib/site_perl/5.18.2
DESTINSTALLSITELIB = $(DESTDIR)$(INSTALLSITELIB)
INSTALLVENDORLIB = $(VENDORPREFIX)/lib
DESTINSTALLVENDORLIB = $(DESTDIR)$(INSTALLVENDORLIB)
INSTALLARCHLIB = $(PERLPREFIX)/lib/5.18.2/x86_64-linux
DESTINSTALLARCHLIB = $(DESTDIR)$(INSTALLARCHLIB)
INSTALLSITEARCH = $(SITEPREFIX)/lib/site_perl/5.18.2/x86_64-linux
DESTINSTALLSITEARCH = $(DESTDIR)$(INSTALLSITEARCH)
INSTALLVENDORARCH = $(VENDORPREFIX)/lib/5.18.2/x86_64-linux
DESTINSTALLVENDORARCH = $(DESTDIR)$(INSTALLVENDORARCH)
INSTALLBIN = $(PERLPREFIX)/bin
DESTINSTALLBIN = $(DESTDIR)$(INSTALLBIN)
INSTALLSITEBIN = $(SITEPREFIX)/bin
DESTINSTALLSITEBIN = $(DESTDIR)$(INSTALLSITEBIN)
INSTALLVENDORBIN = $(VENDORPREFIX)/bin
DESTINSTALLVENDORBIN = $(DESTDIR)$(INSTALLVENDORBIN)
INSTALLSCRIPT = $(PERLPREFIX)/bin
DESTINSTALLSCRIPT = $(DESTDIR)$(INSTALLSCRIPT)
INSTALLSITESCRIPT = $(SITEPREFIX)/bin
DESTINSTALLSITESCRIPT = $(DESTDIR)$(INSTALLSITESCRIPT)
INSTALLVENDORSCRIPT = $(VENDORPREFIX)/bin
DESTINSTALLVENDORSCRIPT = $(DESTDIR)$(INSTALLVENDORSCRIPT)
INSTALLMAN1DIR = $(PERLPREFIX)/man/man1
DESTINSTALLMAN1DIR = $(DESTDIR)$(INSTALLMAN1DIR)
INSTALLSITEMAN1DIR = $(SITEPREFIX)/man/man1
DESTINSTALLSITEMAN1DIR = $(DESTDIR)$(INSTALLSITEMAN1DIR)
INSTALLVENDORMAN1DIR = $(VENDORPREFIX)/man/man1
DESTINSTALLVENDORMAN1DIR = $(DESTDIR)$(INSTALLVENDORMAN1DIR)
INSTALLMAN3DIR = $(PERLPREFIX)/man/man3
DESTINSTALLMAN3DIR = $(DESTDIR)$(INSTALLMAN3DIR)
INSTALLSITEMAN3DIR = $(SITEPREFIX)/man/man3
DESTINSTALLSITEMAN3DIR = $(DESTDIR)$(INSTALLSITEMAN3DIR)
INSTALLVENDORMAN3DIR = $(VENDORPREFIX)/man/man3
DESTINSTALLVENDORMAN3DIR = $(DESTDIR)$(INSTALLVENDORMAN3DIR)
PERL_LIB = /local/perl/5.18.2/lib/5.18.2
PERL_ARCHLIB = /local/perl/5.18.2/lib/5.18.2/x86_64-linux
LIBPERL_A = libperl.a
FIRST_MAKEFILE = Makefile
MAKEFILE_OLD = Makefile.old
MAKE_APERL_FILE = Makefile.aperl
PERLMAINCC = $(CC)
PERL_INC = /local/perl/5.18.2/lib/5.18.2/x86_64-linux/CORE
PERL = /softs/local/perl/5.18.2/bin/perl
FULLPERL = /softs/local/perl/5.18.2/bin/perl
ABSPERL = $(PERL)
PERLRUN = $(PERL)
FULLPERLRUN = $(FULLPERL)
ABSPERLRUN = $(ABSPERL)
PERLRUNINST = $(PERLRUN) "-I$(INST_ARCHLIB)" "-I$(INST_LIB)"
FULLPERLRUNINST = $(FULLPERLRUN) "-I$(INST_ARCHLIB)" "-I$(INST_LIB)"
ABSPERLRUNINST = $(ABSPERLRUN) "-I$(INST_ARCHLIB)" "-I$(INST_LIB)"
PERL_CORE = 0
PERM_DIR = 755
PERM_RW = 644
PERM_RWX = 755

MAKEMAKER   = /local/perl/5.18.2/lib/5.18.2/ExtUtils/MakeMaker.pm
MM_VERSION  = 6.66
MM_REVISION = 66600

# FULLEXT = Pathname for extension directory (eg Foo/Bar/Oracle).
# BASEEXT = Basename part of FULLEXT. May be just equal FULLEXT. (eg Oracle)
# PARENT_NAME = NAME without BASEEXT and no trailing :: (eg Foo::Bar)
# DLBASE  = Basename part of dynamic library. May be just equal BASEEXT.
MAKE = make
FULLEXT = FEELnc
BASEEXT = FEELnc
PARENT_NAME = 
DLBASE = $(BASEEXT)
VERSION_FROM = lib/FEELnc.pm
OBJECT = 
LDFROM = $(OBJECT)
LINKTYPE = dynamic
BOOTDEP = 

# Handy lists of source code files:
XS_FILES = 
C_FILES  = 
O_FILES  = 
H_FILES  = 
MAN1PODS = 
MAN3PODS = lib/5.8.8/x86_64-linux/perllocal.pod \
	lib/Bio/SeqFeature/Extended.pm \
	lib/Bio/SeqFeature/Genic.pm \
	lib/Bio/SeqFeature/InterGenic.pm \
	lib/Bio/SeqFeature/Interaction.pm \
	lib/Bio/SeqFeature/InteractionCollection.pm \
	lib/Bio/SeqFeature/InteractionIterator.pm \
	lib/Bio/SeqFeature/LncRNAs_Factory.pm \
	lib/Bio/SeqFeature/database_part.pm \
	lib/Parser.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Extended.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Genic.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InterGenic.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Interaction.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InteractionCollection.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InteractionIterator.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/LncRNAs_Factory.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/database_part.pm \
	lib/site_perl/5.8.8/Parser.pm

# Where is the Config information that we are using/depend on
CONFIGDEP = $(PERL_ARCHLIB)$(DFSEP)Config.pm $(PERL_INC)$(DFSEP)config.h

# Where to build things
INST_LIBDIR      = $(INST_LIB)
INST_ARCHLIBDIR  = $(INST_ARCHLIB)

INST_AUTODIR     = $(INST_LIB)/auto/$(FULLEXT)
INST_ARCHAUTODIR = $(INST_ARCHLIB)/auto/$(FULLEXT)

INST_STATIC      = 
INST_DYNAMIC     = 
INST_BOOT        = 

# Extra linker info
EXPORT_LIST        = 
PERL_ARCHIVE       = 
PERL_ARCHIVE_AFTER = 


TO_INST_PM = lib/5.8.8/x86_64-linux/perllocal.pod \
	lib/Bio/SeqFeature/Empty.pm \
	lib/Bio/SeqFeature/Extended.pm \
	lib/Bio/SeqFeature/Genic.pm \
	lib/Bio/SeqFeature/InterGenic.pm \
	lib/Bio/SeqFeature/Interaction.pm \
	lib/Bio/SeqFeature/InteractionCollection.pm \
	lib/Bio/SeqFeature/InteractionIterator.pm \
	lib/Bio/SeqFeature/LncRNAs_Factory.pm \
	lib/Bio/SeqFeature/database_part.pm \
	lib/Cpat.pm \
	lib/ExtractFromFeature.pm \
	lib/ExtractFromHash.pm \
	lib/FEELnc.pm \
	lib/Filter.pm \
	lib/Intersect.pm \
	lib/Orf.pm \
	lib/Parser.pm \
	lib/StringUtils.pm \
	lib/Utils.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Empty.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Extended.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Genic.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InterGenic.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Interaction.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InteractionCollection.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InteractionIterator.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/LncRNAs_Factory.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/database_part.pm \
	lib/site_perl/5.8.8/Cpat.pm \
	lib/site_perl/5.8.8/ExtractFromFeature.pm \
	lib/site_perl/5.8.8/ExtractFromHash.pm \
	lib/site_perl/5.8.8/FEELnc.pm \
	lib/site_perl/5.8.8/Filter.pm \
	lib/site_perl/5.8.8/Intersect.pm \
	lib/site_perl/5.8.8/Orf.pm \
	lib/site_perl/5.8.8/Parser.pm \
	lib/site_perl/5.8.8/StringUtils.pm \
	lib/site_perl/5.8.8/Utils.pm \
	lib/site_perl/5.8.8/x86_64-linux/auto/FEELnc/.packlist

PM_TO_BLIB = lib/5.8.8/x86_64-linux/perllocal.pod \
	blib/lib/5.8.8/x86_64-linux/perllocal.pod \
	lib/site_perl/5.8.8/x86_64-linux/auto/FEELnc/.packlist \
	blib/lib/site_perl/5.8.8/x86_64-linux/auto/FEELnc/.packlist \
	lib/Bio/SeqFeature/Empty.pm \
	blib/lib/Bio/SeqFeature/Empty.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InteractionCollection.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/InteractionCollection.pm \
	lib/Intersect.pm \
	blib/lib/Intersect.pm \
	lib/Bio/SeqFeature/InteractionCollection.pm \
	blib/lib/Bio/SeqFeature/InteractionCollection.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Interaction.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/Interaction.pm \
	lib/ExtractFromFeature.pm \
	blib/lib/ExtractFromFeature.pm \
	lib/Bio/SeqFeature/Extended.pm \
	blib/lib/Bio/SeqFeature/Extended.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Extended.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/Extended.pm \
	lib/Bio/SeqFeature/Genic.pm \
	blib/lib/Bio/SeqFeature/Genic.pm \
	lib/site_perl/5.8.8/Filter.pm \
	blib/lib/site_perl/5.8.8/Filter.pm \
	lib/site_perl/5.8.8/Utils.pm \
	blib/lib/site_perl/5.8.8/Utils.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InterGenic.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/InterGenic.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Genic.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/Genic.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InteractionIterator.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/InteractionIterator.pm \
	lib/FEELnc.pm \
	blib/lib/FEELnc.pm \
	lib/site_perl/5.8.8/ExtractFromHash.pm \
	blib/lib/site_perl/5.8.8/ExtractFromHash.pm \
	lib/Parser.pm \
	blib/lib/Parser.pm \
	lib/site_perl/5.8.8/FEELnc.pm \
	blib/lib/site_perl/5.8.8/FEELnc.pm \
	lib/site_perl/5.8.8/ExtractFromFeature.pm \
	blib/lib/site_perl/5.8.8/ExtractFromFeature.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Empty.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/Empty.pm \
	lib/site_perl/5.8.8/Orf.pm \
	blib/lib/site_perl/5.8.8/Orf.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/database_part.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/database_part.pm \
	lib/Orf.pm \
	blib/lib/Orf.pm \
	lib/StringUtils.pm \
	blib/lib/StringUtils.pm \
	lib/Cpat.pm \
	blib/lib/Cpat.pm \
	lib/Bio/SeqFeature/LncRNAs_Factory.pm \
	blib/lib/Bio/SeqFeature/LncRNAs_Factory.pm \
	lib/Bio/SeqFeature/InteractionIterator.pm \
	blib/lib/Bio/SeqFeature/InteractionIterator.pm \
	lib/site_perl/5.8.8/Cpat.pm \
	blib/lib/site_perl/5.8.8/Cpat.pm \
	lib/site_perl/5.8.8/Parser.pm \
	blib/lib/site_perl/5.8.8/Parser.pm \
	lib/Bio/SeqFeature/database_part.pm \
	blib/lib/Bio/SeqFeature/database_part.pm \
	lib/ExtractFromHash.pm \
	blib/lib/ExtractFromHash.pm \
	lib/Utils.pm \
	blib/lib/Utils.pm \
	lib/site_perl/5.8.8/StringUtils.pm \
	blib/lib/site_perl/5.8.8/StringUtils.pm \
	lib/Filter.pm \
	blib/lib/Filter.pm \
	lib/site_perl/5.8.8/Intersect.pm \
	blib/lib/site_perl/5.8.8/Intersect.pm \
	lib/Bio/SeqFeature/Interaction.pm \
	blib/lib/Bio/SeqFeature/Interaction.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/LncRNAs_Factory.pm \
	blib/lib/site_perl/5.8.8/Bio/SeqFeature/LncRNAs_Factory.pm \
	lib/Bio/SeqFeature/InterGenic.pm \
	blib/lib/Bio/SeqFeature/InterGenic.pm


# --- MakeMaker platform_constants section:
MM_Unix_VERSION = 6.66
PERL_MALLOC_DEF = -DPERL_EXTMALLOC_DEF -Dmalloc=Perl_malloc -Dfree=Perl_mfree -Drealloc=Perl_realloc -Dcalloc=Perl_calloc


# --- MakeMaker tool_autosplit section:
# Usage: $(AUTOSPLITFILE) FileToSplit AutoDirToSplitInto
AUTOSPLITFILE = $(ABSPERLRUN)  -e 'use AutoSplit;  autosplit($$$$ARGV[0], $$$$ARGV[1], 0, 1, 1)' --



# --- MakeMaker tool_xsubpp section:


# --- MakeMaker tools_other section:
SHELL = /bin/sh
CHMOD = chmod
CP = cp
MV = mv
NOOP = $(TRUE)
NOECHO = @
RM_F = rm -f
RM_RF = rm -rf
TEST_F = test -f
TOUCH = touch
UMASK_NULL = umask 0
DEV_NULL = > /dev/null 2>&1
MKPATH = $(ABSPERLRUN) -MExtUtils::Command -e 'mkpath' --
EQUALIZE_TIMESTAMP = $(ABSPERLRUN) -MExtUtils::Command -e 'eqtime' --
FALSE = false
TRUE = true
ECHO = echo
ECHO_N = echo -n
UNINST = 0
VERBINST = 0
MOD_INSTALL = $(ABSPERLRUN) -MExtUtils::Install -e 'install([ from_to => {@ARGV}, verbose => '\''$(VERBINST)'\'', uninstall_shadows => '\''$(UNINST)'\'', dir_mode => '\''$(PERM_DIR)'\'' ]);' --
DOC_INSTALL = $(ABSPERLRUN) -MExtUtils::Command::MM -e 'perllocal_install' --
UNINSTALL = $(ABSPERLRUN) -MExtUtils::Command::MM -e 'uninstall' --
WARN_IF_OLD_PACKLIST = $(ABSPERLRUN) -MExtUtils::Command::MM -e 'warn_if_old_packlist' --
MACROSTART = 
MACROEND = 
USEMAKEFILE = -f
FIXIN = $(ABSPERLRUN) -MExtUtils::MY -e 'MY->fixin(shift)' --


# --- MakeMaker makemakerdflt section:
makemakerdflt : all
	$(NOECHO) $(NOOP)


# --- MakeMaker dist section:
TAR = tar
TARFLAGS = cvf
ZIP = zip
ZIPFLAGS = -r
COMPRESS = gzip --best
SUFFIX = .gz
SHAR = shar
PREOP = $(NOECHO) $(NOOP)
POSTOP = $(NOECHO) $(NOOP)
TO_UNIX = $(NOECHO) $(NOOP)
CI = ci -u
RCS_LABEL = rcs -Nv$(VERSION_SYM): -q
DIST_CP = best
DIST_DEFAULT = tardist
DISTNAME = FEELnc
DISTVNAME = FEELnc-0.1


# --- MakeMaker macro section:


# --- MakeMaker depend section:


# --- MakeMaker cflags section:


# --- MakeMaker const_loadlibs section:


# --- MakeMaker const_cccmd section:


# --- MakeMaker post_constants section:


# --- MakeMaker pasthru section:

PASTHRU = LIBPERL_A="$(LIBPERL_A)"\
	LINKTYPE="$(LINKTYPE)"\
	PREFIX="$(PREFIX)"


# --- MakeMaker special_targets section:
.SUFFIXES : .xs .c .C .cpp .i .s .cxx .cc $(OBJ_EXT)

.PHONY: all config static dynamic test linkext manifest blibdirs clean realclean disttest distdir



# --- MakeMaker c_o section:


# --- MakeMaker xs_c section:


# --- MakeMaker xs_o section:


# --- MakeMaker top_targets section:
all :: pure_all manifypods
	$(NOECHO) $(NOOP)


pure_all :: config pm_to_blib subdirs linkext
	$(NOECHO) $(NOOP)

subdirs :: $(MYEXTLIB)
	$(NOECHO) $(NOOP)

config :: $(FIRST_MAKEFILE) blibdirs
	$(NOECHO) $(NOOP)

help :
	perldoc ExtUtils::MakeMaker


# --- MakeMaker blibdirs section:
blibdirs : $(INST_LIBDIR)$(DFSEP).exists $(INST_ARCHLIB)$(DFSEP).exists $(INST_AUTODIR)$(DFSEP).exists $(INST_ARCHAUTODIR)$(DFSEP).exists $(INST_BIN)$(DFSEP).exists $(INST_SCRIPT)$(DFSEP).exists $(INST_MAN1DIR)$(DFSEP).exists $(INST_MAN3DIR)$(DFSEP).exists
	$(NOECHO) $(NOOP)

# Backwards compat with 6.18 through 6.25
blibdirs.ts : blibdirs
	$(NOECHO) $(NOOP)

$(INST_LIBDIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_LIBDIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_LIBDIR)
	$(NOECHO) $(TOUCH) $(INST_LIBDIR)$(DFSEP).exists

$(INST_ARCHLIB)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_ARCHLIB)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_ARCHLIB)
	$(NOECHO) $(TOUCH) $(INST_ARCHLIB)$(DFSEP).exists

$(INST_AUTODIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_AUTODIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_AUTODIR)
	$(NOECHO) $(TOUCH) $(INST_AUTODIR)$(DFSEP).exists

$(INST_ARCHAUTODIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_ARCHAUTODIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_ARCHAUTODIR)
	$(NOECHO) $(TOUCH) $(INST_ARCHAUTODIR)$(DFSEP).exists

$(INST_BIN)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_BIN)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_BIN)
	$(NOECHO) $(TOUCH) $(INST_BIN)$(DFSEP).exists

$(INST_SCRIPT)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_SCRIPT)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_SCRIPT)
	$(NOECHO) $(TOUCH) $(INST_SCRIPT)$(DFSEP).exists

$(INST_MAN1DIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_MAN1DIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_MAN1DIR)
	$(NOECHO) $(TOUCH) $(INST_MAN1DIR)$(DFSEP).exists

$(INST_MAN3DIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_MAN3DIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_MAN3DIR)
	$(NOECHO) $(TOUCH) $(INST_MAN3DIR)$(DFSEP).exists



# --- MakeMaker linkext section:

linkext :: $(LINKTYPE)
	$(NOECHO) $(NOOP)


# --- MakeMaker dlsyms section:


# --- MakeMaker dynamic section:

dynamic :: $(FIRST_MAKEFILE) $(INST_DYNAMIC) $(INST_BOOT)
	$(NOECHO) $(NOOP)


# --- MakeMaker dynamic_bs section:

BOOTSTRAP =


# --- MakeMaker dynamic_lib section:


# --- MakeMaker static section:

## $(INST_PM) has been moved to the all: target.
## It remains here for awhile to allow for old usage: "make static"
static :: $(FIRST_MAKEFILE) $(INST_STATIC)
	$(NOECHO) $(NOOP)


# --- MakeMaker static_lib section:


# --- MakeMaker manifypods section:

POD2MAN_EXE = $(PERLRUN) "-MExtUtils::Command::MM" -e pod2man "--"
POD2MAN = $(POD2MAN_EXE)


manifypods : pure_all  \
	lib/Bio/SeqFeature/Genic.pm \
	lib/site_perl/5.8.8/Parser.pm \
	lib/Bio/SeqFeature/database_part.pm \
	lib/Bio/SeqFeature/Interaction.pm \
	lib/Bio/SeqFeature/InterGenic.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/LncRNAs_Factory.pm \
	lib/Bio/SeqFeature/LncRNAs_Factory.pm \
	lib/Bio/SeqFeature/Extended.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Extended.pm \
	lib/Bio/SeqFeature/InteractionIterator.pm \
	lib/Parser.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/database_part.pm \
	lib/Bio/SeqFeature/InteractionCollection.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InteractionCollection.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/Interaction.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InterGenic.pm \
	lib/5.8.8/x86_64-linux/perllocal.pod \
	lib/site_perl/5.8.8/Bio/SeqFeature/Genic.pm \
	lib/site_perl/5.8.8/Bio/SeqFeature/InteractionIterator.pm
	$(NOECHO) $(POD2MAN) --section=3 --perm_rw=$(PERM_RW) \
	  lib/Bio/SeqFeature/Genic.pm $(INST_MAN3DIR)/Bio::SeqFeature::Genic.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Parser.pm $(INST_MAN3DIR)/site_perl::5.8.8::Parser.$(MAN3EXT) \
	  lib/Bio/SeqFeature/database_part.pm $(INST_MAN3DIR)/Bio::SeqFeature::database_part.$(MAN3EXT) \
	  lib/Bio/SeqFeature/Interaction.pm $(INST_MAN3DIR)/Bio::SeqFeature::Interaction.$(MAN3EXT) \
	  lib/Bio/SeqFeature/InterGenic.pm $(INST_MAN3DIR)/Bio::SeqFeature::InterGenic.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Bio/SeqFeature/LncRNAs_Factory.pm $(INST_MAN3DIR)/site_perl::5.8.8::Bio::SeqFeature::LncRNAs_Factory.$(MAN3EXT) \
	  lib/Bio/SeqFeature/LncRNAs_Factory.pm $(INST_MAN3DIR)/Bio::SeqFeature::LncRNAs_Factory.$(MAN3EXT) \
	  lib/Bio/SeqFeature/Extended.pm $(INST_MAN3DIR)/Bio::SeqFeature::Extended.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Bio/SeqFeature/Extended.pm $(INST_MAN3DIR)/site_perl::5.8.8::Bio::SeqFeature::Extended.$(MAN3EXT) \
	  lib/Bio/SeqFeature/InteractionIterator.pm $(INST_MAN3DIR)/Bio::SeqFeature::InteractionIterator.$(MAN3EXT) \
	  lib/Parser.pm $(INST_MAN3DIR)/Parser.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Bio/SeqFeature/database_part.pm $(INST_MAN3DIR)/site_perl::5.8.8::Bio::SeqFeature::database_part.$(MAN3EXT) \
	  lib/Bio/SeqFeature/InteractionCollection.pm $(INST_MAN3DIR)/Bio::SeqFeature::InteractionCollection.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Bio/SeqFeature/InteractionCollection.pm $(INST_MAN3DIR)/site_perl::5.8.8::Bio::SeqFeature::InteractionCollection.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Bio/SeqFeature/Interaction.pm $(INST_MAN3DIR)/site_perl::5.8.8::Bio::SeqFeature::Interaction.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Bio/SeqFeature/InterGenic.pm $(INST_MAN3DIR)/site_perl::5.8.8::Bio::SeqFeature::InterGenic.$(MAN3EXT) \
	  lib/5.8.8/x86_64-linux/perllocal.pod $(INST_MAN3DIR)/5.8.8::x86_64-linux::perllocal.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Bio/SeqFeature/Genic.pm $(INST_MAN3DIR)/site_perl::5.8.8::Bio::SeqFeature::Genic.$(MAN3EXT) \
	  lib/site_perl/5.8.8/Bio/SeqFeature/InteractionIterator.pm $(INST_MAN3DIR)/site_perl::5.8.8::Bio::SeqFeature::InteractionIterator.$(MAN3EXT) 




# --- MakeMaker processPL section:


# --- MakeMaker installbin section:


# --- MakeMaker subdirs section:

# none

# --- MakeMaker clean_subdirs section:
clean_subdirs :
	$(NOECHO) $(NOOP)


# --- MakeMaker clean section:

# Delete temporary files but do not touch installed files. We don't delete
# the Makefile here so a later make realclean still has a makefile to use.

clean :: clean_subdirs
	- $(RM_F) \
	  lib$(BASEEXT).def $(INST_ARCHAUTODIR)/extralibs.ld \
	  pm_to_blib.ts mon.out \
	  $(BASEEXT).exp core.[0-9][0-9][0-9][0-9] \
	  *perl.core $(INST_ARCHAUTODIR)/extralibs.all \
	  *$(LIB_EXT) core.[0-9] \
	  $(MAKE_APERL_FILE) core.*perl.*.? \
	  blibdirs.ts core.[0-9][0-9] \
	  pm_to_blib MYMETA.yml \
	  $(BASEEXT).def perl$(EXE_EXT) \
	  tmon.out $(BOOTSTRAP) \
	  core.[0-9][0-9][0-9] core.[0-9][0-9][0-9][0-9][0-9] \
	  MYMETA.json perlmain.c \
	  *$(OBJ_EXT) so_locations \
	  perl.exe $(BASEEXT).bso \
	  $(BASEEXT).x perl \
	  core 
	- $(RM_RF) \
	  blib 
	- $(MV) $(FIRST_MAKEFILE) $(MAKEFILE_OLD) $(DEV_NULL)


# --- MakeMaker realclean_subdirs section:
realclean_subdirs :
	$(NOECHO) $(NOOP)


# --- MakeMaker realclean section:
# Delete temporary files (via clean) and also delete dist files
realclean purge ::  clean realclean_subdirs
	- $(RM_F) \
	  $(FIRST_MAKEFILE) $(MAKEFILE_OLD) 
	- $(RM_RF) \
	  $(DISTVNAME) 


# --- MakeMaker metafile section:
metafile : create_distdir
	$(NOECHO) $(ECHO) Generating META.yml
	$(NOECHO) $(ECHO) '---' > META_new.yml
	$(NOECHO) $(ECHO) 'abstract: unknown' >> META_new.yml
	$(NOECHO) $(ECHO) 'author:' >> META_new.yml
	$(NOECHO) $(ECHO) '  - '\'' Fabrice Legeai ; Audrey David ; Thomas Derrien'\''' >> META_new.yml
	$(NOECHO) $(ECHO) 'build_requires:' >> META_new.yml
	$(NOECHO) $(ECHO) '  ExtUtils::MakeMaker: 0' >> META_new.yml
	$(NOECHO) $(ECHO) 'configure_requires:' >> META_new.yml
	$(NOECHO) $(ECHO) '  ExtUtils::MakeMaker: 0' >> META_new.yml
	$(NOECHO) $(ECHO) 'dynamic_config: 1' >> META_new.yml
	$(NOECHO) $(ECHO) 'generated_by: '\''ExtUtils::MakeMaker version 6.66, CPAN::Meta::Converter version 2.120921'\''' >> META_new.yml
	$(NOECHO) $(ECHO) 'license: unknown' >> META_new.yml
	$(NOECHO) $(ECHO) 'meta-spec:' >> META_new.yml
	$(NOECHO) $(ECHO) '  url: http://module-build.sourceforge.net/META-spec-v1.4.html' >> META_new.yml
	$(NOECHO) $(ECHO) '  version: 1.4' >> META_new.yml
	$(NOECHO) $(ECHO) 'name: FEELnc' >> META_new.yml
	$(NOECHO) $(ECHO) 'no_index:' >> META_new.yml
	$(NOECHO) $(ECHO) '  directory:' >> META_new.yml
	$(NOECHO) $(ECHO) '    - t' >> META_new.yml
	$(NOECHO) $(ECHO) '    - inc' >> META_new.yml
	$(NOECHO) $(ECHO) 'requires:' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::SeqFeature::Generic: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::Tools::GFF: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Parallel::ForkManager: 1.07' >> META_new.yml
	$(NOECHO) $(ECHO) 'version: 0.1' >> META_new.yml
	-$(NOECHO) $(MV) META_new.yml $(DISTVNAME)/META.yml
	$(NOECHO) $(ECHO) Generating META.json
	$(NOECHO) $(ECHO) '{' > META_new.json
	$(NOECHO) $(ECHO) '   "abstract" : "unknown",' >> META_new.json
	$(NOECHO) $(ECHO) '   "author" : [' >> META_new.json
	$(NOECHO) $(ECHO) '      " Fabrice Legeai ; Audrey David ; Thomas Derrien"' >> META_new.json
	$(NOECHO) $(ECHO) '   ],' >> META_new.json
	$(NOECHO) $(ECHO) '   "dynamic_config" : 1,' >> META_new.json
	$(NOECHO) $(ECHO) '   "generated_by" : "ExtUtils::MakeMaker version 6.66, CPAN::Meta::Converter version 2.120921",' >> META_new.json
	$(NOECHO) $(ECHO) '   "license" : [' >> META_new.json
	$(NOECHO) $(ECHO) '      "unknown"' >> META_new.json
	$(NOECHO) $(ECHO) '   ],' >> META_new.json
	$(NOECHO) $(ECHO) '   "meta-spec" : {' >> META_new.json
	$(NOECHO) $(ECHO) '      "url" : "http://search.cpan.org/perldoc?CPAN::Meta::Spec",' >> META_new.json
	$(NOECHO) $(ECHO) '      "version" : "2"' >> META_new.json
	$(NOECHO) $(ECHO) '   },' >> META_new.json
	$(NOECHO) $(ECHO) '   "name" : "FEELnc",' >> META_new.json
	$(NOECHO) $(ECHO) '   "no_index" : {' >> META_new.json
	$(NOECHO) $(ECHO) '      "directory" : [' >> META_new.json
	$(NOECHO) $(ECHO) '         "t",' >> META_new.json
	$(NOECHO) $(ECHO) '         "inc"' >> META_new.json
	$(NOECHO) $(ECHO) '      ]' >> META_new.json
	$(NOECHO) $(ECHO) '   },' >> META_new.json
	$(NOECHO) $(ECHO) '   "prereqs" : {' >> META_new.json
	$(NOECHO) $(ECHO) '      "build" : {' >> META_new.json
	$(NOECHO) $(ECHO) '         "requires" : {' >> META_new.json
	$(NOECHO) $(ECHO) '            "ExtUtils::MakeMaker" : "0"' >> META_new.json
	$(NOECHO) $(ECHO) '         }' >> META_new.json
	$(NOECHO) $(ECHO) '      },' >> META_new.json
	$(NOECHO) $(ECHO) '      "configure" : {' >> META_new.json
	$(NOECHO) $(ECHO) '         "requires" : {' >> META_new.json
	$(NOECHO) $(ECHO) '            "ExtUtils::MakeMaker" : "0"' >> META_new.json
	$(NOECHO) $(ECHO) '         }' >> META_new.json
	$(NOECHO) $(ECHO) '      },' >> META_new.json
	$(NOECHO) $(ECHO) '      "runtime" : {' >> META_new.json
	$(NOECHO) $(ECHO) '         "requires" : {' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::SeqFeature::Generic" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::Tools::GFF" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Parallel::ForkManager" : "1.07"' >> META_new.json
	$(NOECHO) $(ECHO) '         }' >> META_new.json
	$(NOECHO) $(ECHO) '      }' >> META_new.json
	$(NOECHO) $(ECHO) '   },' >> META_new.json
	$(NOECHO) $(ECHO) '   "release_status" : "stable",' >> META_new.json
	$(NOECHO) $(ECHO) '   "version" : "0.1"' >> META_new.json
	$(NOECHO) $(ECHO) '}' >> META_new.json
	-$(NOECHO) $(MV) META_new.json $(DISTVNAME)/META.json


# --- MakeMaker signature section:
signature :
	cpansign -s


# --- MakeMaker dist_basics section:
distclean :: realclean distcheck
	$(NOECHO) $(NOOP)

distcheck :
	$(PERLRUN) "-MExtUtils::Manifest=fullcheck" -e fullcheck

skipcheck :
	$(PERLRUN) "-MExtUtils::Manifest=skipcheck" -e skipcheck

manifest :
	$(PERLRUN) "-MExtUtils::Manifest=mkmanifest" -e mkmanifest

veryclean : realclean
	$(RM_F) *~ */*~ *.orig */*.orig *.bak */*.bak *.old */*.old 



# --- MakeMaker dist_core section:

dist : $(DIST_DEFAULT) $(FIRST_MAKEFILE)
	$(NOECHO) $(ABSPERLRUN) -l -e 'print '\''Warning: Makefile possibly out of date with $(VERSION_FROM)'\''' \
	  -e '    if -e '\''$(VERSION_FROM)'\'' and -M '\''$(VERSION_FROM)'\'' < -M '\''$(FIRST_MAKEFILE)'\'';' --

tardist : $(DISTVNAME).tar$(SUFFIX)
	$(NOECHO) $(NOOP)

uutardist : $(DISTVNAME).tar$(SUFFIX)
	uuencode $(DISTVNAME).tar$(SUFFIX) $(DISTVNAME).tar$(SUFFIX) > $(DISTVNAME).tar$(SUFFIX)_uu

$(DISTVNAME).tar$(SUFFIX) : distdir
	$(PREOP)
	$(TO_UNIX)
	$(TAR) $(TARFLAGS) $(DISTVNAME).tar $(DISTVNAME)
	$(RM_RF) $(DISTVNAME)
	$(COMPRESS) $(DISTVNAME).tar
	$(POSTOP)

zipdist : $(DISTVNAME).zip
	$(NOECHO) $(NOOP)

$(DISTVNAME).zip : distdir
	$(PREOP)
	$(ZIP) $(ZIPFLAGS) $(DISTVNAME).zip $(DISTVNAME)
	$(RM_RF) $(DISTVNAME)
	$(POSTOP)

shdist : distdir
	$(PREOP)
	$(SHAR) $(DISTVNAME) > $(DISTVNAME).shar
	$(RM_RF) $(DISTVNAME)
	$(POSTOP)


# --- MakeMaker distdir section:
create_distdir :
	$(RM_RF) $(DISTVNAME)
	$(PERLRUN) "-MExtUtils::Manifest=manicopy,maniread" \
		-e "manicopy(maniread(),'$(DISTVNAME)', '$(DIST_CP)');"

distdir : create_distdir distmeta 
	$(NOECHO) $(NOOP)



# --- MakeMaker dist_test section:
disttest : distdir
	cd $(DISTVNAME) && $(ABSPERLRUN) Makefile.PL "PREFIX=/home/genouest/umr6061/recomgen/tderrien/bin/perl/FEELnc_v0.1"
	cd $(DISTVNAME) && $(MAKE) $(PASTHRU)
	cd $(DISTVNAME) && $(MAKE) test $(PASTHRU)



# --- MakeMaker dist_ci section:

ci :
	$(PERLRUN) "-MExtUtils::Manifest=maniread" \
	  -e "@all = keys %{ maniread() };" \
	  -e "print(qq{Executing $(CI) @all\n}); system(qq{$(CI) @all});" \
	  -e "print(qq{Executing $(RCS_LABEL) ...\n}); system(qq{$(RCS_LABEL) @all});"


# --- MakeMaker distmeta section:
distmeta : create_distdir metafile
	$(NOECHO) cd $(DISTVNAME) && $(ABSPERLRUN) -MExtUtils::Manifest=maniadd -e 'exit unless -e q{META.yml};' \
	  -e 'eval { maniadd({q{META.yml} => q{Module YAML meta-data (added by MakeMaker)}}) }' \
	  -e '    or print "Could not add META.yml to MANIFEST: $$$${'\''@'\''}\n"' --
	$(NOECHO) cd $(DISTVNAME) && $(ABSPERLRUN) -MExtUtils::Manifest=maniadd -e 'exit unless -f q{META.json};' \
	  -e 'eval { maniadd({q{META.json} => q{Module JSON meta-data (added by MakeMaker)}}) }' \
	  -e '    or print "Could not add META.json to MANIFEST: $$$${'\''@'\''}\n"' --



# --- MakeMaker distsignature section:
distsignature : create_distdir
	$(NOECHO) cd $(DISTVNAME) && $(ABSPERLRUN) -MExtUtils::Manifest=maniadd -e 'eval { maniadd({q{SIGNATURE} => q{Public-key signature (added by MakeMaker)}}) } ' \
	  -e '    or print "Could not add SIGNATURE to MANIFEST: $$$${'\''@'\''}\n"' --
	$(NOECHO) cd $(DISTVNAME) && $(TOUCH) SIGNATURE
	cd $(DISTVNAME) && cpansign -s



# --- MakeMaker install section:

install :: pure_install doc_install
	$(NOECHO) $(NOOP)

install_perl :: pure_perl_install doc_perl_install
	$(NOECHO) $(NOOP)

install_site :: pure_site_install doc_site_install
	$(NOECHO) $(NOOP)

install_vendor :: pure_vendor_install doc_vendor_install
	$(NOECHO) $(NOOP)

pure_install :: pure_$(INSTALLDIRS)_install
	$(NOECHO) $(NOOP)

doc_install :: doc_$(INSTALLDIRS)_install
	$(NOECHO) $(NOOP)

pure__install : pure_site_install
	$(NOECHO) $(ECHO) INSTALLDIRS not defined, defaulting to INSTALLDIRS=site

doc__install : doc_site_install
	$(NOECHO) $(ECHO) INSTALLDIRS not defined, defaulting to INSTALLDIRS=site

pure_perl_install :: all
	$(NOECHO) $(MOD_INSTALL) \
		read $(PERL_ARCHLIB)/auto/$(FULLEXT)/.packlist \
		write $(DESTINSTALLARCHLIB)/auto/$(FULLEXT)/.packlist \
		$(INST_LIB) $(DESTINSTALLPRIVLIB) \
		$(INST_ARCHLIB) $(DESTINSTALLARCHLIB) \
		$(INST_BIN) $(DESTINSTALLBIN) \
		$(INST_SCRIPT) $(DESTINSTALLSCRIPT) \
		$(INST_MAN1DIR) $(DESTINSTALLMAN1DIR) \
		$(INST_MAN3DIR) $(DESTINSTALLMAN3DIR)
	$(NOECHO) $(WARN_IF_OLD_PACKLIST) \
		$(SITEARCHEXP)/auto/$(FULLEXT)


pure_site_install :: all
	$(NOECHO) $(MOD_INSTALL) \
		read $(SITEARCHEXP)/auto/$(FULLEXT)/.packlist \
		write $(DESTINSTALLSITEARCH)/auto/$(FULLEXT)/.packlist \
		$(INST_LIB) $(DESTINSTALLSITELIB) \
		$(INST_ARCHLIB) $(DESTINSTALLSITEARCH) \
		$(INST_BIN) $(DESTINSTALLSITEBIN) \
		$(INST_SCRIPT) $(DESTINSTALLSITESCRIPT) \
		$(INST_MAN1DIR) $(DESTINSTALLSITEMAN1DIR) \
		$(INST_MAN3DIR) $(DESTINSTALLSITEMAN3DIR)
	$(NOECHO) $(WARN_IF_OLD_PACKLIST) \
		$(PERL_ARCHLIB)/auto/$(FULLEXT)

pure_vendor_install :: all
	$(NOECHO) $(MOD_INSTALL) \
		read $(VENDORARCHEXP)/auto/$(FULLEXT)/.packlist \
		write $(DESTINSTALLVENDORARCH)/auto/$(FULLEXT)/.packlist \
		$(INST_LIB) $(DESTINSTALLVENDORLIB) \
		$(INST_ARCHLIB) $(DESTINSTALLVENDORARCH) \
		$(INST_BIN) $(DESTINSTALLVENDORBIN) \
		$(INST_SCRIPT) $(DESTINSTALLVENDORSCRIPT) \
		$(INST_MAN1DIR) $(DESTINSTALLVENDORMAN1DIR) \
		$(INST_MAN3DIR) $(DESTINSTALLVENDORMAN3DIR)

doc_perl_install :: all
	$(NOECHO) $(ECHO) Appending installation info to $(DESTINSTALLARCHLIB)/perllocal.pod
	-$(NOECHO) $(MKPATH) $(DESTINSTALLARCHLIB)
	-$(NOECHO) $(DOC_INSTALL) \
		"Module" "$(NAME)" \
		"installed into" "$(INSTALLPRIVLIB)" \
		LINKTYPE "$(LINKTYPE)" \
		VERSION "$(VERSION)" \
		EXE_FILES "$(EXE_FILES)" \
		>> $(DESTINSTALLARCHLIB)/perllocal.pod

doc_site_install :: all
	$(NOECHO) $(ECHO) Appending installation info to $(DESTINSTALLARCHLIB)/perllocal.pod
	-$(NOECHO) $(MKPATH) $(DESTINSTALLARCHLIB)
	-$(NOECHO) $(DOC_INSTALL) \
		"Module" "$(NAME)" \
		"installed into" "$(INSTALLSITELIB)" \
		LINKTYPE "$(LINKTYPE)" \
		VERSION "$(VERSION)" \
		EXE_FILES "$(EXE_FILES)" \
		>> $(DESTINSTALLARCHLIB)/perllocal.pod

doc_vendor_install :: all
	$(NOECHO) $(ECHO) Appending installation info to $(DESTINSTALLARCHLIB)/perllocal.pod
	-$(NOECHO) $(MKPATH) $(DESTINSTALLARCHLIB)
	-$(NOECHO) $(DOC_INSTALL) \
		"Module" "$(NAME)" \
		"installed into" "$(INSTALLVENDORLIB)" \
		LINKTYPE "$(LINKTYPE)" \
		VERSION "$(VERSION)" \
		EXE_FILES "$(EXE_FILES)" \
		>> $(DESTINSTALLARCHLIB)/perllocal.pod


uninstall :: uninstall_from_$(INSTALLDIRS)dirs
	$(NOECHO) $(NOOP)

uninstall_from_perldirs ::
	$(NOECHO) $(UNINSTALL) $(PERL_ARCHLIB)/auto/$(FULLEXT)/.packlist

uninstall_from_sitedirs ::
	$(NOECHO) $(UNINSTALL) $(SITEARCHEXP)/auto/$(FULLEXT)/.packlist

uninstall_from_vendordirs ::
	$(NOECHO) $(UNINSTALL) $(VENDORARCHEXP)/auto/$(FULLEXT)/.packlist


# --- MakeMaker force section:
# Phony target to force checking subdirectories.
FORCE :
	$(NOECHO) $(NOOP)


# --- MakeMaker perldepend section:


# --- MakeMaker makefile section:
# We take a very conservative approach here, but it's worth it.
# We move Makefile to Makefile.old here to avoid gnu make looping.
$(FIRST_MAKEFILE) : Makefile.PL $(CONFIGDEP)
	$(NOECHO) $(ECHO) "Makefile out-of-date with respect to $?"
	$(NOECHO) $(ECHO) "Cleaning current config before rebuilding Makefile..."
	-$(NOECHO) $(RM_F) $(MAKEFILE_OLD)
	-$(NOECHO) $(MV)   $(FIRST_MAKEFILE) $(MAKEFILE_OLD)
	- $(MAKE) $(USEMAKEFILE) $(MAKEFILE_OLD) clean $(DEV_NULL)
	$(PERLRUN) Makefile.PL "PREFIX=/home/genouest/umr6061/recomgen/tderrien/bin/perl/FEELnc_v0.1"
	$(NOECHO) $(ECHO) "==> Your Makefile has been rebuilt. <=="
	$(NOECHO) $(ECHO) "==> Please rerun the $(MAKE) command.  <=="
	$(FALSE)



# --- MakeMaker staticmake section:

# --- MakeMaker makeaperl section ---
MAP_TARGET    = perl
FULLPERL      = /softs/local/perl/5.18.2/bin/perl

$(MAP_TARGET) :: static $(MAKE_APERL_FILE)
	$(MAKE) $(USEMAKEFILE) $(MAKE_APERL_FILE) $@

$(MAKE_APERL_FILE) : $(FIRST_MAKEFILE) pm_to_blib
	$(NOECHO) $(ECHO) Writing \"$(MAKE_APERL_FILE)\" for this $(MAP_TARGET)
	$(NOECHO) $(PERLRUNINST) \
		Makefile.PL DIR= \
		MAKEFILE=$(MAKE_APERL_FILE) LINKTYPE=static \
		MAKEAPERL=1 NORECURS=1 CCCDLFLAGS= \
		PREFIX=/home/genouest/umr6061/recomgen/tderrien/bin/perl/FEELnc_v0.1


# --- MakeMaker test section:

TEST_VERBOSE=0
TEST_TYPE=test_$(LINKTYPE)
TEST_FILE = test.pl
TEST_FILES = 
TESTDB_SW = -d

testdb :: testdb_$(LINKTYPE)

test :: $(TEST_TYPE) subdirs-test

subdirs-test ::
	$(NOECHO) $(NOOP)

	$(NOECHO) $(ECHO) 'No tests defined for $(NAME) extension.'

test_dynamic :: pure_all

testdb_dynamic :: pure_all
	PERL_DL_NONLAZY=1 $(FULLPERLRUN) $(TESTDB_SW) "-I$(INST_LIB)" "-I$(INST_ARCHLIB)" $(TEST_FILE)

test_ : test_dynamic

test_static :: test_dynamic
testdb_static :: testdb_dynamic


# --- MakeMaker ppd section:
# Creates a PPD (Perl Package Description) for a binary distribution.
ppd :
	$(NOECHO) $(ECHO) '<SOFTPKG NAME="$(DISTNAME)" VERSION="$(VERSION)">' > $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '    <ABSTRACT></ABSTRACT>' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '    <AUTHOR> Fabrice Legeai ; Audrey David ; Thomas Derrien</AUTHOR>' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '    <IMPLEMENTATION>' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::SeqFeature::Generic" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::Tools::GFF" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE VERSION="1.07" NAME="Parallel::ForkManager" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <ARCHITECTURE NAME="x86_64-linux-5.18" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <CODEBASE HREF="" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '    </IMPLEMENTATION>' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '</SOFTPKG>' >> $(DISTNAME).ppd


# --- MakeMaker pm_to_blib section:

pm_to_blib : $(FIRST_MAKEFILE) $(TO_INST_PM)
	$(NOECHO) $(ABSPERLRUN) -MExtUtils::Install -e 'pm_to_blib({@ARGV}, '\''$(INST_LIB)/auto'\'', q[$(PM_FILTER)], '\''$(PERM_DIR)'\'')' -- \
	  lib/5.8.8/x86_64-linux/perllocal.pod blib/lib/5.8.8/x86_64-linux/perllocal.pod \
	  lib/site_perl/5.8.8/x86_64-linux/auto/FEELnc/.packlist blib/lib/site_perl/5.8.8/x86_64-linux/auto/FEELnc/.packlist \
	  lib/Bio/SeqFeature/Empty.pm blib/lib/Bio/SeqFeature/Empty.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/InteractionCollection.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/InteractionCollection.pm \
	  lib/Intersect.pm blib/lib/Intersect.pm \
	  lib/Bio/SeqFeature/InteractionCollection.pm blib/lib/Bio/SeqFeature/InteractionCollection.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/Interaction.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/Interaction.pm \
	  lib/ExtractFromFeature.pm blib/lib/ExtractFromFeature.pm \
	  lib/Bio/SeqFeature/Extended.pm blib/lib/Bio/SeqFeature/Extended.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/Extended.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/Extended.pm \
	  lib/Bio/SeqFeature/Genic.pm blib/lib/Bio/SeqFeature/Genic.pm \
	  lib/site_perl/5.8.8/Filter.pm blib/lib/site_perl/5.8.8/Filter.pm \
	  lib/site_perl/5.8.8/Utils.pm blib/lib/site_perl/5.8.8/Utils.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/InterGenic.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/InterGenic.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/Genic.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/Genic.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/InteractionIterator.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/InteractionIterator.pm \
	  lib/FEELnc.pm blib/lib/FEELnc.pm \
	  lib/site_perl/5.8.8/ExtractFromHash.pm blib/lib/site_perl/5.8.8/ExtractFromHash.pm \
	  lib/Parser.pm blib/lib/Parser.pm \
	  lib/site_perl/5.8.8/FEELnc.pm blib/lib/site_perl/5.8.8/FEELnc.pm \
	  lib/site_perl/5.8.8/ExtractFromFeature.pm blib/lib/site_perl/5.8.8/ExtractFromFeature.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/Empty.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/Empty.pm \
	  lib/site_perl/5.8.8/Orf.pm blib/lib/site_perl/5.8.8/Orf.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/database_part.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/database_part.pm \
	  lib/Orf.pm blib/lib/Orf.pm \
	  lib/StringUtils.pm blib/lib/StringUtils.pm \
	  lib/Cpat.pm blib/lib/Cpat.pm \
	  lib/Bio/SeqFeature/LncRNAs_Factory.pm blib/lib/Bio/SeqFeature/LncRNAs_Factory.pm \
	  lib/Bio/SeqFeature/InteractionIterator.pm blib/lib/Bio/SeqFeature/InteractionIterator.pm \
	  lib/site_perl/5.8.8/Cpat.pm blib/lib/site_perl/5.8.8/Cpat.pm \
	  lib/site_perl/5.8.8/Parser.pm blib/lib/site_perl/5.8.8/Parser.pm \
	  lib/Bio/SeqFeature/database_part.pm blib/lib/Bio/SeqFeature/database_part.pm \
	  lib/ExtractFromHash.pm blib/lib/ExtractFromHash.pm \
	  lib/Utils.pm blib/lib/Utils.pm \
	  lib/site_perl/5.8.8/StringUtils.pm blib/lib/site_perl/5.8.8/StringUtils.pm \
	  lib/Filter.pm blib/lib/Filter.pm \
	  lib/site_perl/5.8.8/Intersect.pm blib/lib/site_perl/5.8.8/Intersect.pm \
	  lib/Bio/SeqFeature/Interaction.pm blib/lib/Bio/SeqFeature/Interaction.pm \
	  lib/site_perl/5.8.8/Bio/SeqFeature/LncRNAs_Factory.pm blib/lib/site_perl/5.8.8/Bio/SeqFeature/LncRNAs_Factory.pm \
	  lib/Bio/SeqFeature/InterGenic.pm blib/lib/Bio/SeqFeature/InterGenic.pm 
	$(NOECHO) $(TOUCH) pm_to_blib


# --- MakeMaker selfdocument section:


# --- MakeMaker postamble section:


# End.
