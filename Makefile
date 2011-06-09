# these values filled in by    yorick -batch make.i
Y_MAKEDIR=/home/brujo/yorick-2.2/relocate
Y_EXE=/home/brujo/yorick-2.2/relocate/bin/yorick
Y_EXE_PKGS=

Y_EXE_HOME=/home/brujo/yorick-2.2/relocate
Y_EXE_SITE=/home/brujo/yorick-2.2/relocate

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)
# ------------------------------------------------ macros for this package

PKG_NAME=yoga_ao
PKG_I= yoga_aolib.i

OBJS= yoga_phase.o yoga_turbu.o yoga_target.cuo yoga_target.o yoga_phase.cuo yoga_wfs.cuo yoga_wfs.o yoga_ao.o #yoga_turbu.cuo idx_test.o yoga_svd.o yoga_svd.cuo 

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS= -L/usr/local/cuda/lib64 -lcudart -lcublas -lcufft -lstdc++ -lcurand
#-lcurand -lcula -L/usr/local/cula/lib 
PKG_DEPLIBS += ../cudpp/lib/libcudpp_x86_64.a
PKG_DEPLIBS += $(Y_HOME)/lib/libyoga.so #-L../ -lyoga

# set compiler (or rarely loader) flags specific to this package

PKG_CFLAGS = -I.. -I../cudpp/include
#PKG_CPPFLAGS = -I.. -I../cudpp/
#-I/usr/local/cula/include #-DDEBUG
PKG_LDFLAGS= 

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
PKG_I_EXTRA = yoga_ao.i yoga_ao_utils.i yoga_turbu.i iterkolmo.i yoga_ao_ystruct.i

# -------------------------------- standard targets and rules (in Makepkg)

#--compiler-options -fpermissive 

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package
# color for printf
#31=red, 32=green, 33=yellow,34=blue, 35=pink, 36=cyan, 37=white

GCC ?= $(CC)
NVCC = /usr/local/cuda/bin/nvcc
NVCCFLAGS = -m64 -Xcompiler -fPIC -c --compiler-options -fpermissive -I. -I.. -I../cudpp/include # --ptxas-options=-v #-I/usr/local/cula/include #-arch sm_13
ifneq ($(CUDPP_PATH),)
 NVCCFLAGS +=-I$(CUDPP_PATH)/include 
endif

%.cuo: %.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(NVCCFLAGS) -o $@ $< 

%.o: %.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) $(PKG_CFLAGS) -c -o $@ $<

$(PKG_NAME): $(PKG_NAME).o
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(PKG_NAME)
	@$(GCC) $(PKG_CPPFLAGS) $(PKG_DEPLIBS) -o $@ $^ 

configure:
	yorick -batch make.i

# -------------------------------------------------------- end of Makefile
