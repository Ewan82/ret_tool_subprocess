#======================================================================
#
#   file:               Makefile 
#
#   purpose:            environment to build, run, and test the
#                       Retrieval System developped within the framework of
#                       the Sentinel Synergy Study (S^3)
#
#   created:            2018/03
#
#   author:             Michael.Vossbeck(at)Inversion-Lab.com
#
#======================================================================
.DEFAULT_GOAL := help

#-----------------------------
#     command shortcuts
#
RM := rm -f


#-----------------------------
#     d i r e c t o r y   s e t t i n g s
#
BINDIR       := bin
DRVDIR       := driver
SRCDIR       := src
LIBLBFGSBDIR := lbfgsb
LIBEISDIR    := eispack


#-----------------------------
#     retrieval experiment configuration
#
RS_SETUP_TOOL   := $(BINDIR)/rs_setup.py
RS_PRE_TOOL     := $(BINDIR)/rs_pre.py
PRE_MODE        := pre_general
RETRCONFIG_FILE := retrconfig.nc
RETRPRIOR_FILE  := retrprior.nc
RETRMODEL_FILE  := retrmodel.nc
RETRCTL_FILE    := retrctl.nml
SITE_FILE       := site.nml
BOUNDS_FILE     := retrctlvecbounds.nml
INPUT_DIR       := input
EXPV            := V1.4.3
TIMSTART        := 20170501
TIMEND          := 20170531
TGTSCHEDULE_FILE := $(INPUT_DIR)/target_schedule_2017.txt
MISSION_LST     = S1A S1B S2A S2B


#-----------------------------
#     path to source directories
#
SRCDIR     := src
SEMID_SDIR := $(SRCDIR)/semid
SENSE_SDIR := $(SRCDIR)/sense
MDIR       := $(SRCDIR)/simulator
RETR_SDIR  := $(SRCDIR)/retrieval
DVSRCDIR   := $(SRCDIR)/derivatives


#-- various utility modules/routines/...
SUPPORT_SDIR  := $(SRCDIR)/support
SUPPORT_OBJS += $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$(shell ls $(SUPPORT_SDIR)/*.f90))))


#-----------------------------
# - all generated library files will reside in $(LIBDIR)
# - all generated object files will reside in $(OBJDIR)
#
LIBDIR := lib
SCRATCH_DIRS += $(LIBDIR)
OBJDIR := objs
ifndef OBJDIR
  $(error OBJDIR *must* be defined in Makefile (use '.' if object codes shall be generated in current working directory))
endif
ifndef LIBDIR
  $(error LIBDIR *must* be defined in Makefile (use '.' if libraries shall be generated in current working directory))
endif

SCRATCH_DIRS += $(OBJDIR)

#-----------------------------
# - where all generated libraries will reside
#
LIBSEMID   = $(LIBDIR)/libsemidis.a
LIBSENSE   = $(LIBDIR)/libsense.a
LIBMODEL   = $(LIBDIR)/libmodel.a
LIBRETR    = $(LIBDIR)/libretr.a
LIBDERIV   = $(LIBDIR)/libderiv.a
LIBSUPPORT = $(LIBDIR)/libsupport.a

#-----------------------------
#
#     external libraries
#
LIBEIS    = $(LIBDIR)/libeis.a
LIBLBFGSB = $(LIBDIR)/liblbfgsb-$(FC).a

#-----------------------------
#
#     SEMIDISCRETE (Fortran) sources
#
#
SEMID_ASRCS = \
$(SEMID_SDIR)/mo_sensimul_s2.f90 \
$(SEMID_SDIR)/price_soil_ftn.f90 \
$(SEMID_SDIR)/prospect_interface_ftn.f90 \
$(SEMID_SDIR)/prospect.redux.f \
$(SEMID_SDIR)/semiDbrf.f \
$(SEMID_SDIR)/semiDbrfe.f \
$(SEMID_SDIR)/semiDtools.f \
$(SEMID_SDIR)/semiD_ftn.f90
SEMID_OBJS = $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f,%.o,$(filter %.f,$(SEMID_ASRCS)))))
SEMID_OBJS += $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$(filter %.f90,$(SEMID_ASRCS)))))


#-----------------------------
#
#     SENSE (Fortran) sources
#
#
SENSE_ASRCS = \
$(SENSE_SDIR)/canopy.f90 \
$(SENSE_SDIR)/core.f90 \
$(SENSE_SDIR)/dobson85.f90 \
$(SENSE_SDIR)/ground.f90 \
$(SENSE_SDIR)/mo_sensimul_s1.f90 \
$(SENSE_SDIR)/Oh92.f90 \
$(SENSE_SDIR)/scatterer.f90 \
$(SENSE_SDIR)/soil.f90 \
$(SENSE_SDIR)/util.f90 \
$(SENSE_SDIR)/sense.f90
SENSE_OBJS = $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$(filter %.f90,$(SENSE_ASRCS)))))

MODEL_ASRCS = \
$(MDIR)/mo_sensimul.f90 \
$(MDIR)/simulate.f90 \
$(MDIR)/mapping.f90 \
$(MDIR)/model.f90 \
$(MDIR)/target_operator.f90

RETR_ASRCS = \
$(RETR_SDIR)/cost.f90 \
$(RETR_SDIR)/misfit.f90 \
$(RETR_SDIR)/prior.f90 \
$(RETR_SDIR)/mo_retrieval.f90

RETR_PSRCS = \
$(RETR_SDIR)/retrieval_control.f90 \
$(RETR_SDIR)/jacobian.f90 \
$(RETR_SDIR)/ua.f90

ACTIVE_SRCS += $(MODEL_ASRCS) $(SEMID_ASRCS) $(SENSE_ASRCS)

PASSIVE_SRCS = \
$(MDIR)/init-finish.f90 \
$(MDIR)/sensimul_io.f90



#-----------------------------
#
#     object file naming scheme
#
ACTIVE_OBJS = $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f,%.o,$(filter %.f,$(ACTIVE_SRCS)))))
ACTIVE_OBJS += $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$(filter %.f90,$(ACTIVE_SRCS)))))

PASSIVE_OBJS = $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$(filter %.f90,$(PASSIVE_SRCS)))))
PASSIVE_OBJS += $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.F90,%.o,$(filter %.F90,$(PASSIVE_SRCS)))))

RETR_OBJS = $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$(filter %.f90,$(RETR_ASRCS)))))
RETR_OBJS += $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.f90,%.o,$(filter %.f90,$(RETR_PSRCS)))))


#-----------------------------
#
#     derivative code sources
#
#
FWEXT = _fw
BWEXT = _bw
FWVEXT = $(FWEXT)v
DVFWVSRC = $(DVSRCDIR)/cost$(FWVEXT).f90
DVBWSRC  = $(DVSRCDIR)/cost$(BWEXT).f90
DVFWVOBJ = $(addprefix $(OBJDIR)/, $(notdir $(DVFWVSRC:.f90=.o)) )
DVBWOBJ  = $(addprefix $(OBJDIR)/, $(notdir $(DVBWSRC:.f90=.o)) )
DVOBJS   = $(DVBWOBJ) $(DVFWVOBJ) $(OBJDIR)/ad_support.o $(OBJDIR)/DIFFSIZES.o $(OBJDIR)/adBuffer.o $(OBJDIR)/adStack.o


#=======================================
#
#          compiler settings
#
FC  = gfortran

ifeq "$(FC)" "gfortran"
  FFLAGS = -Wall -ffixed-line-length-132 -fcheck=all -fbacktrace -ffpe-trap=invalid,overflow,underflow,zero,invalid
  FFLAGS = -Wall -ffixed-line-length-132 -fcheck=all -fbacktrace -ffpe-trap=invalid
  FFLAGS = -Wall -ffixed-line-length-132
  FFOPTI = -O3 -funroll-loops -floop-block -ftree-loop-distribution -march=native
  ifdef DBG
    FFLAGS = -Wall -Wextra -Wpedantic
    FFLAGS += -ggdb -Warray-temporaries -Wconversion -fimplicit-none -ffree-line-length-0 -ffixed-line-length-132 -fcheck=all -fbacktrace -finit-real=snan -ffpe-trap=invalid,overflow,zero,invalid
    # FFLAGS += -ggdb -Warray-temporaries -Wconversion -ffree-line-length-0 -ffixed-line-length-132 -fcheck=all -fbacktrace -ffpe-trap=invalid,overflow,zero,invalid
    FFOPTI = -Og
  endif
  ifdef PROFILE
    FFLAGS += -ggdb -pg
    FFOPTI = -O2
  endif
endif

#-------------------
# NOTE::inversion-interface uses
#       double-precision (kind=8)
#
R8 := 1
ifeq "$(FC)" "gfortran"
  ifdef R8
    FFLAGS += -fdefault-real-8 -fdefault-double-8
  endif
endif

#-------------------
#     final compile-flag assembly
#
# Note:
# '-I' also defines the search-path for modules
# (at least with gfortran)
FFINC  = -I/usr/include -I$(OBJDIR)
# Note:
# - all generated '.mod' files will reside in OBJDIR
# - potentially expecting above 2GB global memory
#
FFLAGS += $(FFINC) -J$(OBJDIR) $(FFOPTI)


#=======================================
#
#          libraries to link
#
LDNCLIB = -lnetcdff -lnetcdf




##################################################
#
#     l o a d   l o c a l   s e t t i n g s
#
-include Make-local.inc
-include Make-debug.inc
##################################################


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#                    T A R G E T S
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
help:
	@echo
	@printf "%s\n" "List of targets:"
	@echo
	@printf " %s\n" "Main targets:"
	@echo
	@printf "  %-25s: %s\n" "make setup" \
                "finalise retrieval tool installation"
	@printf "  %-25s  %s\n" "" \
                "Finalising The installation of the retrieval tool requires that the"
	@printf "  %-25s  %s\n" "" \
	        "signaturesimulator is properly installed on the system."
	@echo
	@printf "  %-25s: %s\n" "make pre" \
                "run preprocessing step"
	@echo
	@printf "  %-25s: %s\n" "make retrieval" \
                "run core retrieval and target steps"
	@echo
	@printf "  %-25s: %s\n" "make mba" \
                "run analysis of uncertainty only"
	@echo
	@printf "  %-25s: %s\n" "make synthetic" \
                "generate pseudo observations"
	@echo
	@printf " %s\n" "Main targets:"
	@echo
	@printf "  %-25s: %s\n" "make clean" \
                "remove many reconstructable files (as e.g. object files,...)"
	@echo
	@printf "  %-25s: %s\n" "make scratch" \
                "should remove all reconstructable files"

x_help:
	@echo
	@printf " ***Running pre-configured experiments***\n"
	@printf " ----------------------------------------\n"
	@echo
	@printf "  %-25s: %s\n" "make complete" \
                "Retrieval system will be run with observations from both S1 and S2"
	@printf "  %-25s  %s\n" "" \
                "and with dynamical model."
	@printf "  %-25s  %s\n" "" \
                "Result files will reside in reresdir_complete_$(EXPV)"
	@printf "  %-25s  %s\n" "" \
                "upon successful completion."
	@echo
	@printf "  %-25s: %s\n" "make nomodel" \
                "Retrieval system will be run as in case complete but with the"
	@printf "  %-25s  %s\n" "" \
                "dynamical model excluded."
	@printf "  %-25s  %s\n" "" \
                "Result files will reside in reresdir_nomodel_$(EXPV)"
	@printf "  %-25s  %s\n" "" \
                "upon successful completion."
	@echo
	@printf "  %-25s: %s\n" "make nos1" \
                "Retrieval system will be run as in case complete but with"
	@printf "  %-25s  %s\n" "" \
                "S1 observations excluded."
	@printf "  %-25s  %s\n" "" \
                "Result files will reside in reresdir_nos1_$(EXPV)"
	@printf "  %-25s  %s\n" "" \
                "upon successful completion."
	@echo
	@printf "  %-25s: %s\n" "make nos2" \
                "Retrieval system will be run as in case complete but with"
	@printf "  %-25s  %s\n" "" \
                "S2 observations excluded."
	@printf "  %-25s  %s\n" "" \
                "Result files will reside in reresdir_nos2_$(EXPV)"
	@printf "  %-25s  %s\n" "" \
                "upon successful completion."
	@echo
	@printf "  %-25s: %s\n" "make reference" \
                "Retrieval system will be run as in case complete but without"
	@printf "  %-25s  %s\n" "" \
                "information provided by S1 and S2 and with the dynamic model excluded."
	@printf "  %-25s  %s\n" "" \
                "(this means that only the prior information is used)"
	@printf "  %-25s  %s\n" "" \
                "Result files will reside in reresdir_reference_$(EXPV)"
	@printf "  %-25s  %s\n" "" \
                "upon successful completion."
	@echo
	@echo "=================================================="



#===========================================================
#
#               c r e a t i o n  /  r u n n i n g   o f   d r i v e r s
#
#--
EXPCLEAN_FILES += $(RETRCONFIG_FILE) $(RETRPRIOR_FILE) $(RETRMODEL_FILE) $(RETRCTL_FILE) $(SITE_FILE) $(BOUNDS_FILE)
EXPCLEAN_FILES += x.b sx.b
EXPCLEAN_FILES += fort.[0-9] fort.[0-9][0-9]

PRGS += simulate.x
simulate.x: $(OBJDIR)/runsim.o $(LIBMODEL) $(LIBSUPPORT)
	$(FC) -o $@ $(FFOPTI) $^ $(LDNCLIB)
EXPCLEAN_FILES += sensimul_s?.nc sensimul_s?.b sensimul.b

PRGS += runinv.x
runinv.x: $(OBJDIR)/runlbfgsb.o $(LIBRETR) $(LIBDERIV) $(LIBMODEL) $(LIBSUPPORT) $(LIBLBFGSB)
	$(FC) -o $@ $(FFOPTI) $^ $(LDNCLIB)
EXPCLEAN_FILES += x.b sx.b iterate.dat

#-------------------
#     posterior uncertainty analysis
#
PRGS += runua.x
runua.x: $(OBJDIR)/runua.o $(LIBRETR) $(LIBDERIV) $(LIBMODEL) $(LIBSUPPORT) $(LIBEIS)
	$(FC) $(FFOPTI) -o $@ $^ $(LDNCLIB)
EXPCLEAN_FILES += xrunua.log sigma.dat

PRGS += target_processor.x
target_processor.x: $(OBJDIR)/target_processor.o $(LIBMODEL) $(LIBSUPPORT)
	$(FC) -o $@ $(FFOPTI) $^ $(LDNCLIB)
EXPCLEAN_FILES += sensimul_s?.nc sensimul_s?.b sensimul.b

#===========================================================
#
#               S1 + S2   T A R G E T S
#
#--
setup_old:
	@if [ "x$(SENTINEL_SIMULATOR)" = "x" ]; \
         then \
	   echo; \
           echo "$@::INSTALLATION_ERROR: ***SENTINEL_SIMULATOR*** not defined!"; \
           echo; echo; \
           exit 1; \
         fi; \
         srf_file=$(SENTINEL_SIMULATOR)/srfData/s2a.srf; \
         tle_file=$(SENTINEL_SIMULATOR)/TLE/norad_resource_tle.txt; \
         jules_file=$(SENTINEL_SIMULATOR)/sentinel_simulator/jules/output/sensitivity_runs/crp_g_77_6.8535_0.17727_0.000573.3_hourly.nc; \
         srf_input=$(INPUT_DIR)/s2a.srf; \
         tle_input=$(INPUT_DIR)/tle_file.txt; \
         jules_output=$(INPUT_DIR)/jules_output.nc; \
         mkdir -p $(INPUT_DIR); \
         ln -fs `readlink -f $${srf_file}` $${srf_input}; \
         ln -fs `readlink -f $${tle_file}` $${tle_input}; \
         ln -fs `readlink -f $${jules_file}` $${jules_output}; \
         rm -f $(INPUT_DIR)/.sentinel_simulator; echo $(SENTINEL_SIMULATOR)/sentinel_simulator > $(INPUT_DIR)/.sentinel_simulator

setup:
	$(RS_SETUP_TOOL)

pre:
	( args="";   \
          if [ "x$(PRE_MODE)" = "xpre_synthetic" ]; then \
            args="--time_start $(TIMSTART) --time_end $(TIMEND)"; \
            if [ "x$(NOXTGT)" = "x" ]; \
            then \
              args="$${args} --target_schedule $(TGTSCHEDULE_FILE)"; \
            fi; \
	    if [ "x$(MISSION_LST)" != "x" ]; \
            then \
	      args="$${args} --mission_lst $(MISSION_LST)"; \
            fi; \
          fi; \
          args="$${args} $(PREARGXTRA)"; \
	  $(RS_PRE_TOOL) $(PRE_MODE) $${args} |& tee mk.$@.log; \
        )
.PHONY: pre
EXPCLEAN_FILES += mk.pre.log

retrieval: runinv.x $(RETRCONFIG_FILE) $(RETRPRIOR_FILE) $(RETRMODEL_FILE) $(RETRCTL_FILE) $(SITE_FILE)
	./runinv.x |& tee mk.retrieval.log
.PHONY: retrieval
EXPCLEAN_FILES += iterate.dat mk.retrieval.log

mba: runua.x $(RETRCONFIG_FILE) $(RETRPRIOR_FILE) $(RETRMODEL_FILE) $(RETRCTL_FILE) $(SITE_FILE)
	./runua.x |& tee mk.mba.log
.PHONY: mba
EXPCLEAN_FILES += sigma.dat controlvector_prior.nc controlvector_post.nc mk.mba.log

target: target_processor.x $(RETRCONFIG_FILE) $(RETRPRIOR_FILE) $(RETRMODEL_FILE) $(RETRCTL_FILE) $(SITE_FILE)
	./$< --target fapar
.PHONY: target

synthetic: simulate.x $(RETRCONFIG_FILE) $(RETRPRIOR_FILE) $(SITE_FILE)
	./$< $(SYNXTRA) |& tee mk.synthetic.log
	if [ -f sensimul_s1.nc ]; then mv sensimul_s1.nc obs_s1.nc; fi
	if [ -f sensimul_s2.nc ]; then mv sensimul_s2.nc obs_s2.nc; fi
.PHONY: synthetic
EXPCLEAN_FILES += obs_s1.nc obs_s2.nc simulation.nc mk.synthetic.log


#===========================================================
#
#               S1 + S2   E X P E R I M E N T S
#
#--
reference: runinv.x runua.x
	$(MAKE) clean_exp
	$(MAKE) pre PRE_MODE=pre_synthetic PREARGXTRA="--no_use_model"
	$(MAKE) retrieval
	$(MAKE) mba
	$(MAKE) save_exp SAVDIR=resdir_$@_$(EXPV)
.PHONY: reference

complete: simulate.x runinv.x runua.x
	$(MAKE) clean_exp
	$(MAKE) pre PRE_MODE=pre_synthetic	
	$(MAKE) synthetic
	$(MAKE) retrieval
	$(MAKE) mba
	$(MAKE) save_exp SAVDIR=resdir_$@_$(EXPV)
.PHONY: complete

nomodel: simulate.x runinv.x runua.x
	$(MAKE) clean_exp
	$(MAKE) pre PRE_MODE=pre_synthetic PREARGXTRA="--no_use_model"
	$(MAKE) synthetic
	$(MAKE) retrieval
	$(MAKE) mba
	$(MAKE) save_exp SAVDIR=resdir_$@_$(EXPV)
.PHONY: nomodel

nos1: simulate.x runinv.x runua.x
	$(MAKE) clean_exp
	$(MAKE) pre PRE_MODE=pre_synthetic MISSION_LST="S2A S2B"
	$(MAKE) synthetic
	$(MAKE) retrieval
	$(MAKE) mba
	$(MAKE) save_exp SAVDIR=resdir_$@_$(EXPV)
.PHONY: nos1

nos2: simulate.x runinv.x runua.x
	$(MAKE) clean_exp
	$(MAKE) pre PRE_MODE=pre_synthetic MISSION_LST="S1A S1B"
	$(MAKE) synthetic
	$(MAKE) retrieval
	$(MAKE) mba
	$(MAKE) save_exp SAVDIR=resdir_$@_$(EXPV)
.PHONY: nos2

experiments_all:
	@( \
           mkdir -p $@-$(EXPV); \
           targets="reference complete nomodel nos1 nos2"; \
           for t in `echo $${targets}`; do \
             $(MAKE) $$t; \
             mv resdir_$${t}_$(EXPV) $@-$(EXPV); \
           done; \
           tar zcf $@-$(EXPV).tgz $@-$(EXPV); \
         )

#===========================================================
#
#               a u x i l i a r y   t a r g e t s
#
#--
# $(RETRCONFIG_FILE) $(RETRPRIOR_FILE) $(RETRMODEL_FILE) $(RETRCTL_FILE) $(SITE_FILE) $(BOUNDS_FILE): pre
$(RETRCTL_FILE):
	cp $(INPUT_DIR)/retrctl.nml $@

save_exp:
	if [ "x$(SAVDIR)" != "x" ]; then \
          mkdir -p $(SAVDIR); \
          $(RM) $(SAVDIR)/*; \
          for f in `echo $(EXPCLEAN_FILES)`; \
          do \
	    if [ -f $${f} ]; \
            then \
              cp -a $${f} $(SAVDIR); \
            fi; \
          done; \
        fi;

#===========================================================
#
#               l i b r a r i e s
#
#-------------------
#-----------------------------
#
#     SEMIDISCRETE model libraray
#
$(LIBSEMID): $(SEMID_OBJS)
	@mkdir -p $(LIBDIR)
	ar -r $@ $?
libsemid:
	$(MAKE) -B $(LIBSEMID)
.PHONY: libsemid

$(LIBSENSE): $(SENSE_OBJS)
	@mkdir -p $(LIBDIR)
	ar -r $@ $?
libsense:
	$(MAKE) -B $(LIBSENSE)
.PHONY: libsense

$(LIBMODEL): $(ACTIVE_OBJS) $(PASSIVE_OBJS)
	@mkdir -p $(LIBDIR)
	ar -r $@ $?
libmodel:
	$(MAKE) -B $(LIBMODEL)
.PHONY: libmodel

$(LIBRETR): $(RETR_OBJS)
	@mkdir -p $(LIBDIR)
	ar -r $@ $?
libretr:
	$(MAKE) -B $(LIBRETR)
.PHONY: libretr

$(LIBDERIV): $(DVOBJS)
	@mkdir -p $(LIBDIR)
	ar -r $@ $?
libderiv:
	$(MAKE) -B $(LIBDERIV)
.PHONY: libderiv

$(LIBSUPPORT): $(SUPPORT_OBJS)
	@mkdir -p $(LIBDIR)
	ar -r $@ $?
libsupport:
	$(MAKE) -B $(LIBSUPPORT)
.PHONY: libsupport

$(LIBLBFGSB): FORCE
	@mkdir -p $(LIBDIR)
	$(MAKE) -C $(LIBLBFGSBDIR) ../$@ FC="$(FC)" FFLAGS="$(FFLAGS)" LIBNAME="../$@"

liblbfgsb:
	$(MAKE) -B $(LIBLBFGSB)
.PHONY: liblbfgsb

$(LIBEIS): FORCE
	@mkdir -p $(LIBDIR)
	$(MAKE) -C $(LIBEISDIR) ../$@ FC=$(FC) FFLAGS="$(FFLAGS)" LIBNAME="../$@"
libeis:
	$(MAKE) -B $(LIBEIS)
.PHONY: libeis


#===========================================================
#
#               p a t t e r n   r u l e s
#
#-------------------
$(OBJDIR)/%.o: $(SEMID_SDIR)/%.f
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(OBJDIR)/%.o: $(SEMID_SDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(OBJDIR)/%.o: $(SENSE_SDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(OBJDIR)/%.o: $(MDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(OBJDIR)/%.o: $(RETR_SDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(OBJDIR)/%.o: $(DRVDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(OBJDIR)/%.o: $(SUPPORT_SDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(MDIR)/%.f90: $(MDIR)/%.F90
	$(FC) -cpp -E -P $(CPPFLAGS) -o $@ $<

$(OBJDIR)/%.o: $(DVSRCDIR)/%.f
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(OBJDIR)/%.o: $(DVSRCDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

$(OBJDIR)/adStack.o: $(DVSRCDIR)/adStack.c
	@mkdir -p $(OBJDIR)
	$(CC) -fPIC -I$(DVSRCDIR) $(CFLAGS) $(CFXTRA) -c -o $@ $<


#===========================================================
#
#               e x t r a   d e p e n d e n c i e s
#
#-------------------
$(OBJDIR)/canopy.o $(OBJDIR)/core.o $(OBJDIR)/ground.o $(OBJDIR)/Oh92.o $(OBJDIR)/scatterer.o $(OBJDIR)/sense_tests.o: $(OBJDIR)/mo_sensimul_s1.o
$(OBJDIR)/mo_sensimul.o: $(OBJDIR)/mo_sensimul_s1.o $(OBJDIR)/mo_sensimul_s2.o
$(OBJDIR)/init-finish.o: $(OBJDIR)/mo_sensimul.o $(OBJDIR)/mo_sensimul_s1.o $(OBJDIR)/mo_sensimul_s2.o
$(OBJDIR)/sensimul_io.o: $(OBJDIR)/mo_sensimul.o $(OBJDIR)/mo_timing.o $(OBJDIR)/init-finish.o
$(OBJDIR)/evalf.o: $(OBJDIR)/mo_sensimul.o $(OBJDIR)/mo_sensimul_s1.o $(OBJDIR)/mo_sensimul_s2.o
$(OBJDIR)/model.o: $(OBJDIR)/mo_sensimul.o
$(DVBWOBJ): $(OBJDIR)/mo_sensimul.o $(OBJDIR)/mo_sensimul_s1.o $(OBJDIR)/mo_sensimul_s2.o $(OBJDIR)/model.o $(OBJDIR)/price_soil_ftn.o
$(DVFWVOBJ): $(OBJDIR)/DIFFSIZES.o
$(OBJDIR)/retrieval_control.o: $(OBJDIR)/mo_retrieval.o $(OBJDIR)/mo_sensimul.o
$(OBJDIR)/cost.o $(DVBWOBJ) $(DVFWVOBJ): $(OBJDIR)/mo_retrieval.o
$(OBJDIR)/misfit.o: $(OBJDIR)/mo_sensimul.o
$(OBJDIR)/jacobian.o: $(OBJDIR)/DIFFSIZES.o $(OBJDIR)/retrieval_control.o
$(OBJDIR)/runsim.o: $(OBJDIR)/mo_sensimul.o
$(OBJDIR)/target_processor.o: $(OBJDIR)/mo_sensimul.o


#===========================================================
#
#               f o r c e   t a r g e t
#
#-------------------
FORCE:


#===========================================================
#
#               c l e a n i n g
#
#--
clean_exp:
	$(RM) $(EXPCLEAN_FILES)

clean:
	find -name ".*~" -delete
	find -name "*~" -delete
	find -name "fort.*" -delete
	$(RM) $(OBJDIR)/* $(LIBDIR)/*
	$(RM) $(CLEAN_FILES)
	$(RM) -r $(CLEAN_DIRS)
	$(MAKE) -C $(LIBLBFGSBDIR) clean
	$(MAKE) -C $(LIBEISDIR) clean
.PHONY: clean


scratch: clean clean_exp
	$(RM) $(PRGS)
	$(RM) -r $(OBJDIR) $(LIBDIR)
	$(RM) $(SCRATCH_FILES)
	$(RM) -r $(SCRATCH_DIRS)
.PHONY: scratch
