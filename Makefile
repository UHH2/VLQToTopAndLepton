LIBRARY := SUHH2VLQToTopAndLepton
DICT := include/BprimeContainer.h include/BprimeGenContainer.h include/SUHH2VLQToTopAndLepton_LinkDef.h
USERLDFLAGS := -lSUHH2core -lSUHH2common -lGenVector -lMinuit
# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
include ../Makefile.common
