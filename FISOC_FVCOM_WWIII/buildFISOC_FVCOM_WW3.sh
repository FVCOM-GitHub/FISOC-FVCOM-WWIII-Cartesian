#!/usr/bin/env bash
                                                                                                            
export CPPFLAGS="$CPPFLAGS -D FISOC_MPI"

export FISOC_EXE="FISOC_caller_FVCOM_WW3"

#export FFLAGS="$FFLAGS -fbacktrace -g -O0 -fbounds-check -Wall"
export FFLAGS=" -O0 -g -traceback"
#export FFLAGS=" -g -check all -fpe0 -warn -traceback -debug extended"
#export FFLAGS=" -O3 -xHost " #-ipo"
export ELMER_HOME="/opt/Elmer_test"
## optionally over-write the default executable name:
#export FISOC_EXE="FISOC_caller"
#export PROJLIBS="-L/opt/proj-oneapi-2024/lib -lfproj4 -lproj -lm"
#export PROJINCS="-I/opt/proj-oneapi-2024/include"

export FISOC_ISM="dummy"
export FISOC_ISM_LIBS=""
export FISOC_ISM_INCLUDE="-I/usr/local/include"
export FISOC_ISM_LIBPATH=""
export FISOC_ISM_GEOM=""


#export FISOC_ISM="Elmer"
#export FISOC_ISM_LIBS="-lelmersolver -lmatc -lfhuti -larpack -lparpack"
#export FISOC_ISM_INCLUDE="$ELMER_HOME/share/elmersolver/include"
#export FISOC_ISM_LIBPATH="$ELMER_HOME/lib/elmersolver/"
#export FISOC_ISM_GEOM="FISOC_ISM_MESH"

export FISOC_AM="dummy"
export FISOC_AM_LIBS=""
export FISOC_AM_LIBPATH=""
export FISOC_AM_INCLUDE=""

#-L/home/jzhge/libs_install/lib -lmetis

#export FISOC_INCPATHS="-I/opt/model_libs/include -I/opt/netcdf-oneapi-2024/include"
export FISOC_INCPATHS="-I/usr/local/share/x86_64/libjulian/1.3.3/intel/2023.2.0/include -I/usr/local/share/x86_64/metis/5.1.0/intel/2023.2.0/include -I/usr/local/share/x86_64/parmetis/4.0.3/intel/2023.2.0/include"
#export FISOC_LIBPATHS="-L/opt/model_libs/lib"
export FISOC_LIBPATHS="-L/usr/local/share/x86_64/libjulian/1.3.3/intel/2023.2.0/lib64 -L/usr/local/share/x86_64/metis/5.1.0/intel/2023.2.0/lib -L/usr/local/share/x86_64/parmetis/4.0.3/intel/2023.2.0/lib"
export FISOC_LIBS="-ljulian -lmetis -lparmetis"


export FISOC_OM="FVCOM"
export FISOC_OM_LIBS="-lfvcom"
export FISOC_OM_INCLUDE="../FVCOM51_src"
export FISOC_OM_LIBPATH="../FVCOM51_src"
export FISOC_OM_GEOM="FISOC_OM_MESH"


export FISOC_WM="WW3"
export FISOC_WM_LIBS="-lww3"
export FISOC_WM_INCLUDE="../WW3/built/install/mod"
export FISOC_WM_LIBPATH="../WW3/built/install/lib"
export FISOC_WM_GEOM="FISOC_WM_MESH"
#export FISOC_WM_PARMETIS="/usr/local/share/x86_64/parmetis/4.0.3/intel/2023.2.0/lib"
#export FISOC_WM_METIS_LIBS="-lparmetis"

#export FISOC_WM="SWAN"
#export FISOC_WM_LIBS="-lswan"
#export FISOC_WM_INCLUDE="../swan4145_pun"
#export FISOC_WM_LIBPATH="../swan4145_pun"
#export FISOC_WM_GEOM="FISOC_WM_MESH"

#export  ESMFMKFILE=/opt/esmf-8.6.0-mpich4.2/lib/libO/Linux.intel.64.mpich.default/esmf.mk
export ESMFMKFILE="/usr/local/share/x86_64/esmf/8.5.0/intel/2023.2.0/lib/libO/Linux.intel.64.mvapich2.default/esmf.mk"
#export ESMFMKFILE="/opt/esmf-8.0.1/lib/libO/Linux.intel.64.mpich3.default/esmf.mk"
#export ESMFMKFILE="/opt/esmf_7.1.0r/lib/libO/Linux.intel.64.mpich3.default/esmf.mk"

make clean
make install

#rm PET*Log

# to run, for example:
# mpirun -np 4 FISOC_caller

