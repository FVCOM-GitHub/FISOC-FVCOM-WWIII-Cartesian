# FISOC-FVCOM-WWIII-XY
FISOC_FVCOM_WWIII_XY is a coupled modeling system that integrates FISOC, FVCOM, and WaveWatch III (WW3) using Cartesian coordinates.


FISOC_FVCOM_WWIII_XY

Cartesian-Coordinate FISOC–FVCOM–WW3 Coupler Model
Modified from FISOC_FVCOM_SWAN by Jianhua Qi, SMAST, UMassD (2025)

Overview

FISOC_FVCOM_WWIII_XY is a coupled modeling system that integrates FISOC, FVCOM, and WaveWatch III (WW3) using Cartesian coordinates.

The package includes:

* FISOC_FVCOM_WWIII — Main FISOC ESMF coupling framework

* FVCOM51_src — FVCOM source code

* WW3 — WaveWatch III model package (unstructured grid; Cartesian coordinates; modifications for Cartesian coordinates are marked with !JQI in the WW3 source code)

1. Compilation Instructions

FVCOM and WW3 must be compiled separately.

Compile FVCOM

(1) Edit make.inc as needed.

(2) Compile using:

    make

Compile WW3

(1) Read the instructions in:

    WW3/note_ww3_li.txt

(2) Follow the steps provided there to compile WW3.

Compile FISOC-FVCOM-WW3

    In the folder FISOC_FVCOM_WWIII, run:

    bash ./buildFISOC_FVCOM_WW3.sh

2. WW3 Preprocessing Workflow

Run the following utilities in the folder ww3_input:

    Mesh generation
    Run ww3_grid in:

    run_mod_def/

    Wind forcing
    Run ww3_prnc in:

    run_prnc/

    Nesting boundary generation
    Run ww3_bounc in:

    run_bound/

3. Running the FISOC–FVCOM–WW3 Coupled Model

(1) In the run directory (for example Ex6_Inlets_realtime), edit the configuration file:

    FISOC_config.rc

(2) Link the executable FISOC_caller_FVCOM_WW3 in the run directory.

(3) Run the coupled model driver:

    mpiexec -np 4 ./FISOC_caller_FVCOM_WW3

4. WW3 Postprocessing (for large spherical coordinates WW3 domain used for small-domain nesting; not included in this package)

   In WW3, within the folder run_shel, two postprocessing utilities are provided after model run:

   * field_netcdf/ -- Creates whole-domain NetCDF output.

   * point_netcdf/ -- Creates point-based NetCDF output for small-domain nesting.

   Steps for nesting preparation:

   (1) Copy the large-domain file, for example:

    ww3.200001_spec.nc

    to:

    FISOC_FVCOM_WW3/Ex6_Inlets_realtime(small domain)/ww3_input/run_bound/

    (2) Add x and y coordinates to the file ww3.200001_spec.nc.

    (3) Run ww3_bounc to generate nest.ww3 for the small domain model.

    Note: The Ex6_Inlets_realtime example does not include a nesting boundary.

5. FISOC–FVCOM–WW3 Postprocessing

   To generate small-domain WW3 output for the FISOC–FVCOM–WW3 run, use the postprocessing utility in:

   ww3_input/run_shel/field_netcdf/

   or find the output file of FVCOM in the folder FISOC_FVCOM_WWIII_XY/FISOC_FVCOM_WWIII/Ex6_Inlets_realtime/output/

