# FVCOM–WAVEWATCH III water-level and current coupling fixes

## Scope

This note documents the changes from `FISOC_FVCOM_WWIII_org` to
`FISOC_FVCOM_WWIII_wslcur_261117`, based on a direct comparison of the
`FISOC-devel/src` source files. This is the version preceding the sea-ice work.
AICE exchange, IC1 attenuation, Cartesian ice initialization, and ice-dynamics
changes are outside the scope of this update.

## Problem observed

During coupling tests, water levels and currents supplied by FVCOM were not
being reflected as expected in the WW3 calculation. Matching ESMF export/import
fields alone did not establish that the imported values reached WW3's native
forcing arrays and were used at the intended time.

The changes therefore address both the mapping from distributed coupling fields
to native WW3 arrays and the forcing-time handling in the WW3 driver.

## Modified source files

| File | Main changes | Reason |
| --- | --- | --- |
| `FISOC_WM_Wrapper_WW3.f90` | Collect owned-node inputs by global ID; copy into native water-level/current arrays; implement bathymetry assignment; add checks; adjust selected wave-output conversions. | Correct array layout, pointer lifetime, and forcing assignment; expose inconsistent inputs. |
| `FISOC_WW3.f90` | Set coupled forcing times and bypass ordinary input-file updates for coupled water levels/currents. | Make native WW3 updates use the fields provided for the current coupling interval. |
| `FISOC_utils.f90` | Make the destination ESMF redistribution array caller-owned and reusable. | Avoid repeated array creation and provide explicit lifetime management. |
| `FISOC_OM_Wrapper_FVCOM.f90` | Use and finalize the reusable ESMF redistribution array. | Keep the FVCOM receiver consistent with the revised utility interface. |
| `FISOC_caller.f90` | Use `ESMF_LOGKIND_MULTI_ON_ERROR`. | Reduce routine per-PET logging while retaining error diagnostics. |

All file paths below are relative to `FISOC-devel/src`.

## 1. Native WW3 input assignment

**File:** `FISOC_WM_Wrapper_WW3.f90`  
**Routine:** `sendFieldDataToWM`

### Previous implementation

The receiver allocated temporary arrays with shape `(NPA,2)`, filled them with
local imported values, reassociated WW3 pointers such as `WLEV`, `CXN`, and `CYN`
with these arrays, and then deallocated the temporary arrays.

This had two important problems:

- `NPA` describes a local partition, including halo nodes. It is not the layout
  of the native global forcing arrays used by this unstructured WW3 setup.
- Reassociating a pointer with a temporary array and then deallocating that
  array leaves the pointer without a valid target.

### Revised implementation

1. Redistribute the imported ESMF field to the WW3 local layout.
2. Check that the imported pointer has `NPA` entries.
3. Use owned nodes `1:NP`, with `ivertp(j)=iplg(j)`, to place local values at their
   global node IDs in a zero-initialized buffer.
4. Collect the buffer with `MPI_Allreduce` on the WW3 component communicator.
   Halo nodes do not contribute a second time.
5. Copy the resulting global field into WW3-owned arrays without replacing their
   pointer targets.

| FVCOM field | WW3 destination | Treatment |
| --- | --- | --- |
| `OM_ETA` | `WLEV(:,1)` | Copy global water levels. |
| `OM_UU` | `CX0(:,1)` and `CXN(:,1)` | Set both current endpoints to the received x component. |
| `OM_VV` | `CY0(:,1)` and `CYN(:,1)` | Set both current endpoints to the received y component. |
| `OM_HH` | `ZB(ii)` | Use `MAPSF(ii,1)` to find the global node and assign `-depth`. |

The assignment targets the unstructured configuration with forcing-array shape
`(NX,1)`. Explicit shape and global-ID checks were added. This change should not
be assumed to support arbitrary structured-grid layouts without adaptation.

The current endpoints have identical values, so currents are held constant
within each coupling interval. No interpolation between successive FVCOM
coupling snapshots was introduced.

Water-level and current imports are gated by `FLAGSC(1)` and `FLAGSC(2)`.
Configured coupled inputs are checked for receipt before returning. `OM_HH`
requires coupled water levels because this implementation uses the water-level
update path to refresh WW3 depth-related quantities.

## 2. Coupled forcing times and input-file handling

**File:** `FISOC_WW3.f90`  
**Routine:** `FISOC_WW3_RUN`

The following assignments were added before advancing WW3:

```fortran
IF (FLAGSC(1)) TLN = TIME0
IF (FLAGSC(2)) THEN
  TC0 = TIME0
  TCN = TIMEN
END IF
```

They associate the imported fields with the integration window. The driver also
skips the ordinary forcing-file update loop for water levels and currents when
those fields are configured as coupled inputs:

```fortran
IF (J >= 1 .AND. J <= 2) THEN
  IF (FLAGSC(J)) CYCLE
END IF
```

This prevents the normal file-input path from being used for fields already
supplied by the coupler.

The corresponding run configuration must enable:

```fortran
INPUT%FORCING%WATER_LEVELS = 'C'
INPUT%FORCING%CURRENTS     = 'C'
```

The FISOC requested/exported variable lists must also include the desired fields.
Source changes alone do not enable a field omitted from the run configuration.

## 3. Reusable ESMF redistribution arrays

**Files:** `FISOC_utils.f90`, `FISOC_OM_Wrapper_FVCOM.f90`,
`FISOC_WM_Wrapper_WW3.f90`

`FISOC_ArrayRedistFromField` previously created a destination ESMF array on every
call. Its interface now accepts `destArray` and `destArrayReady` from the caller.
The array is created on first use and reused for subsequent redistribution calls.

Both model wrappers retain the array as module state and destroy it during
finalization. This removes repeated creation of destination ESMF arrays and
makes their ownership explicit. It is a memory/lifetime improvement, separate
from the local-to-global WW3 mapping fix.

Because the utility interface changed, the utility and both wrappers must be
rebuilt together. The reuse assumes the destination distribution remains fixed.

## 4. Additional wave-variable conversions and remaining limitations

**File:** `FISOC_WM_Wrapper_WW3.f90`  
**Routine:** `getFieldDataFromWM`

- `WM_DIR`: convert `THM` from radians to degrees and normalize to `[0,360)` for
  the FVCOM receiver's direction convention. This is an angular-unit conversion;
  it does not add a geographic rotation or a from/to direction reversal.
- `WM_TPEAK`: replace direct transfer of peak frequency `FP0` with
  `1.0/(FP0+1.0E-8)` to represent a period in seconds.

The peak-period expression is an interim implementation: zero or invalid peak
frequency, including dry/inactive nodes, still needs explicit treatment. The
small denominator offset is not a physical solution for those cases.

The existing placeholder mappings for `WM_TMBOT`, `WM_UBOT`, `WM_DIRBOT`, and
`WM_QB1` remain unresolved and are marked for further work. They should not be
presented as validated wave inputs for FVCOM physics that requires those fields.

## 5. Logging

**File:** `FISOC_caller.f90`

Change the ESMF log kind from `ESMF_LOGKIND_MULTI` to
`ESMF_LOGKIND_MULTI_ON_ERROR`. This is a diagnostic-output change and does not
alter the coupling physics. Error-related per-PET logs may still be generated.
Some node-specific debug writes in the replaced WW3 import code were also
commented out.

## Validation status

The author reports successful water-level/current coupling tests after these
changes. This README does not claim an independent rerun or a comprehensive
regression/physics validation. Quantitative results, compiler/library versions,
MPI process counts, and the exact test configuration should accompany any formal
validation report.

For review or reproduction, distinguish these checks:

1. Compare FVCOM exports and WW3 imports at corresponding coordinates and times.
2. Confirm native WW3 water-level and current output (`WLV`, `CUR`) reflects the
   received inputs.
3. Compare wave solutions with matched forcing and ocean feedback enabled versus
   disabled.

Agreement at step 1 alone does not demonstrate steps 2 or 3.

`FISOC_coupler.f90` and `FISOC_OM.f90` are unchanged between the two source
snapshots compared here. Separate diagnostic-output or FVCOM-output changes
discussed during debugging are not included in this source-difference summary.
