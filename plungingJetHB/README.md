# plungingJetHB

Production-oriented OpenFOAM Foundation v14 case for an **axisymmetric circular plunging jet of a Herschel–Bulkley liquid** entering a quiescent liquid pool through air.

## Verified OpenFOAM Foundation v14 references

This case is assembled from verified OpenFOAM Foundation v14 patterns:

- **Solver module:** `incompressibleVoF`, executed through `foamRun` with `solver incompressibleVoF;` in `system/controlDict`.
- **Closest geometry/boundary reference:** `tutorials/incompressibleVoF/nozzleFlow2D` for the wedge mesh, inlet `alpha`/`U`, and open boundaries.
- **Canonical incompressibleVoF syntax:** `tutorials/incompressibleVoF/damBreakLaminar` for `phaseProperties`, `physicalProperties.<phase>`, `momentumTransport`, `fvSchemes`, and `fvSolution` patterns.
- **Two-phase transport split:** `tutorials/incompressibleVoF/climbingRod` for the `twoPhaseTransport` / `momentumTransport.<phase>` layout.
- **Monitoring examples:** `tutorials/incompressibleVoF/waterChannel` and `tutorials/incompressibleVoF/wave` for runtime function-object patterns.

OpenFOAM Foundation v14 no longer uses a single solver-side `constant/transportProperties` file for this class of case. The actual runtime inputs are:

- `constant/phaseProperties`
- `constant/physicalProperties.liquid`
- `constant/physicalProperties.air`
- `constant/momentumTransport`
- `constant/momentumTransport.liquid`
- `constant/momentumTransport.air`

The included `constant/transportProperties` and `constant/turbulenceProperties` files are compatibility notes for users expecting the older layout.

## Directory contents

```text
plungingJetHB/
├── 0/
├── constant/
├── system/
├── scripts/renderCase.py
├── Allclean
├── Allmesh
├── Allrun
├── Allrun.parallel
├── Allpost
└── README.md
```

Custom rheology support is provided alongside the case in the repository root:

```text
applications/regularisedHerschelBulkley/
build.sh
```

## User-editable inputs

Edit **only** `constant/caseParameters`.

### Parameter conventions

- `jetDiameter` is the physical jet diameter **D** in metres.
- `jetStandOffDistance`, `poolDepth`, `gasHeight`, and `domainRadius` are specified in multiples of **D**.
- `wedgeAngle` is in degrees.
- Herschel–Bulkley inputs are dynamic-property inputs in SI units.

## Geometry and initial condition

The mesh is a **narrow 3D wedge** about the symmetry axis, not a planar 2D mesh.

- axial direction: `z`
- radial direction: `r`
- axis: `r = 0`
- inlet: circular top patch for `r <= D/2`
- atmosphere: top annulus outside the inlet
- outer boundary: open radial patch
- bottom: no-slip wall

The initial pool occupies `0 <= z <= poolDepth`.

The top boundary is the nozzle supply boundary. To let `jetStandOffDistance` remain independently editable while keeping the inlet on the top patch, the case optionally pre-fills a short liquid column below the inlet. Its lower tip starts `jetStandOffDistance` above the initial pool surface. Therefore `gasHeight` must be greater than or equal to `jetStandOffDistance`.

## Herschel–Bulkley implementation

OpenFOAM Foundation v14 ships native generalized-Newtonian support, but its stock `HerschelBulkley` model is limited to the built-in formulation and does not natively write the required `gammaDot` and dynamic `muEff` fields for this two-phase VOF use-case. This repository therefore adds a minimal custom runtime viscosity model:

- library: `libregularisedHBViscosityModel.so`
- source: `applications/regularisedHerschelBulkley/`
- selected in: `constant/momentumTransport.liquid`

### Implemented constitutive law

For the **liquid phase only**:

```text
gammaDot = sqrt(2 D:D)
gammaReg = max(gammaDot, regularizationParameter)
muRaw    = tauY/gammaReg + K*gammaReg^(n - 1)
muEff    = clip(muRaw, muMin, muMax)
nuEff    = muEff/rhoLiquid
```

Notes:

- `regularizationParameter > 0` removes the zero-shear singularity.
- `muEff` is clipped to `[muMin, muMax]`.
- the **air phase remains Newtonian** with constant `muAir`.
- the VOF momentum equation uses the standard incompressibleVoF two-phase mixture treatment; interface cells blend phase properties through the native phase-fraction formulation, while the liquid rheology is evaluated from the liquid-phase velocity field and supplied as `nuEff` for the liquid stress model.
- runtime output fields are written as:
  - `gammaDot.liquid`
  - `muEff.liquid`
  - native generalized-Newtonian `nu.liquid`

## Mesh generation

`Allmesh` and both run scripts call `scripts/renderCase.py` first. That script reads `constant/caseParameters` and regenerates all parameter-dependent dictionaries, including:

- `system/blockMeshDict`
- `system/setFieldsDict`
- `system/controlDict`
- `constant/phaseProperties`
- `constant/physicalProperties.*`
- `constant/momentumTransport*`

This keeps the case parameterized from one file.

## Boundary and initial conditions

### `alpha.liquid`

- `inlet`: `fixedValue 1`
- `bottom`: `zeroGradient`
- `outer`, `atmosphere`: `inletOutlet` with `0` inflow
- `front`, `back`: `wedge`
- `axis`: `empty`

### `U`

- `inlet`: prescribed downward jet velocity `(0 0 -jetVelocity)`
- `bottom`: `noSlip`
- `outer`, `atmosphere`: `pressureInletOutletVelocity`
- `front`, `back`: `wedge`
- `axis`: `empty`

### `p_rgh`

- `bottom`: `fixedFluxPressure`
- `inlet`: `zeroGradient`
- `outer`, `atmosphere`: `prghTotalPressure`
- `front`, `back`: `wedge`
- `axis`: `empty`

## Time stepping and robustness defaults

- adaptive time stepping: `adjustTimeStep yes`
- `maxCo <= 0.5`
- `maxAlphaCo <= 0.25`
- bounded interface advection via `Gauss interfaceCompression vanLeer 1`
- PIMPLE with 2 outer correctors, 3 pressure correctors, and MULES correction

## Runtime monitoring

The case writes/monitors:

- residuals (`U`, `p_rgh`, `alpha.liquid`)
- inlet liquid flow rate
- liquid volume
- max velocity magnitude
- alpha min/max
- max shear rate
- min/max effective viscosity
- interface height on the axis (for cavity-depth tracking)
- Courant number and alpha extrema in the solver log

Post-process fields of interest:

- `alpha.liquid`
- `U`
- `p_rgh`
- `gammaDot.liquid`
- `muEff.liquid`
- `nu.liquid`

## Serial run

```bash
cd plungingJetHB
./Allrun
```

Behavior:

1. regenerate dictionaries from `constant/caseParameters`
2. attempt to build the custom viscosity library via `../build.sh`
3. run `blockMesh`
4. run `setFields`
5. run `foamRun`

If OpenFOAM utilities are not currently available in `PATH`, the scripts stop after rendering and print a clear message.

## Parallel run

```bash
cd plungingJetHB
./Allrun.parallel 8
```

The script rewrites `system/decomposeParDict` with the requested subdomain count, then runs:

- `blockMesh`
- `setFields`
- `decomposePar -force`
- `mpirun -np N foamRun -parallel`
- `reconstructPar`

## Post-processing helper

```bash
cd plungingJetHB
./Allpost
```

This creates:

- `postProcessing/paraview/loadCase.py`
- `postProcessing/summary.txt`

The ParaView helper requests the key interface, velocity, pressure, shear-rate, and viscosity fields.

## Sensitivity checks to perform

Before using the case for production data generation, check sensitivity to:

1. `meshLevel`
2. `domainRadius`
3. `gasHeight`
4. `poolDepth`
5. `wedgeAngle`
6. `maxCo`, `maxAlphaCo`, and `maxDeltaT`
7. rheology limits `muMin`, `muMax`, and `regularizationParameter`

## Limitations

- This is an **axisymmetric wedge approximation**, so it cannot capture true 3D asymmetric splashing, azimuthal breakup, or off-axis air entrainment structures.
- The custom rheology is applied to the liquid-phase laminar stress model; interface cells still rely on the native incompressibleVoF phase-fraction transport and mixture treatment.
- The current cavity-depth monitor is axis-based. If the cavity becomes strongly off-axis, add extra `interfaceHeight` sample locations.
- Because the inlet is represented on the top boundary, `jetStandOffDistance` is realized through the initial jet-column setup described above instead of an internal nozzle patch.

## Troubleshooting

### Divergence or pressure drift

- reduce `maxDeltaT`
- reduce `maxCo` and `maxAlphaCo`
- increase `domainRadius`
- check that the outer boundary is far enough from the jet core

### Unbounded `alpha`

- refine `meshLevel`
- reduce `maxAlphaCo`
- shorten `writeInterval` and inspect the interface near the inlet impact region

### Excessive viscosity / apparent stagnation

- check `muMax`
- check `regularizationParameter`
- verify `tauY`, `K`, and `n` are in SI units

### Wedge errors

- keep `wedgeAngle` between 2 and 5 degrees
- do not edit generated wedge faces manually; re-run `scripts/renderCase.py` or `Allmesh`

### Build issues for the custom rheology library

- source OpenFOAM Foundation v14 first
- run `../build.sh` from the case directory or `./build.sh` from the repository root
- if your site policies forbid user libraries, port `applications/regularisedHerschelBulkley/` into your local user application tree and rebuild there
