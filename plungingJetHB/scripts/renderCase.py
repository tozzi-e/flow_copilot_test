#!/usr/bin/env python3
from __future__ import annotations

import math
import re
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parents[1]
PARAM_FILE = CASE_DIR / "constant" / "caseParameters"

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  14
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
"""

PARAM_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s+([^;]+);")


def parse_parameters(path: Path) -> dict[str, float]:
    params: dict[str, float] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("//", 1)[0].strip()
        if not line:
            continue
        match = PARAM_RE.match(line)
        if match:
            key, raw_value = match.groups()
            try:
                params[key] = float(raw_value.strip())
            except ValueError:
                continue
    return params


def foam_dict(object_name: str, location: str, body: str, cls: str = "dictionary") -> str:
    return (
        HEADER
        + f"FoamFile\n{{\n    format      ascii;\n    class       {cls};\n    location    \"{location}\";\n    object      {object_name};\n}}\n"
        + "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
        + body.strip()
        + "\n\n// ************************************************************************* //\n"
    )


def fmt_scalar(value: float) -> str:
    return f"{value:.12g}"


def fmt_int(value: int) -> str:
    return str(int(value))


def fmt_point(pt: tuple[float, float, float]) -> str:
    return f"({fmt_scalar(pt[0])} {fmt_scalar(pt[1])} {fmt_scalar(pt[2])})"


def write(path: Path, content: str) -> None:
    path.write_text(content)


def main() -> None:
    p = parse_parameters(PARAM_FILE)
    required = [
        "jetDiameter", "jetVelocity", "jetStandOffDistance", "poolDepth",
        "gasHeight", "domainRadius", "wedgeAngle", "meshLevel", "rhoLiquid",
        "rhoAir", "muAir", "surfaceTension", "tauY", "K", "n",
        "regularizationParameter", "muMin", "muMax", "maxCo", "maxAlphaCo",
        "maxDeltaT", "endTime", "writeInterval"
    ]
    missing = [key for key in required if key not in p]
    if missing:
        raise SystemExit(f"Missing caseParameters entries: {', '.join(missing)}")

    D = p["jetDiameter"]
    inlet_radius = 0.5*D
    pool_depth = p["poolDepth"]*D
    gas_height = p["gasHeight"]*D
    total_height = pool_depth + gas_height
    stand_off = p["jetStandOffDistance"]*D
    outer_radius = p["domainRadius"]*D
    theta = math.radians(p["wedgeAngle"]/2.0)
    jet_tip_z = pool_depth + stand_off
    stub_start_z = min(jet_tip_z, total_height)
    eps = max(D, total_height)*1.0e-6

    if not (0 < p["wedgeAngle"] <= 5.0):
        raise SystemExit("wedgeAngle must lie in (0, 5] degrees")
    if gas_height < stand_off:
        raise SystemExit("gasHeight [D] must be >= jetStandOffDistance [D]")
    if outer_radius <= inlet_radius:
        raise SystemExit("domainRadius must be greater than 0.5D")
    if p["meshLevel"] < 4:
        raise SystemExit("meshLevel must be >= 4")
    if p["muMin"] <= 0 or p["muMax"] < p["muMin"]:
        raise SystemExit("Require 0 < muMin <= muMax")
    if p["regularizationParameter"] <= 0:
        raise SystemExit("regularizationParameter must be positive")

    def wedge_point(r: float, z: float, sign: float) -> tuple[float, float, float]:
        return (r*math.cos(theta), sign*r*math.sin(theta), z)

    vertices = [
        wedge_point(0.0, 0.0, -1.0),
        wedge_point(inlet_radius, 0.0, -1.0),
        wedge_point(outer_radius, 0.0, -1.0),
        wedge_point(0.0, pool_depth, -1.0),
        wedge_point(inlet_radius, pool_depth, -1.0),
        wedge_point(outer_radius, pool_depth, -1.0),
        wedge_point(0.0, total_height, -1.0),
        wedge_point(inlet_radius, total_height, -1.0),
        wedge_point(outer_radius, total_height, -1.0),
        wedge_point(0.0, 0.0, 1.0),
        wedge_point(inlet_radius, 0.0, 1.0),
        wedge_point(outer_radius, 0.0, 1.0),
        wedge_point(0.0, pool_depth, 1.0),
        wedge_point(inlet_radius, pool_depth, 1.0),
        wedge_point(outer_radius, pool_depth, 1.0),
        wedge_point(0.0, total_height, 1.0),
        wedge_point(inlet_radius, total_height, 1.0),
        wedge_point(outer_radius, total_height, 1.0),
    ]

    nr_inner = max(6, round(0.5*p["meshLevel"]))
    nr_outer = max(10, round((p["domainRadius"] - 0.5)*p["meshLevel"]))
    nz_pool = max(12, round(p["poolDepth"]*p["meshLevel"]))
    nz_gas = max(12, round(p["gasHeight"]*p["meshLevel"]))

    write(
        CASE_DIR / "system" / "blockMeshDict",
        foam_dict(
            "blockMeshDict",
            "system",
            f"""units [m];

vertices
(
{chr(10).join(f'    {fmt_point(v)}' for v in vertices)}
);

blocks
(
    hex (0 1 4 3 9 10 13 12) ({fmt_int(nr_inner)} {fmt_int(nz_pool)} 1) simpleGrading (1 1 1)
    hex (1 2 5 4 10 11 14 13) ({fmt_int(nr_outer)} {fmt_int(nz_pool)} 1) simpleGrading (1.2 1 1)
    hex (3 4 7 6 12 13 16 15) ({fmt_int(nr_inner)} {fmt_int(nz_gas)} 1) simpleGrading (1 1 1)
    hex (4 5 8 7 13 14 17 16) ({fmt_int(nr_outer)} {fmt_int(nz_gas)} 1) simpleGrading (1.2 1 1)
);

boundary
(
    axis
    {{
        type empty;
        faces
        (
            (0 3 12 9)
            (3 6 15 12)
        );
    }}

    inlet
    {{
        type patch;
        faces
        (
            (6 7 16 15)
        );
    }}

    bottom
    {{
        type wall;
        faces
        (
            (0 1 10 9)
            (1 2 11 10)
        );
    }}

    outer
    {{
        type patch;
        faces
        (
            (2 11 14 5)
            (5 14 17 8)
        );
    }}

    atmosphere
    {{
        type patch;
        faces
        (
            (7 8 17 16)
        );
    }}

    front
    {{
        type wedge;
        faces
        (
            (0 1 4 3)
            (1 2 5 4)
            (3 4 7 6)
            (4 5 8 7)
        );
    }}

    back
    {{
        type wedge;
        faces
        (
            (9 12 13 10)
            (10 13 14 11)
            (12 15 16 13)
            (13 16 17 14)
        );
    }}
);
"""
        )
    )

    write(
        CASE_DIR / "system" / "setFieldsDict",
        foam_dict(
            "setFieldsDict",
            "system",
            f"""defaultFieldValues
(
    volScalarFieldValue alpha.liquid 0
);

regions
(
    boxToCell
    {{
        box ({fmt_scalar(-outer_radius)} {fmt_scalar(-outer_radius)} {fmt_scalar(-eps)}) ({fmt_scalar(outer_radius)} {fmt_scalar(outer_radius)} {fmt_scalar(pool_depth)});
        fieldValues
        (
            volScalarFieldValue alpha.liquid 1
        );
    }}

    cylinderToCell
    {{
        p1 (0 0 {fmt_scalar(stub_start_z)});
        p2 (0 0 {fmt_scalar(total_height + eps)});
        radius {fmt_scalar(inlet_radius)};
        fieldValues
        (
            volScalarFieldValue alpha.liquid 1
        );
    }}
);
"""
        )
    )

    write(
        CASE_DIR / "constant" / "phaseProperties",
        foam_dict(
            "phaseProperties",
            "constant",
            f"""phases          (liquid air);
sigma           {fmt_scalar(p['surfaceTension'])};"""
        )
    )

    write(
        CASE_DIR / "constant" / "physicalProperties.liquid",
        foam_dict(
            "physicalProperties.liquid",
            "constant",
            f"""viscosityModel  constant;
nu              {fmt_scalar(p['muMin']/p['rhoLiquid'])};
rho             {fmt_scalar(p['rhoLiquid'])};"""
        )
    )

    write(
        CASE_DIR / "constant" / "physicalProperties.air",
        foam_dict(
            "physicalProperties.air",
            "constant",
            f"""viscosityModel  constant;
nu              {fmt_scalar(p['muAir']/p['rhoAir'])};
rho             {fmt_scalar(p['rhoAir'])};"""
        )
    )

    write(
        CASE_DIR / "constant" / "momentumTransport",
        foam_dict(
            "momentumTransport",
            "constant",
            "simulationType  twoPhaseTransport;"
        )
    )

    write(
        CASE_DIR / "constant" / "momentumTransport.air",
        foam_dict(
            "momentumTransport.air",
            "constant",
            "simulationType  laminar;"
        )
    )

    write(
        CASE_DIR / "constant" / "momentumTransport.liquid",
        foam_dict(
            "momentumTransport.liquid",
            "constant",
            f"""simulationType  laminar;

laminar
{{
    model                   generalisedNewtonian;
    viscosityModel          regularisedHerschelBulkley;

    regularisedHerschelBulkleyCoeffs
    {{
        rho                     {fmt_scalar(p['rhoLiquid'])};
        tauY                    {fmt_scalar(p['tauY'])};
        K                       {fmt_scalar(p['K'])};
        n                       {fmt_scalar(p['n'])};
        regularizationParameter {fmt_scalar(p['regularizationParameter'])};
        muMin                   {fmt_scalar(p['muMin'])};
        muMax                   {fmt_scalar(p['muMax'])};
    }}
}}"""
        )
    )

    write(
        CASE_DIR / "constant" / "transportProperties",
        foam_dict(
            "transportProperties",
            "constant",
            """// Informational compatibility file for users migrating from older cases.
// Foundation v14 incompressibleVoF reads phaseProperties,
// physicalProperties.<phase>, and momentumTransport.<phase> instead.
includeFiles
(
    phaseProperties
    physicalProperties.liquid
    physicalProperties.air
    momentumTransport
    momentumTransport.liquid
    momentumTransport.air
);"""
        )
    )

    write(
        CASE_DIR / "constant" / "turbulenceProperties",
        foam_dict(
            "turbulenceProperties",
            "constant",
            """// Compatibility placeholder only.
// Foundation v14 uses the momentumTransport dictionaries.
simulationType  laminar;"""
        )
    )

    write(
        CASE_DIR / "system" / "controlDict",
        foam_dict(
            "controlDict",
            "system",
            f"""libs
(
    \"libutilityFunctionObjects.so\"
    \"libfieldFunctionObjects.so\"
    \"libregularisedHBViscosityModel.so\"
);

solver          incompressibleVoF;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {fmt_scalar(p['endTime'])};
deltaT          1e-5;

writeControl    adjustableRunTime;
writeInterval   {fmt_scalar(p['writeInterval'])};
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable yes;

adjustTimeStep  yes;
maxCo           {fmt_scalar(p['maxCo'])};
maxAlphaCo      {fmt_scalar(p['maxAlphaCo'])};
maxDeltaT       {fmt_scalar(p['maxDeltaT'])};

functions
{{
    residuals
    {{
        type            residuals;
        libs            (\"libutilityFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        fields          (U p_rgh alpha.liquid);
    }}

    inletFlowRate
    {{
        type            surfaceFieldValue;
        libs            (\"libfieldFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     false;
        patch           inlet;
        operation       sum;
        fields          (phi);
        weightField     alpha.liquid;
    }}

    liquidVolume
    {{
        type            volFieldValue;
        libs            (\"libfieldFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     false;
        cellZone        all;
        operation       volIntegrate;
        fields          (alpha.liquid);
    }}

    maxVelocity
    {{
        type            volFieldValue;
        libs            (\"libfieldFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     false;
        cellZone        all;
        operation       maxMag;
        fields          (U);
    }}

    alphaMin
    {{
        type            volFieldValue;
        libs            (\"libfieldFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     false;
        cellZone        all;
        operation       min;
        fields          (alpha.liquid);
    }}

    alphaMax
    {{
        type            volFieldValue;
        libs            (\"libfieldFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     false;
        cellZone        all;
        operation       max;
        fields          (alpha.liquid);
    }}

    gammaDotMax
    {{
        type            volFieldValue;
        libs            (\"libfieldFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     false;
        cellZone        all;
        operation       max;
        fields          (gammaDot.liquid);
    }}

    muEffMin
    {{
        type            volFieldValue;
        libs            (\"libfieldFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     false;
        cellZone        all;
        operation       min;
        fields          (muEff.liquid);
    }}

    muEffMax
    {{
        type            volFieldValue;
        libs            (\"libfieldFunctionObjects.so\");
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     false;
        cellZone        all;
        operation       max;
        fields          (muEff.liquid);
    }}

    cavityDepthAxis
    {{
        type            interfaceHeight;
        libs            (\"libfieldFunctionObjects.so\");
        alpha           alpha.liquid;
        locations       ((0 0 {fmt_scalar(total_height)}));
    }}
}}"""
        )
    )

    u_path = CASE_DIR / "0" / "U"
    write(
        u_path,
        foam_dict(
            "U",
            "0",
            f"""dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);

boundaryField
{{
    axis
    {{
        type            empty;
    }}

    inlet
    {{
        type            fixedValue;
        value           uniform (0 0 {fmt_scalar(-p['jetVelocity'])});
    }}

    bottom
    {{
        type            noSlip;
    }}

    outer
    {{
        type            pressureInletOutletVelocity;
        value           uniform (0 0 0);
    }}

    atmosphere
    {{
        type            pressureInletOutletVelocity;
        value           uniform (0 0 0);
    }}

    front
    {{
        type            wedge;
    }}

    back
    {{
        type            wedge;
    }}
}}""",
            cls="volVectorField"
        )
    )

    summary = (
        f"Rendered plungingJetHB for D={fmt_scalar(D)} m, Uj={fmt_scalar(p['jetVelocity'])} m/s, "
        f"R={fmt_scalar(outer_radius)} m, H={fmt_scalar(total_height)} m, wedge={fmt_scalar(p['wedgeAngle'])} deg.\n"
        f"Initial pool depth = {fmt_scalar(pool_depth)} m, gas height = {fmt_scalar(gas_height)} m, jet tip stand-off = {fmt_scalar(stand_off)} m.\n"
        f"HB limits: mu in [{fmt_scalar(p['muMin'])}, {fmt_scalar(p['muMax'])}] Pa.s, gamma_reg >= {fmt_scalar(p['regularizationParameter'])} 1/s.\n"
    )
    write(CASE_DIR / "renderedCaseSummary.txt", summary)
    print(summary, end="")


if __name__ == "__main__":
    main()
