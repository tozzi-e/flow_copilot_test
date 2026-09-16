/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenAI
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
\*---------------------------------------------------------------------------*/

#include "regularisedHerschelBulkley.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{
    defineTypeNameAndDebug(regularisedHerschelBulkley, 0);

    addToRunTimeSelectionTable
    (
        generalisedNewtonianViscosityModel,
        regularisedHerschelBulkley,
        dictionary
    );
}
}
}

Foam::laminarModels::generalisedNewtonianViscosityModels::regularisedHerschelBulkley::
regularisedHerschelBulkley
(
    const dictionary& viscosityProperties,
    const Foam::viscosity& viscosity,
    const volVectorField& U
)
:
    strainRateViscosityModel(viscosityProperties, viscosity, U),
    rho_("rho", dimDensity, 0),
    n_("n", dimless, 1),
    K_("K", dimDynamicViscosity*pow(dimTime, n_ - 1), 0),
    tauY_("tauY", dimPressure, 0),
    regularizationParameter_("regularizationParameter", dimless/dimTime, 0),
    muMin_("muMin", dimDynamicViscosity, 0),
    muMax_("muMax", dimDynamicViscosity, 0),
    gammaDot_
    (
        IOobject
        (
            IOobject::groupName(typedName("gammaDot"), U.group()),
            U.time().name(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimless/dimTime, 0)
    ),
    muEff_
    (
        IOobject
        (
            IOobject::groupName(typedName("muEff"), U.group()),
            U.time().name(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimDynamicViscosity, 0)
    )
{
    read(viscosityProperties);
    correct();
}

bool Foam::laminarModels::generalisedNewtonianViscosityModels::
regularisedHerschelBulkley::read
(
    const dictionary& viscosityProperties
)
{
    strainRateViscosityModel::read(viscosityProperties);

    const dictionary& coeffs =
        viscosityProperties.optionalTypeDict(typeName);

    rho_.read(coeffs);
    n_.read(coeffs);
    K_.dimensions().reset(dimDynamicViscosity*pow(dimTime, n_ - 1));
    K_.read(coeffs);
    tauY_.read(coeffs);
    regularizationParameter_.read(coeffs);
    muMin_.read(coeffs);
    muMax_.read(coeffs);

    return true;
}

Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalisedNewtonianViscosityModels::
regularisedHerschelBulkley::nu
(
    const volScalarField&,
    const volScalarField& strainRate
) const
{
    gammaDot_ = max
    (
        strainRate,
        dimensionedScalar(dimless/dimTime, regularizationParameter_.value())
    );

    muEff_ =
        tauY_/gammaDot_
      + K_*pow(gammaDot_, n_ - scalar(1));

    muEff_ = max(muEff_, muMin_);
    muEff_ = min(muEff_, muMax_);

    return muEff_/rho_;
}

// ************************************************************************* //
