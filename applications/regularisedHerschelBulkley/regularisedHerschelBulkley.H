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

#ifndef regularisedHerschelBulkley_H
#define regularisedHerschelBulkley_H

#include "strainRateViscosityModel.H"

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{

class regularisedHerschelBulkley
:
    public strainRateViscosityModel
{
    dimensionedScalar rho_;
    dimensionedScalar n_;
    dimensionedScalar K_;
    dimensionedScalar tauY_;
    dimensionedScalar regularizationParameter_;
    dimensionedScalar muMin_;
    dimensionedScalar muMax_;

    mutable volScalarField gammaDot_;
    mutable volScalarField muEff_;

public:

    TypeName("regularisedHerschelBulkley");

    regularisedHerschelBulkley
    (
        const dictionary& viscosityProperties,
        const Foam::viscosity& viscosity,
        const volVectorField& U
    );

    virtual ~regularisedHerschelBulkley()
    {}

    virtual bool read(const dictionary& viscosityProperties);

    virtual tmp<volScalarField> nu
    (
        const volScalarField& nu0,
        const volScalarField& strainRate
    ) const;
};

} // End namespace generalisedNewtonianViscosityModels
} // End namespace laminarModels
} // End namespace Foam

#endif

// ************************************************************************* //
