/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson

Description

SourceFiles
    JohnsonJacksonFrictionalStress.C

\*---------------------------------------------------------------------------*/

#ifndef JohnsonJackson_H
#define JohnsonJackson_H

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{

/*---------------------------------------------------------------------------*\
               Class JohnsonJackson Declaration
\*---------------------------------------------------------------------------*/

class JohnsonJackson
:
    public frictionalStressModel
{

public:

    //- Runtime type information
    TypeName("JohnsonJackson");


    // Constructors

        //- Construct from components
        JohnsonJackson(const dictionary& dict);


    //- Destructor
    virtual ~JohnsonJackson();


    // Member functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& eta,
            const dimensionedScalar& p
        ) const;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& eta,
            const dimensionedScalar& p
        ) const;

        virtual tmp<volScalarField> muf
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax,
            const volScalarField& pf,
            const volSymmTensorField& D,
            const dimensionedScalar& phi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionalStressModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
