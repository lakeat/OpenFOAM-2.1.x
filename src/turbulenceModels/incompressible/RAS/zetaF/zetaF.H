/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    zetaF

Description
    Durbin zetaF turbulence model.
    Modified by Daniel Wei, @ NatHaz, 20120425

SourceFiles
    zetaF.C

\*---------------------------------------------------------------------------*/

#ifndef zetaF_H
#define zetaF_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class zetaF Declaration
\*---------------------------------------------------------------------------*/

class zetaF
:
    public RASModel
{
    // Private data

        dimensionedScalar Cmu;
        dimensionedScalar CmuKE;
        dimensionedScalar Ceps2;
        dimensionedScalar C1;
        dimensionedScalar C2;
        dimensionedScalar CL;
        dimensionedScalar CEta;
        dimensionedScalar CT;	
        dimensionedScalar SigmaK;
        dimensionedScalar SigmaEps;
        dimensionedScalar SigmaZeta;
        dimensionedScalar alpha;

        dimensionedScalar f0_;

        volScalarField k_;
        volScalarField epsilon_;
        volScalarField zeta_;
        volScalarField f_;
        volScalarField nut_;


    // Private member functions

        tmp<volScalarField> T() const;
        tmp<volScalarField> L() const;

public:

    //- Runtime type information
    TypeName("zetaF");

    // Constructors

        //- from components
        zetaF
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    // Destructor

        ~zetaF(){}


    // Member Functions

        tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DkEff", nut_ / SigmaK + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DepsilonEff",nut_/ SigmaEps + nu())
            );
        }

        //- Return the effective diffusivity for zeta
        tmp<volScalarField> DzetaEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DzetaEff", nut_/SigmaZeta + nu())
            );
        }

        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Note that epsilonTilda is returned as epsilon.
        //  This is the appropriate variable for most purposes.
        tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        tmp<volScalarField> zeta() const
        {
            return zeta_;
        }

        tmp<volScalarField> f() const
        {
            return f_;
        }

        tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
