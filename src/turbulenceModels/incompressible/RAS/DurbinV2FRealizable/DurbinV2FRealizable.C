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

\*---------------------------------------------------------------------------*/

#include "DurbinV2FRealizable.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DurbinV2FRealizable, 0);
addToRunTimeSelectionTable(RASModel, DurbinV2FRealizable, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> DurbinV2FRealizable::T() const
{
return min ( max(k_/epsilon_ ,6.0*sqrt(nu()/epsilon_ )),0.6*k_/(sqrt(6.0)*Cmu*v2_*mag(symm(fvc::grad(U_)))+VSMALL ));
//return max(k_/(epsilon_ + VSMALL),6.0*sqrt(nu()/(epsilon_ + VSMALL)));
}

tmp<volScalarField> DurbinV2FRealizable::L() const
{
return CL*max( min(pow(k_,1.5)/epsilon_ ,pow(k_,1.5)/(sqrt(6.0)*Cmu*v2_*mag(symm(fvc::grad(U_)))+VSMALL)),CEta*pow(pow(nu(),3.0)/epsilon_ ,0.25));
//return  CL*max(pow(k_,1.5)/(epsilon_ + VSMALL),CEta*pow(pow(nu(),3.0)/(epsilon_ + VSMALL),0.25));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
DurbinV2FRealizable::DurbinV2FRealizable
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.22
        )
    ),
    
    
    CmuKE
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuKE",
            coeffDict_,
            0.09
        )
    ),
    Ceps10
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps10",
            coeffDict_,
            1.40
        )
    ),
    
    Ceps11
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps11",
            coeffDict_,
            0.05
        )
    ),
    
    Ceps2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.90
        )
    ),
    
    C1
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.40
        )
    ),
    C2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.30
        )
    ),
    CL
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.23
        )
    ),
    CEta
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            coeffDict_,
            70.0
        )
    ),
    
    oneOnSigmaK
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "oneOnSigmaK",
            coeffDict_,
            1.0
        )
    ),
    
    oneOnSigmaEps
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "oneOnSigmaEps",
            coeffDict_,
            0.76923
        )
    ),

    alpha
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha",
            coeffDict_,
            0.6
        )
    ),

    f0_("f0small", dimless/dimTime, SMALL),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    v2_
    (
        IOobject
        (
            "v2",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),


    // Calculate viscosity (with Davidson correction - 2003)
    nut_(min(CmuKE*sqr(k_)/(epsilon_ + VSMALL), Cmu*v2_*T()))

{

    printCoeffs();
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> DurbinV2FRealizable::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*2*symm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> DurbinV2FRealizable::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> DurbinV2FRealizable::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool DurbinV2FRealizable::read()
{
    if (RASModel::read())
    {
        Cmu.readIfPresent(coeffDict_);
        CmuKE.readIfPresent(coeffDict_);
        Ceps10.readIfPresent(coeffDict_);
        Ceps11.readIfPresent(coeffDict_);
        Ceps2.readIfPresent(coeffDict_);
        C1.readIfPresent(coeffDict_);
        C2.readIfPresent(coeffDict_);
        CL.readIfPresent(coeffDict_);
        CEta.readIfPresent(coeffDict_);
        oneOnSigmaK.readIfPresent(coeffDict_);
        oneOnSigmaEps.readIfPresent(coeffDict_);
        alpha.readIfPresent(coeffDict_);
        
        return true;
    }
    else
    {
        return false;
    }
}


void DurbinV2FRealizable::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    RASModel::correct();

    volScalarField S2 = 2*magSqr(symm(fvc::grad(U_)));

    volScalarField G = nut_*S2;

    volScalarField T_ = T();
    volScalarField L_ = L();
    volScalarField Ceps1 = Ceps10*(scalar(1.0)+Ceps11*sqrt(k_/v2_));

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp(1.0/T_, k_)
    );
    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);

    // Dissipation rate equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        Ceps1*G/T_
      - fvm::Sp(Ceps2/T_, epsilon_)
    );
    epsEqn().relax();
#   include "epsilonV2FWallI.H"
#   include "wallDissipationV2FI.H"
    solve(epsEqn);
    bound(epsilon_, epsilonMin_);

    // v2 equation
    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
      - fvm::laplacian(DkEff(), v2_)
     ==
    // Davidson correction - 2003
        min(k_*f_,-((C1-scalar(6.0))*v2_ - 2.0/3.0*k_*(C1-scalar(1.0)))/T_ + C2*G)
      - fvm::Sp(6.0/T_, v2_)
    );
    v2Eqn().relax();
    solve(v2Eqn);
    bound(v2_, kMin_);

    // f equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(L_),f_)
      - ((C1-scalar(6.0))*v2_/k_ - 2.0/3.0*(C1-scalar(1.0)))/(sqr(L_)*T_)
      + C2*G/(k_*sqr(L_))
    );
    fEqn().relax();
    solve(fEqn);
    bound(f_, f0_);


    // Re-calculate viscosity (with Davidson correction - 2003)
    nut_ = min(CmuKE*sqr(k_)/(epsilon_ + VSMALL),Cmu*v2_*T());

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
