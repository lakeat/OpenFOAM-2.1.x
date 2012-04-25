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

#include "zetaF.H"
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

defineTypeNameAndDebug(zetaF, 0);
addToRunTimeSelectionTable(RASModel, zetaF, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> zetaF::T() const
{
return min ( 
             max(
             k_/epsilon_ ,CT*sqrt(nu()/epsilon_ )
                 ),
            alpha/( sqrt(6.0)*Cmu*(mag(symm(fvc::grad(U_)))+f0_)*zeta_ )
            );

//return max (k_/epsilon_ , CT*sqrt(nu()/epsilon_));

}

tmp<volScalarField> zetaF::L() const
{
    return CL*max(
                      min(
                              pow(k_,1.5)/epsilon_ ,
                              sqrt(k_)/(sqrt(6.0)*Cmu*(mag(symm(fvc::grad(U_)))+f0_)*zeta_)
                          ),
                      CEta*pow(pow(nu(),3)/epsilon_,0.25)
                  );

//return CL*max( pow(k_,1.5)/epsilon_ ,CEta*pow(pow(nu(),3)/epsilon_,0.25));

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
zetaF::zetaF
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

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
    
    Ceps2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.9
        )
    ),
    
    C1
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            0.40
        )
    ),
    C2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.65
        )
    ),
    CL
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.36
        )
    ),
    CEta
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            coeffDict_,
            85.0
        )
    ),
    CT
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT",
            coeffDict_,
            6.0
        )
    ),
    
    SigmaK
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "SigmaK",
            coeffDict_,
            1.0
        )
    ),
    
    SigmaEps
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "SigmaEps",
            coeffDict_,
            1.3
        )
    ),

    SigmaZeta
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "SigmaZeta",
            coeffDict_,
            1.2
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

    zeta_
    (
        IOobject
        (
            "zeta",
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


    nut_(min(CmuKE*sqr(k_)/(epsilon_ + epsilonSmall_), Cmu*zeta_*k_*T()))
    //nut_(Cmu*zeta_*k_*T())

{

    printCoeffs();
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> zetaF::R() const
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


tmp<volSymmTensorField> zetaF::devReff() const
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


tmp<fvVectorMatrix> zetaF::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool zetaF::read()
{
    if (RASModel::read())
    {
        Cmu.readIfPresent(coeffDict_);
	CmuKE.readIfPresent(coeffDict_);
	Cmu.readIfPresent(coeffDict_);
        Ceps2.readIfPresent(coeffDict_);
        C1.readIfPresent(coeffDict_);
        C2.readIfPresent(coeffDict_);
        CL.readIfPresent(coeffDict_);
        CEta.readIfPresent(coeffDict_);
        CT.readIfPresent(coeffDict_);
        SigmaK.readIfPresent(coeffDict_);
        SigmaEps.readIfPresent(coeffDict_);
	SigmaZeta.readIfPresent(coeffDict_);
        alpha.readIfPresent(coeffDict_);
        
        return true;
    }
    else
    {
        return false;
    }
}


void zetaF::correct()
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
    volScalarField Ceps1 = 1.4*(1.0+0.012/zeta_);

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
    bound(k_, k0_);

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
#   include "epsilonWallI.H"
#   include "wallDissipation.H"
    solve(epsEqn);
    bound(epsilon_, epsilon0_);

    // zeta equation
    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(zeta_)
      + fvm::div(phi_, zeta_)
      - fvm::laplacian(DzetaEff(), zeta_)
     ==
	min(-(1.0/T_)*(C1+(C2*G/epsilon_))*(zeta_-(2.0/3.0)),f_)
      - fvm::Sp(G/k_, zeta_)
    );
    zetaEqn().relax();
    solve(zetaEqn);
    bound(zeta_, 1e-20);
    /*zeta_=-zeta_;
    bound(zeta_, -0.666666666);
    zeta_=-zeta_;*/




    // f equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(L_),f_)
      - (1.0/T_)*(C1+(C2*G/epsilon_))*(zeta_-(2.0/3.0))*(1.0/sqr(L_))
    );
    fEqn().relax();
#   include "FWallI.H"
#   include "wallF.H"
    solve(fEqn);


   // nut_=Cmu*zeta_*k_*T();
    nut_=min(CmuKE*sqr(k_)/(epsilon_ + epsilonSmall_), Cmu*zeta_*k_*T());


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
