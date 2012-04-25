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

#include "V2F_original.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{
#include "/home/rom/OpenFOAM/OpenFOAM-1.6/applications/solvers/incompressible/boundaryFoam2/scalars.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(V2F_original, 0);
addToRunTimeSelectionTable(RASModel, V2F_original, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> V2F_original::T() const
{
return min ( 
             max(
             (k_+k0_)/(epsilon_+epsilonSmall_) ,CT*sqrt(nu()/(epsilon_ +epsilonSmall_))
                 ),
            alpha*(k_+k0_)/( sqrt(6.0)*Cmu*(mag(symm(fvc::grad(U_)))+f0_)*(v2_+k0_) )
            );

/*return max ((k_+k0_)/(epsilon_+epsilonSmall_), CT*sqrt(nu()/(epsilon_+epsilonSmall_)));*/

}

tmp<volScalarField> V2F_original::L() const
{
    return CL*max(
                      min(
                              pow((k_+k0_),1.5)/(epsilon_ +epsilonSmall_),
                              pow((k_+k0_),1.5)/(sqrt(6.0)*Cmu*(mag(symm(fvc::grad(U_)))+f0_)*(v2_+k0_))
                          ),
                      CEta*pow(pow(nu(),3)/(epsilon_+epsilonSmall_),0.25)
                  );

/*return CL*max( pow(k_,1.5)/epsilon_ ,CEta*pow(pow(nu(),3)/epsilon_,0.25));*/

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
V2F_original::V2F_original
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
            0.19
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
            1.40
        )
    ),
    C2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.3
        )
    ),
    CL
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.3
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

    Sigmav2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmav2",
            coeffDict_,
            1.0
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

    nut_(Cmu*v2_*k_*T())

{

    printCoeffs();
}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> V2F_original::R() const
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


tmp<volSymmTensorField> V2F_original::devReff() const
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


tmp<fvVectorMatrix> V2F_original::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool V2F_original::read()
{
    if (RASModel::read())
    {

	Cmu.readIfPresent(coeffDict());
        Ceps2.readIfPresent(coeffDict_);
        C1.readIfPresent(coeffDict_);
        C2.readIfPresent(coeffDict_);
        CL.readIfPresent(coeffDict_);
        CEta.readIfPresent(coeffDict_);
        CT.readIfPresent(coeffDict_);
        SigmaK.readIfPresent(coeffDict_);
        SigmaEps.readIfPresent(coeffDict_);
	Sigmav2.readIfPresent(coeffDict_);
        alpha.readIfPresent(coeffDict_);
        
        return true;
    }
    else
    {
        return false;
    }
}


void V2F_original::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    RASModel::correct();

    //volScalarField y=mesh_.C().component(vector::Y);

    volScalarField S2 = 2*magSqr(symm(fvc::grad(U_)));

    volScalarField G = nut_*S2;

    volScalarField T_ = T();
    volScalarField L_ = L();
    volScalarField Ceps1 = 1.4*(1+0.05*sqrt(k_/(v2_+k0_)));
   // Ceps1 = 1.3+0.25/(1+pow(y/(2*L_/CL),8));
    
    //bound f and k for stability
    dimensionedScalar k1_("k1", dimVelocity*dimVelocity, -100.0); //upper bound for k
    dimensionedScalar f1_("f1small", dimless/dimTime, -50000.0);  //upper and lower bound for f 

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
    k_=-k_;
    bound(k_, k1_);
    k_=-k_;

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

    bound(f_, f1_);
    // v2 equation
    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
      - fvm::laplacian(Dv2Eff(), v2_)
     ==
	k_*f_
      - fvm::Sp(1.0/T_, v2_)
    );
    v2Eqn().relax();
    solve(v2Eqn);
    bound(v2_, k0_);


    // f equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(L_),f_)
      - ((1.0-C1)*(2/3-v2_/(k_+k0_))/T_-C2*G/(k_+k0_))/sqr(L_)
    );
    fEqn().relax();
#   include "FWallI.H"
#   include "wallF.H"
    solve(fEqn);
    bound(f_, f1_);
    f_=-f_;
    bound(f_, f1_);
    f_=-f_;



    nut_=Cmu*v2_*k_*T();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
