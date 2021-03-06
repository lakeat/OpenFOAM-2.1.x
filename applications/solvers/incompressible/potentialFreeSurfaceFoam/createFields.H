    Info<< "Reading field p (kinematic)\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    #include "readGravitationalAcceleration.H"

    Info<< "\nReading freeSurfaceProperties\n" << endl;

    IOdictionary freeSurfaceProperties
    (
        IOobject
        (
            "freeSurfaceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word freeSurfacePatch(freeSurfaceProperties.lookup("freeSurfacePatch"));
    label freeSurfacePatchI = mesh.boundaryMesh().findPatchID(freeSurfacePatch);
    if (freeSurfacePatchI < 0)
    {
        FatalErrorIn(args.executable())
            << "Patch " << freeSurfacePatch << " not found. "
            << "Available patches are:" << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    Info<< "Creating field refLevel\n" << endl;
    volVectorField refLevel
    (
        IOobject
        (
            "refLevel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    refLevel.boundaryField()[freeSurfacePatchI]
        == mesh.C().boundaryField()[freeSurfacePatchI];

    Info<< "Creating field zeta\n" << endl;
    volVectorField zeta
    (
        IOobject
        (
            "zeta",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    Info<< "Creating field p_gh\n" << endl;
    volScalarField p_gh
    (
        IOobject
        (
            "p_gh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Force p_gh to be consistent with p
    // Height is made relative to field 'refLevel'
    p_gh =  p - (g & (mesh.C() + zeta - refLevel));


    label p_ghRefCell = 0;
    scalar p_ghRefValue = 0.0;
    setRefCell(p_gh, pimple.dict(), p_ghRefCell, p_ghRefValue);


    IObasicSourceList sources(mesh);
