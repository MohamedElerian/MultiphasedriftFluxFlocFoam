Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
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

multiphaseDriftMixture mixture(U, phi);

// Mixture density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);
rho.oldTime();

// Construct incompressible turbulence model
// autoPtr<incompressible::turbulenceModel> turbulence
// (
//     incompressible::turbulenceModel::New(U, phi, mixture)
// );
const surfaceScalarField& rhoPhi = mixture.rhoPhi();

autoPtr
<
    CompressibleTurbulenceModel<multiphaseDriftMixture>
> turbulence
(
    CompressibleTurbulenceModel<multiphaseDriftMixture>::New(rho, U, rhoPhi, mixture)
);

// Set turbulence in the mixture
mixture.setTurbulence(turbulence());

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);



if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
