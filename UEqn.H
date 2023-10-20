MRF.correctBoundaryVelocity(U);

const surfaceScalarField phiRho(fvc::interpolate(rho)*phi);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(phiRho, U)
  + MRF.DDt(rho, U)
  + fvc::div(mixture.tauDm(), "div(tauDm)")
  + turbulence->divDevRhoReff(U)
    ==
    fvOptions(rho, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
        ==
        fvc::reconstruct
        (
            (
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
