/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

Application
    mulitphaseDriftFluxFoam

Group
    grpMultiphaseSolvers

Description
    Solver for n incompressible fluids using the mixture approach with the
    drift-flux approximation for relative motion of the phases.

    The mixture is composed of a continous phase and n phases with
    corresponding density and diameter for each one.

    Used for simulating the settling of the dispersed phase and other similar
    separation problems.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiphaseDriftMixture.H"
//#include "multiphaseRelativeVelocityModel.H"
#include "turbulenceModel.H"
//#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
//#include "localMin.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "CourantNo.H"
    
    //#include "alphaCourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
       // #include "alphaCourantNo.H"
        #include "CourantNo.H" //
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mixture.correctUkm();   // COMPUTE THE DIFFUSION VELOCITY  (DISPERSED PHASE RELATIVE TO MIXTURE)
        mixture.solve();       //  SOLVE THE CONTINUITY EQUATION; UPDATE THE ALPHA
        rho = mixture.rho();   //  COMPUTE MIXTURE VELOCITY
        mixture.correct();     //  UPDATE VISCOSITY

       

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = mixture.rho();

        mixture.correct();

        runTime.write();
        //mixture.writePhase(); 

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
