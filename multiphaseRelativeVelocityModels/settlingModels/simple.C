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

\*---------------------------------------------------------------------------*/

#include "simple.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(simple, 0);
    addToRunTimeSelectionTable(settlingModel, simple, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModels::simple::simple
(
    const terminalVelocitySettlingModel& UInf,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    settlingModel(UInf, dict, mesh),
    alphatMax_("alphatMax", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModels::simple::~simple()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::settlingModels::simple::Vk
(
    const phaseDrift& phasec,
    const phaseDrift& phasek
) const
{
    tmp<volScalarField> tVk
    (
        new volScalarField
        (
            IOobject
            (
                "Vk",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("zero", dimless, Zero)
        )
    );
    volScalarField& Vk = tVk.ref();

    const volScalarField& nu = phasec.nu()();

    const volScalarField& V0(mag(UInfModel().V0(phasec, phasek)));

    const volScalarField& alphat  =
        mesh_.lookupObject<volScalarField>("alpha.solids");
    //const volScalarField alphat(1.0 - phasec);

    forAll (Vk, i)
    {
        if (alphat[i] < alphatMax_.value())
        {
            const scalar RepInf34
            (
                pow
                (
                    V0[i]
                  * phasek.d().value()/nu[i],
                    scalar(0.75)
                )
            );

            const scalar nuk
            (
                2.35*(2 + 0.175*RepInf34)/(1 + 0.175*RepInf34) - 1.0
            );

            Vk[i] = pow((1- alphat[i]), nuk);

        }
    }

    return tVk;

}


// ************************************************************************* //
