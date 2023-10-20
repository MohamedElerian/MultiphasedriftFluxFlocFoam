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

#include "multiphaseRelativeVelocityModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseRelativeVelocityModel::multiphaseRelativeVelocityModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    UInfModel_
    (
        terminalVelocitySettlingModel::New
        (
            dict,
            mesh
        )
    ),
    settlingModel_
    (
        settlingModel::New
        (
            UInfModel_(),
            dict,
            mesh
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseRelativeVelocityModel::~multiphaseRelativeVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::volVectorField>
Foam::multiphaseRelativeVelocityModel::Ukr
(
    const phaseDrift& phaseCont,
    const phaseDrift& phasek
) const
{
    return
    (
        UInfModel_->V0(phaseCont, phasek)
      * settlingModel_->Vk(phaseCont, phasek)
    );
}


// ************************************************************************* //
