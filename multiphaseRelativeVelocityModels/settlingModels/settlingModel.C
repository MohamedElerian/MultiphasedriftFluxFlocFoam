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

#include "settlingModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(settlingModel, 0);
    defineRunTimeSelectionTable(settlingModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModels::settlingModel::settlingModel
(
    const terminalVelocitySettlingModel& UInf,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    UInfModel_(UInf),
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::settlingModels::settlingModel>
Foam::settlingModels::settlingModel::New
(
    const terminalVelocitySettlingModel& UInf,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    const word modelType(dict.lookup(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown time scale model type "
            << modelType << nl << nl
            << "Valid time scale model types :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return
        autoPtr<settlingModel>
        (
            cstrIter()
            (
                UInf,
                dict.optionalSubDict(modelType + "Coeffs"),
                mesh
            )
        );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModels::settlingModel::~settlingModel()
{}


// ************************************************************************* //
