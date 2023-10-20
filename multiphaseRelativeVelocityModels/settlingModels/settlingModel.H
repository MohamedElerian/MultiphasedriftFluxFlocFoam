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

Class
    Foam::settlingModel

Description
    Base class for hindered settling functions

SourceFiles
    settlingModel.C

\*---------------------------------------------------------------------------*/

#ifndef settlingModel_H
#define settlingModel_H

#include "dictionary.H"
#include "phaseDrift.H"
#include "runTimeSelectionTables.H"
#include "terminalVelocitySettlingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam::terminalVelocitySettlingModels;

namespace Foam
{
namespace settlingModels
{

/*---------------------------------------------------------------------------*\
                          Class settlingModel Declaration
\*---------------------------------------------------------------------------*/

class settlingModel
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        settlingModel(const settlingModel&);

        //- Disallow default bitwise assignment
        void operator=(const settlingModel&);


    // Private data

        //- Terminal settling velocity model
        const terminalVelocitySettlingModel& UInfModel_;

        //- Model dictinary
        const dictionary dict_;


protected:

    // Protected data

        //- fvMesh reference
        const fvMesh& mesh_;


public:


    //- Runtime type information
    TypeName("settlingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        settlingModel,
        dictionary,
        (
            const terminalVelocitySettlingModel& UInf,
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (UInf, dict, mesh)
    );


    // Constructors

        //- Construct from components
        settlingModel
        (
            const terminalVelocitySettlingModel& UInf,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Selector
    static autoPtr<settlingModel> New
    (
        const terminalVelocitySettlingModel& UInf,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~settlingModel();


    // Member Functions

        //- Return UInfModel model
        const terminalVelocitySettlingModel& UInfModel() const
        {
            return UInfModel_;
        }

        //- Return the settling velocity between continuos phase c and phase k
        virtual tmp<volScalarField> Vk
        (
            const phaseDrift& phasec,
            const phaseDrift& phasek
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace settlingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
