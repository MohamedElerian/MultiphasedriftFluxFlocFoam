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
    Foam::terminalVelocitySettlingModel

Description


SourceFiles
    terminalVelocitySettlingModel.C

\*---------------------------------------------------------------------------*/

#ifndef terminalVelocitySettlingModel_H
#define terminalVelocitySettlingModel_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "phaseDrift.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace terminalVelocitySettlingModels
{
/*---------------------------------------------------------------------------*\
                Class terminalVelocitySettlingModel Declaration
\*---------------------------------------------------------------------------*/

class terminalVelocitySettlingModel
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        terminalVelocitySettlingModel(const terminalVelocitySettlingModel&);

        //- Disallow default bitwise assignment
        void operator=(const terminalVelocitySettlingModel&);


protected:

    // Protected data

        //- Model dictinary
        const dictionary dict_;

        //- fvMesh reference
        const fvMesh& mesh_;


public:


    //- Runtime type information
    TypeName("terminalVelocitySettlingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        terminalVelocitySettlingModel,
        dictionary,
        (const dictionary& dict, const fvMesh& mesh),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        terminalVelocitySettlingModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Selector
    static autoPtr<terminalVelocitySettlingModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~terminalVelocitySettlingModel();


    // Member Functions


        //- Return the terminal settling velocity given continuos phase and
        // phase k
        virtual tmp<volVectorField> V0
        (
            const phaseDrift& phasec,
            const phaseDrift& phasek
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace terminalVelocitySettlingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
