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
    Foam::simple

Description
    General relative velocity model following Rowe

SourceFiles
    simple.C

\*---------------------------------------------------------------------------*/

#ifndef simple_H
#define simple_H

#include "settlingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam::terminalVelocitySettlingModels;

namespace Foam
{
namespace settlingModels
{

/*---------------------------------------------------------------------------*\
                            Class simple Declaration
\*---------------------------------------------------------------------------*/

class simple
:
    public settlingModel
{
    // Private data

        //- Maximum phase fraction
        dimensionedScalar alphatMax_;


public:

    //- Runtime type information
    TypeName("simple");


    // Constructors

        //- Construct from components
        simple
        (
            const terminalVelocitySettlingModel& UInf,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~simple();


    // Member Functions

        //- Return the settling velocity between continuos phase and phase k
        virtual tmp<volScalarField> Vk
        (
            const phaseDrift& phasec,
            const phaseDrift& phasek
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace settlingModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
