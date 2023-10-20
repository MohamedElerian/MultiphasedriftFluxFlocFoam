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
    Foam::Ferguson

Description
    Ferguson and Church terminal settling velocity model

SourceFiles
    Ferguson.C

\*---------------------------------------------------------------------------*/

#ifndef Ferguson_H
#define Ferguson_H

#include "terminalVelocitySettlingModel.H"
#include "phaseDrift.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace terminalVelocitySettlingModels
{

/*---------------------------------------------------------------------------*\
                         Class Ferguson Declaration
\*---------------------------------------------------------------------------*/

class Ferguson
:
    public terminalVelocitySettlingModel
{

private:

    // Private data


        //- C1 coefficient
        dimensionedScalar C1_;

        //- C1 coefficient
        dimensionedScalar C2_;

        //- Fractal dimension
        dimensionedScalar dff_;

        //- Fractal dimension
        dimensionedScalar d_p_;
       




public:

    //- Runtime type information
    TypeName("Ferguson");


    // Constructors

        //- Construct from components
        Ferguson(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~Ferguson();


    // Member Functions

        //- Return the terminal settling velocity given continuos phase c and
        // phase k
        virtual tmp<volVectorField> V0
        (
            const phaseDrift& phasec,
            const phaseDrift& phasek
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace terminalVelocitySettlingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
