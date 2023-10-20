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
    Foam::multiphaseRelativeVelocityModel

Description
    This is a based class that described a model for the terminal settling
    velocity following the relationship described by Richarson and Zaki.
    /verbatin
    J.F. Richarson, W.F. Zaki, Chem Eng Sci 1954, 3, 65
    /endverbatin

SourceFiles
    multiphaseRelativeVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseRelativeVelocityModel_H
#define multiphaseRelativeVelocityModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "phaseDrift.H"

#include "terminalVelocitySettlingModel.H"
#include "settlingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam::terminalVelocitySettlingModels;
using namespace Foam::settlingModels;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multiphaseRelativeVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class multiphaseRelativeVelocityModel
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiphaseRelativeVelocityModel(const multiphaseRelativeVelocityModel&);

        //- Disallow default bitwise assignment
        void operator=(const multiphaseRelativeVelocityModel&);


    // Private data

        //- Terminal settling velocity model
        autoPtr<terminalVelocitySettlingModel> UInfModel_;

        //- Hindered settling function model
        autoPtr<settlingModel> settlingModel_;


public:


    // Constructors

        //- Construct from components
        multiphaseRelativeVelocityModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~multiphaseRelativeVelocityModel();


    // Member Functions

        //- Return the relative velocity between phaseCont and phasek
        tmp<volVectorField> Ukr
        (
            const phaseDrift& phaseCont,
            const phaseDrift& phasek
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
