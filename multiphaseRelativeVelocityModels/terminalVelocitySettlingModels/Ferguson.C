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

#include "Ferguson.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace terminalVelocitySettlingModels
{
    defineTypeNameAndDebug(Ferguson, 0);
    addToRunTimeSelectionTable
    (
        terminalVelocitySettlingModel,
        Ferguson,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::terminalVelocitySettlingModels::Ferguson::Ferguson
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    terminalVelocitySettlingModel(dict, mesh),
    C1_("C1", dimless, dict),
    C2_("C2", dimless, dict),
    dff_("df2", dimless, dict),
    d_p_("dp", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::terminalVelocitySettlingModels::Ferguson::~Ferguson()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::terminalVelocitySettlingModels::Ferguson::V0
(
    const phaseDrift& phasec,
    const phaseDrift& phasek
) const
{
    tmp<volVectorField> tV0
    (
        new volVectorField
        (
            IOobject
            (
                "V0",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("zero", dimVelocity, vector::zero)
        )
    );
    volVectorField& V0 = tV0.ref();
    
    //const uniformDimensionedVectorField& gField(mesh_.lookupObject<uniformDimensionedVectorField>("g"));
     
    const  uniformDimensionedVectorField gField
    (
        IOobject
        (
        "g",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        )
    );  

    const vector g(gField.value());

    const volScalarField& nuf = phasec.nu()();

    const dimensionedScalar R((phasek.rho() - phasec.rho())/phasec.rho());


    



//This calculation is not based on Frguson. It is based on the paper of "A simple method for calculating in situ floc settling velocities based on effective density functions" of the flocculation. 

    forAll (V0, i)  
    {
        V0[i] =
            R.value()*g*pow(phasek.d().value(), dff_.value()-1)
         / (
             C1_.value()*nuf[i] * pow(d_p_.value(), dff_.value()-3 )
           + (
                C2_.value()*Foam::sqrt
             (
                R.value()*mag(g)*pow(phasek.d().value(),dff_.value())* pow(d_p_.value(), dff_.value()-3 )
             )
	     )
           );
    /*   V0[i] =
            R.value()*g*sqr(phasek.d().value())
         / (
             C1_.value()*nuf[i]
           + Foam::sqrt
             (
                 0.75*C2_.value()*R.value()*mag(g)*pow3(phasek.d().value())
             )
           );*/
    }

    return tV0;
}


// ************************************************************************* //
