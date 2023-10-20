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

#include "phaseDrift.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Constructor for the phaseDrift class
// Initializes a phaseDrift object with provided parameters and configuration settings.

Foam::phaseDrift::phaseDrift
(
    const surfaceScalarField& phi,  // Reference to a surfaceScalarField object for phase information (rho*U)
    const volVectorField& U,        // Reference to a volVectorField object for velocity information
    const word& phaseDriftName,     // Name of the phaseDrift
    const dictionary& dict         // Dictionary containing the simulation configuration settings (transportproprteis)
)
:
    volScalarField
    (
        // Initialize the volScalarField base class with I/O settings
        IOobject
        (
            IOobject::groupName("alpha", phaseDriftName),   // Field name and group name
            U.mesh().time().timeName(),                     // Current time for I/O
            U.mesh(),                                       // Mesh associated with the field
            IOobject::MUST_READ,                            // Must read from input
            IOobject::AUTO_WRITE                           // Automatically write to output
        ),
        U.mesh()                                           // Mesh associated with the field
    ),
    phi_(phi),               // Initialize the phi_ member variable with the provided surfaceScalarField
    U_(U),                   // Initialize the U_ member variable with the provided volVectorField
    name_(phaseDriftName),   // Initialize the name_ member variable with the provided phaseDriftName
    dict_(dict),             // Initialize the dict_ member variable with the provided dictionary
    rho_("rho", dimDensity, dict_),   // Initialize the rho_ member variable with the keyword "rho," dimension of density, and dictionary
    d_("d", dimLength, -1)           // Initialize the d_ member variable with the keyword "d," dimension of length, and default value -1
{
    // Check if the "transportModel" keyword is found in the dictionary
    if (dict_.found("transportModel"))
    {
        // If found, create a viscosity model and associate it with "nu" field
        nuModel_.reset
        (
            viscosityModel::New
            (
                IOobject::groupName("nu", phaseDriftName),  // Group name for nu field
                dict_,                                      // Dictionary with configuration settings
                U_,                                         // Velocity field U
                phi_                                        // Phase field phi
            ).ptr()
        );
    }

    // Set the value of the d_ member variable to the value found in the dictionary with the key "d"
    // If not found, set it to the default value of -1 with the dimension of length (dimLength).
    d_ = d_.lookupOrDefault("d", dict_, dimLength);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseDrift> Foam::phaseDrift::clone() const
{
    NotImplemented;
    return autoPtr<phaseDrift>(nullptr);
}

// Correct the phaseDrift properties
// If a viscosity model (nuModel_) is available, call its correct() function to update properties.

void Foam::phaseDrift::correct()
{
    if (nuModel_.valid())
    {
        nuModel_->correct();
    }
}

// Read phaseDrift properties from a dictionary (transportProprties)
// Reads properties such as density (rho) and diameter (d) from the provided phaseDriftDict dictionary.
// Returns true if the properties were successfully read, otherwise returns false.

bool Foam::phaseDrift::read(const dictionary& phaseDriftDict)
{
    bool readOK = true;

    if (readOK)
    {
        dict_.lookup("rho") >> rho_;
        dict_.lookup("d") >> d_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
