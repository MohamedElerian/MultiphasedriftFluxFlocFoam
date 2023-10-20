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
    Foam::phaseDrift

Description
    Single incompressible phaseDrift derived from volScalarField and
    containing information of the phaseDrift (rho, d, nuModel and U).

    Used as part of the multiPhaseDriftMixture to hold phaseDrift data

SourceFiles
    phaseDrift.C

\*---------------------------------------------------------------------------*/

#ifndef phaseDrift_H
#define phaseDrift_H

#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseDrift Declaration
\*---------------------------------------------------------------------------*/

//The phaseDrift class is defined, which is derived from volScalarField. This indicates that phaseDrift is a type of scalar field used for modeling phase properties
class phaseDrift
:
    public volScalarField  
{

    // Private data


         //- reference to mixture phi
        const surfaceScalarField& phi_;

        //- reference to mixture U
        const volVectorField& U_;

        //- Name of the phaseDrift
        word name_;

        //- Phase dictionary
        dictionary dict_;

        //- Phase density
        dimensionedScalar rho_;

        //- contiuos Phase density
        //dimensionedScalar rhocc_;



        //- Phase diameter for discrete phases
        dimensionedScalar d_;

        //- Viscosity model for continous phase
        autoPtr<viscosityModel> nuModel_;


public:

    // Constructors

        //- Construct from components
        phaseDrift
        (
            const surfaceScalarField& phi,
            const volVectorField& U,
            const word& name,
            const dictionary& dict
        );

        //- Return clone
        autoPtr<phaseDrift> clone() const;

        //- Return a pointer to a new phaseDrift created on freestore
        //  from Istream
        class iNew
        {

            const fvMesh& mesh_;
            const volVectorField& U_;
            const surfaceScalarField& phi_;

        public:

            iNew(const volVectorField& U, const surfaceScalarField& phi)
            :
                mesh_(U.mesh()),
                U_(U),
                phi_(phi)
            {}

            autoPtr<phaseDrift> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<phaseDrift>
                (
                    new phaseDrift(phi_, U_, ent.keyword(), ent)
                );
            }
        };


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
        }

        //- Return const-access to viscosityModel
        const viscosityModel& nuModel() const
        {
            return nuModel_();
        }

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nuModel_->nu();
        }

        //- Return the kinematic laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nuModel_->nu(patchi);
        }

        //- Return the dynamic viscosity
        tmp<volScalarField> mu() const
        {
            return nuModel_->nu()*rho_;
        }

        //- Return the dynamic laminar viscosity for patch
        tmp<scalarField> mu(const label patchi) const
        {
            return nuModel_->nu(patchi)*rho_.value();
        }

        //- Return const-access to phaseDrift1 density
        const dimensionedScalar& rho() const
        {
            return rho_;
        }




        //- Return const-access to phaseDrift diameter
        const dimensionedScalar& d() const
        {
            return d_;
        }

        //- Correct the phaseDrift properties
        void correct();

        //- Read base transportProperties dictionary
        bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
