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
    Foam::multiphaseDriftMixture

Description
    Incompressible multi-phase mixture with built in solution for the
    phase fractions



SourceFiles
    multiphaseDriftMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseDriftMixture_H
#define multiphaseDriftMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "phaseDrift.H"
#include "PtrDictionary.H"
#include "UPtrList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mixtureViscosityModel.H"
#include "multiphaseRelativeVelocityModel.H"
#include "CompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiphaseDriftMixture Declaration
\*---------------------------------------------------------------------------*/

class multiphaseDriftMixture
:
    public IOdictionary,
    public transportModel
{

private:


    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Ptr for solid phases
        PtrDictionary<phaseDrift> phasesk_; 

        //- Name of phases
        wordList phaseskName_;

        //- Ptr for continuos phase
        PtrDictionary<phaseDrift> phasec_;

        //- UPtrList for all phases
        UPtrList<phaseDrift> phases_;

        //- PtrList for all phase relative velocities with the mixture (Ukm)
        PtrList<volVectorField> UkmPtr_;

        //- Reference to mixture U
        const volVectorField& U_;

        //- Reference to volumetric mixture flux
        const surfaceScalarField& phi_;

        //- Reference to mass mixture flux
        surfaceScalarField rhoPhi_;

        //- Total solid phase fraction
        volScalarField alphas_;

        //- Total Solid Phase Limiter
        volScalarField kLimit_;

        //- Total Soild Phase Limiter
        volScalarField klim_;

        //- Total Solid Phase Diffusion Coefficient
        volScalarField alphaDiffusion_1_;

        //- Total Solid Phase to Ensure it is 1
        volScalarField sumAlpham_;

        //- alphas_ before Diffusion
        volScalarField alphaBeforeDiffusion_;

        //- alphas_ after Diffsusion
        volScalarField alphaAfterDiffusion_;

        //- alphas_ amount diffused
        volScalarField alphaAmountDiffused_;

        //- turbulence kinetick energy
        volScalarField kk_;

        //- Reference to the turbulence model
        compressibleTurbulenceModel* turbulencePtr_;

        //- Viscosity model of the mixture
        autoPtr<mixtureViscosityModel> muModel_;

        //- Viscosity of the mixture
        volScalarField mu_;

        //- Maximum dispersed phase-fraction
        scalar alphaMax_;

        //- Maximum dispersed phase-fraction
        scalar C11_;

	//Mass frac- tal dimension ofthe flocs, which is an indicator of the space occupied by the particles in a given aggregate and its value is between 1 (a line of particles) and 3 (a solid sphere).
	scalar df_;

	//Proportionality constant related to the particle packing density, which is generally taken as 1 
	scalar kc_;

	//Packing coefficent, which is generally taken as 1 
	scalar packc_;

        //- Numericial Diffision 
        scalar alphaDiffusion_;

        //- Diffusion Unit
        dimensionedScalar alphaDif_;

        //- Model for relative velocities for phases
        multiphaseRelativeVelocityModel relativeUModel_;

        //- Total volume concentration of sediment
	    scalarField alphaPhases_;

        //- Collision Efficiency
        scalar Alpha_;

        //- A constant parameter for break-up process
        scalar Eb_;

        //- Floc Strength
        scalar Fy_;


        // Collision Frequency
        PtrList<volScalarField>  Beta_;

        // Collision Frequency
        PtrList<volScalarField>  Alphaa_;

        // Breakage Rate
        PtrList<volScalarField>  S_;

        // A list to save phase transition terms, which is volScalarField for each dispersed phase
        PtrList<volScalarField> SrcList;

        // The summation of all phase transition terms, which should be zero theoretically
        volScalarField SumSource_;

        // Breakage Distribution Function
        scalar Gamma_;

        // Shear Rate
        volScalarField G_;
        // Temperature
        scalar T_;

        // DiameterList
        scalarField diameters_;

        // collision diameter list R_ci
        scalarField colldiameters_;

        // Particle volume
        scalarField volumes_;

        // Radius of the aggregates taking into account the fractal shape
        scalarField aggradius_;

        // porositiy of the aggregates
        scalarField aggporositiy_;

        // permeability of the aggregates
        scalarField aggpermeability_;

        // Zeta coeffient
        scalarField zetacoef_;

        // J coeffient
        scalarField Jcoef_;

        // c coeffient
        scalarField ccoef_;

        // d coeffient
        scalarField dcoef_;

        // omega coeffient:the ratio of the force exerted by the fluid on a permeable aggregate to the force exerted by the fluid on an impervious sphere of equivalent size
        scalarField omegacoef_;
        // eta coeffient: the ratio of the flow moving through an aggregate to the total flow approaching the aggregate
        scalarField etacoef_;

        // For Test
        volScalarField Source1_;
        volScalarField Source2_;
        volScalarField Source3_;
        volScalarField Source4_;
        volScalarField Source5_;
        volScalarField Source6_;
        //
        volScalarField alphasSolids_;
        volScalarField alphasPhase1_;
        volScalarField alphasPhase2_;
        volScalarField alphasPhase3_;
        volScalarField alphasPhase4_;
        volScalarField alphasPhase5_;
        volScalarField alphasPhase6_;
        volScalarField alphasPhase7_;
        volScalarField alphasPhase8_;
        volScalarField alphasPhase9_;
        volScalarField alphasPhase10_;
        volScalarField alphasPhase11_;
        volScalarField alphasPhase12_;
        volScalarField alphasPhase13_;
        volScalarField alphasPhase14_;
        volScalarField alphasPhase15_;
        volScalarField alphasPhase16_;
        volScalarField alphasPhase17_;
        volScalarField alphasPhase18_;
        volScalarField alphasPhase19_;
        volScalarField alphasPhase20_;
        volScalarField alphasPhase21_;
        volScalarField alphasPhase22_;
        volScalarField alphasPhase23_;
        volScalarField alphasPhase24_;
        volScalarField alphasPhase25_;
        volScalarField alphasPhase26_;
        volScalarField alphasPhase27_;
	//- Height of cell centers
	volScalarField CellHeight_;
	volScalarField CellVolume_;


    // Private member functions

        //- Calculate alphas
        void calcAlphas();

        //- Return the list of patchFieldTypes for Udm derived from U
        wordList UdmPatchFieldTypes() const;

        //- Calculate alphaPhi for all phases
        void calculateAlphaPhis
        (
            PtrList<surfaceScalarField>& alphaPhiCorrs,
            const UPtrList<phaseDrift>& phases
        );

        //- Calculate alphaPhi for each phase
        void calculateAlphaPhi
        (
            surfaceScalarField& alphaPhiCorr,
            const phaseDrift& alpha,
            label phasei
        );

        //- fix flux value at patch
        void fixedFluxOnPatches
        (
            surfaceScalarField& alphaPhiCorr,
            const volScalarField & alpha
        );

        //- Solve phase transport Equations
        void solveAlphas();

        //- Calculate source terms
            tmp<volScalarField> calSource(label i);

         // Inter phase functions

            //- Return mass concentration of phase k
            tmp<volScalarField> massFraction(const phaseDrift& phasek) const;

            //- Return the sum of relative velocities between phase k and the
            // mixture
            tmp<volVectorField> sumCkUkr() const;


public:

    // Constructors

        //- Construct from components
        multiphaseDriftMixture  //this is a consturctor of the multiphaseDriftMixture Class 
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~multiphaseDriftMixture()
    {}


    // Member Functions

        // Access functions

            //- Return the phases
            const PtrDictionary<phaseDrift>& phasesk() const
            {
                return phasesk_;
            }

            const PtrDictionary<phaseDrift>& phasec() const
            {
                return phasec_;
            }

            //- Return the mixture velocity
            const volVectorField& U() const
            {
                return U_;
            }

            //- Return the volumetric flux
            const surfaceScalarField& phi() const
            {
                return phi_;
            }

            const surfaceScalarField& rhoPhi() const
            {
                return rhoPhi_;
            }

             //- Set turbulence
            void setTurbulence(compressibleTurbulenceModel& turbModel)
            {
                turbulencePtr_ = &turbModel;
            }

        // Inter phase functions

            //- Return the stress tensor due to the phase transport
            tmp<volSymmTensorField> tauDm() const;



        // Transport

            //- Return the mixture density
            tmp<volScalarField> rho() const;

            //- Return the solid density
            tmp<volScalarField> rhos() const;

            //- Return the  continous density
            tmp<volScalarField> rhoc() const;

            //- Return the mixture density for patch
            tmp<scalarField> rho(const label patchi) const;

            //- Return the dynamic laminar viscosity
            virtual tmp<volScalarField> mu() const;

            //- Return the dynamic laminar viscosity for patch
            virtual tmp<scalarField> mu(const label patchi) const;

            //- Return the kinematic laminar viscosity
            virtual tmp<volScalarField> nu() const;

            //- Return the laminar viscosity for patch
            virtual tmp<scalarField> nu(const label patchi) const;

            //- Correct the mixture mu
            virtual void correct();

            //- Correct phase velocities Ukm
            void correctUkm();

            //- Solve for the mixture phase-fractions
            void solve();

            //- Diameter List 
            void copyDiameters();
            //- collision Diameter List R_ci
            void calccollDiameters();


            //- Calculate Volumes
            void CalcVolumes();

            //- Radius of the aggreagates List 
            void Calcaggradius();

            //- Radius of the aggreagates List 
            void Calcaggporositiy();
            //- permeability of the aggreagates List
	    void Calcaggpermeability();	

            //- zeta coefficent
	    void Calczetacoef();		
	    //- J coefficent
	    void CalcJcoef();
	    //- c coefficent
	    void Calcccoef();
    	    //- d coefficent
	    void Calcdcoef();
    	    //- omega coefficent
	    void Calcomegacoef();
    	    //- eta coefficent
	    void Calcetacoef();

            //- Calculate Shear Rate
            void CalcG();

            //- Construct Collision Frequency
            void ColFreq();
            //- Construct Collision Frequency
            void Coleff();

            //- Breakage Rate 
            void BrkRate();
            //
            void testFun();
        // I-O
            //- wirte phase
            void writePhase();

            //- Read base transportProperties dictionary
            bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
