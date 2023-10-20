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

#include "multiphaseDriftMixture.H"
#include "Time.H"
#include "subCycle.H"
#include "fvmDdt.H"
#include "CMULES.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"
#include "fvcFlux.H"
#include "fvmLaplacian.H"
#include "fvCFD.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "LESdelta.H"
#include "cyclicFvPatch.H"
#include "fvcLaplacian.H"
#include "laplacianScheme.H"
#include <iostream>
#include <fstream>
using namespace std;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseDriftMixture::calcAlphas()
{
    alphas_ == 0.0;

    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
        alphas_ += iter();
    }
    volScalarField unitalphaCorr 
    (
        IOobject
        (
            "unitalphaCorr ",
             mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar( "unitalphaCorr ", dimless, scalar(1.0))
    );    
    alphasSolids_=  unitalphaCorr *alphas_.weightedAverage(mesh_.V()).value();


}

Foam::wordList Foam::multiphaseDriftMixture::UdmPatchFieldTypes() const
{
    wordList UdmTypes
    (
        U_.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );

    forAll(U_.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
        )
        {
            UdmTypes[i] = fixedValueFvPatchVectorField::typeName;
        }
    }

    return UdmTypes;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseDriftMixture::multiphaseDriftMixture 
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:  //the colon is responsible for all the follwing constructors

    IOdictionary //This is a base class
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(U.mesh()),
    phasesk_(lookup("phasesk"), phaseDrift::iNew(U, phi)),  // construct from Istream
    phaseskName_(phasesk_.toc()),                         
    phasec_(lookup("continousPhase"), phaseDrift::iNew(U, phi)),
    phases_(phasesk_.size() + phasec_.size()),
    UkmPtr_(phases_.size()),                                //
    U_(U),                                                  // U is the mixture velocity
    phi_(phi),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("rhoPhi", dimMass/dimTime, 0.0)
    ),
    alphas_
    (
        IOobject
        (
            "alpha.solids",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphas", dimless, 0.0)
    ),
    kLimit_
    (
        IOobject
        (
            "kLimit",
            mesh_.time().timeName(),
            mesh_,
       		IOobject::NO_READ,
        	IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0, 2, -1, 0, 0, 0), Zero)
    ),
    klim_
    (
        IOobject
        (
            "klim",
            mesh_.time().timeName(),
            mesh_,
       		IOobject::NO_READ,
        	IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0, 2, -1, 0, 0, 0), Zero)
    ),
    alphaDiffusion_1_
    (
        IOobject
        (
            "alphaDiffusion_1",
            mesh_.time().timeName(),
            mesh_,
       		IOobject::NO_READ,
        	IOobject::NO_WRITE

        ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0, 2, -1, 0, 0, 0), Zero)
    ),
    sumAlpham_
    (
        IOobject
        (
            "sumAlpham",
            mesh_.time().timeName(),
            mesh_,
       		IOobject::NO_READ,
        	IOobject::NO_WRITE

        ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0), Zero)
    ),
    alphaBeforeDiffusion_
    (
        IOobject
        (
            "alphaBeforeDiffusion",
            mesh_.time().timeName(),
            mesh_,
       		IOobject::NO_READ,
        	IOobject::NO_WRITE

        ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0), Zero) 
    ),
    alphaAfterDiffusion_
        (
        IOobject
        (
            "alphaAfterDiffusion",
            mesh_.time().timeName(),
            mesh_,
       		IOobject::NO_READ,
        	IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0), Zero)
    ),
    alphaAmountDiffused_
    (
        IOobject
        (
            "alphaAmountrDiffused",
            mesh_.time().timeName(),
            mesh_,
       		IOobject::NO_READ,
        	IOobject::NO_WRITE

        ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0), Zero)
        
    ),  
   kk_
    (
        IOobject
        (
            "kk",
            mesh_.time().timeName(),
            mesh_,
       		IOobject::NO_READ,
        	IOobject::AUTO_WRITE

        ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0, 2, -2, 0, 0, 0), Zero)
        
    ),   
    muModel_
    (
        mixtureViscosityModel::New
        (
            "mu",
            subDict("mixtureViscosityModel").subDict("solids"),
            U,
            phi
        )
    ),
    mu_
    (
        IOobject
        (
            "nu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimTime/dimLength, 0.0)
    ),
    // Maximum sediment volume concentration
    alphaMax_(lookupOrDefault<scalar>("alphaMax", 0.45)),
    C11_(lookupOrDefault<scalar>("C_breakfreq", 1)),
    df_(lookupOrDefault<scalar>("df", 1)),
    kc_(lookupOrDefault<scalar>("kc", 1)),
    packc_(lookupOrDefault<scalar>("C", 1)),
    alphaDiffusion_(lookupOrDefault<scalar>("alphaDiffusion", 1.0)),
    alphaDif_(lookupOrDefault<scalar>("alphaDif", 1.0)),
    relativeUModel_(subDict("relativeVelocityModel"), mesh_),

    alphaPhases_(phasesk_.size()),
    // Collision Efficiency 
    Alpha_(lookupOrDefault<scalar>("Alpha", 1.0)),
    // Emprical break-up parameter
    Eb_(lookupOrDefault<scalar>("Eb", 1e-5)),
    // Floc shear strength
    Fy_(lookupOrDefault<scalar>("Fy", 1e-10)),

    // Collision Frequency Functionw
    Beta_(phasesk_.size()*phasesk_.size()),
    // Collision Frequency Function
    Alphaa_(phasesk_.size()*phasesk_.size()),
    // Break-up frequency
    S_(phasesk_.size()),
    // Flocculation source terms for all sediment phases (thus it is a list)
    SrcList(phasesk_.size()),
    // The summation of flocculation source terms (theoretically be zero)
    SumSource_
    (
        IOobject
        (
            "SumSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("SumSource", dimless/dimTime, scalar(0.0))
    ),
    // Breakup Distribtion Function
    Gamma_(lookupOrDefault<scalar>("Gamma", 2.0)),
    // Turbulent Shear
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimless, scalar(0.0))
    ),
    // Temperature
    T_(lookupOrDefault<scalar>("T", 293.0)),
    // Particle Diameter for all sediment phases
    diameters_(phasesk_.size()),
    // Collision Particle Diameter for all sediment phases
     colldiameters_(phasesk_.size()),

    // Particle Volume for all sediment phases
    volumes_(phasesk_.size()),

    // aggregate radius using the fractal dimension
     aggradius_(phasesk_.size()),

    // aggregate radius using the fractal dimension
     aggporositiy_(phasesk_.size()),
    // aggregate permeability
     aggpermeability_(phasesk_.size()),
    // zeta coefficent
     zetacoef_(phasesk_.size()),
    // J coefffient
     Jcoef_(phasesk_.size()),
    // c coefffient
     ccoef_(phasesk_.size()),
    // c coefffient
     dcoef_(phasesk_.size()),
    // omega coefffient
     omegacoef_(phasesk_.size()),
    // eta coefffient
     etacoef_(phasesk_.size()),



   Source1_
    (
        IOobject
        (
            "Source1_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source1_", dimless, scalar(0.0))
    ),
    Source2_
    (
        IOobject
        (
            "Source2_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source2_", dimless, scalar(0.0))
    ),
    Source3_
    (
        IOobject
        (
            "Source3_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source3_", dimless, scalar(0.0))
    ),
    Source4_
    (
        IOobject
        (
            "Source4_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source4_", dimless, scalar(0.0))
    ),
    Source5_
    (
        IOobject
        (
            "Source5_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source5_", dimless, scalar(0.0))
    ),
    Source6_
    (
        IOobject
        (
            "Source6_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source6_", dimless, scalar(0.0))
    ),
    // Total Sediment Volume Concentration
    alphasSolids_
    (
        IOobject
        (
            "alphasSolids",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasSolids", dimless, scalar(0.0))
    ),
    // Volume Concentration of sediment phase 1,2,3 ...
    alphasPhase1_
    (
        IOobject
        (
            "alphasPhase1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase1", dimless, scalar(0.0))
    ),
    alphasPhase2_
    (
        IOobject
        (
            "alphasPhase2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase2", dimless, scalar(0.0))
    ),
    alphasPhase3_
    (
        IOobject
        (
            "alphasPhase3",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase3", dimless, scalar(0.0))
    ),
    alphasPhase4_
    (
        IOobject
        (
            "alphasPhase4",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase4", dimless, scalar(0.0))
    ),
    alphasPhase5_
    (
        IOobject
        (
            "alphasPhase5",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase5", dimless, scalar(0.0))
    ),
    alphasPhase6_
    (
        IOobject
        (
            "alphasPhase6",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase6", dimless, scalar(0.0))
    ),
    alphasPhase7_
    (
        IOobject
        (
            "alphasPhase7",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase7", dimless, scalar(0.0))
    ),
    alphasPhase8_
    (
        IOobject
        (
            "alphasPhase8",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase8", dimless, scalar(0.0))
    ),
    alphasPhase9_
    (
        IOobject
        (
            "alphasPhase9",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase9", dimless, scalar(0.0))
    ),
    alphasPhase10_
    (
        IOobject
        (
            "alphasPhase10",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase10", dimless, scalar(0.0))
    ),
    alphasPhase11_
    (
        IOobject
        (
            "alphasPhase11",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase11", dimless, scalar(0.0))
    ),
    alphasPhase12_
    (
        IOobject
        (
            "alphasPhase12",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase12", dimless, scalar(0.0))
    ),

    alphasPhase13_
    (
        IOobject
        (
            "alphasPhase13",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase13", dimless, scalar(0.0))
    ),

    alphasPhase14_
    (
        IOobject
        (
            "alphasPhase14",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase14", dimless, scalar(0.0))
    ),

    alphasPhase15_
    (
        IOobject
        (
            "alphasPhase15",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase15", dimless, scalar(0.0))
    ),
    alphasPhase16_
    (
        IOobject
        (
            "alphasPhase16",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase16", dimless, scalar(0.0))
    ),

    alphasPhase17_
    (
        IOobject
        (
            "alphasPhase17",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase17", dimless, scalar(0.0))
    ),
    alphasPhase18_
    (
        IOobject
        (
            "alphasPhase18",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase18", dimless, scalar(0.0))
    ),
    alphasPhase19_
    (
        IOobject
        (
            "alphasPhase19",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase19", dimless, scalar(0.0))
    ),
    alphasPhase20_
    (
        IOobject
        (
            "alphasPhase20",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase20", dimless, scalar(0.0))
    ),
    alphasPhase21_
    (
        IOobject
        (
            "alphasPhase21",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase21", dimless, scalar(0.0))
    ),

alphasPhase22_
    (
        IOobject
        (
            "alphasPhase22",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase22", dimless, scalar(0.0))
    ),

    alphasPhase23_
    (
        IOobject
        (
            "alphasPhase23",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase23", dimless, scalar(0.0))
    ),
    alphasPhase24_
    (
        IOobject
        (
            "alphasPhase24",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase24", dimless, scalar(0.0))
    ),
    alphasPhase25_
    (
        IOobject
        (
            "alphasPhase25",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase25", dimless, scalar(0.0))
    ),
    alphasPhase26_
    (
        IOobject
        (
            "alphasPhase26",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase26", dimless, scalar(0.0))
    ),
    alphasPhase27_
    (
        IOobject
        (
            "alphasPhase27",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphasPhase27", dimless, scalar(0.0))
    ),


    CellHeight_
    (
        IOobject
        (
            "CellHeight",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("CellHeight", dimLength, scalar(0.0))
    ),

    CellVolume_
    (
        IOobject
        (
            "CellVolume",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("CellVolume", dimVolume, scalar(0.0))
    )
//...............start...............................////...............start...............................////...............start...............................////...............start....................
// start to define the multiphaseDriftMixture. 
{
    rhoPhi_.setOriented();

    label phasei = 0;
    // Populate all phases List
   // forAllIters (phasesk_, iter)
   // {   phaseDrift& pm =iter();
       // phases_.set(phasei++, &pm);
   // }
   // forAllIters (phasec_, iter)
   // {
    //    phaseDrift& pm =iter();
   //     phases_.set(phasei++, &pm);
   // }

    forAllConstIters (phasesk_, iter)
    {
        phases_.set(phasei++, iter());
    }
    forAllConstIters (phasec_, iter)
    {
        phases_.set(phasei++, iter());
    }


    forAll (phases_, i)
    {
        UkmPtr_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("Ukm", phases_[i].name()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("Uk", dimVelocity, Zero),
                UdmPatchFieldTypes()
            )
        );
    }

    // Height of Cell Center from the center to the bottom of the domain
    CellHeight_ = mesh_.C().component(vector::Z);
    forAll(CellHeight_.boundaryField(),cellI )
    {
        CellHeight_[cellI] = 0;
    }
    CellVolume_.ref() = mesh_.V();

    // Copy particle diameter
    copyDiameters();

    // Copy collison particle diameter
 	calccollDiameters();

    // Calculate particle volume
    CalcVolumes();
    // calculate the aggregated radius r
    Calcaggradius();
    // calculate the aggregates porositiy
    Calcaggporositiy();
    // calculate the aggregates permeability
    Calcaggpermeability();
    // calculate the zeta coefficnet
    Calczetacoef();
    // calculate J coefficent
    CalcJcoef();
    //calculate c coefficent
    Calcccoef();
    //calculate d coefficent
    Calcdcoef();
    //calculate d coefficent
    Calcomegacoef();
    //calculate d coefficent
    Calcetacoef();
    // Collision Frequency List
    forAll(Beta_,i)
    {
        word nameBeta("Beta"+ i);
        Beta_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    nameBeta,
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(nameBeta, dimless, scalar(0.0))
            )
        );
    }
    // Collision effeciency List
    forAll(Alphaa_,i)
    {
        word nameAlphaa("Alphaa"+ i);
        Alphaa_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    nameAlphaa,
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(nameAlphaa, dimless, scalar(0.0))
            )
        );
    }

    // Flocculation Source Terms
    forAll(SrcList,i)
    {
        word nameSrc("SrcPhase."+phaseskName_[i]);
        SrcList.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    nameSrc,
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(nameSrc, dimless/dimTime, scalar(0.0))
            )
        );
    }

    // Break-up frequency list
    forAll(S_,i)
    {
        word nameS("S"+ i);
        S_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    nameS,
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(nameS, dimless, scalar(0.0))
            )
        );
    }
    // Update mu
    correct();
    // Update alpha.solid
    calcAlphas();
    alphas_.write();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //
// Here definiention of the functions outside the class. 

void Foam::multiphaseDriftMixture::fixedFluxOnPatches  //- fix flux value at patch , boundary
(
    surfaceScalarField& alphaPhiCorr,
    const volScalarField & alpha
)
{
    forAll(alphaPhiCorr.boundaryField(), patchi)
    {
        fvsPatchScalarField& phiAlphaCorrp =
            alphaPhiCorr.boundaryFieldRef()[patchi];

        if (!phiAlphaCorrp.coupled()) // coupled() return false
        {
            const scalarField& phi1p = phi_.boundaryField()[patchi]; // Um
            const scalarField& alpha1p = 
                alpha.boundaryField()[patchi];   // 

            forAll(phiAlphaCorrp, facei)
            {
                if (phi1p[facei] < 0)
                {
                    phiAlphaCorrp[facei] = alpha1p[facei]*phi1p[facei]; 
                }
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////
void Foam::multiphaseDriftMixture::calculateAlphaPhi    //- Calculate alphaPhi for each phase
(
    surfaceScalarField& alphaPhiCorr,
    const phaseDrift& phase,
    label phasei
)
{
    word alphaScheme("div(phi,alpha)");   // read from fvSchemes
    word alpharScheme("div(phirb,alpha)");//

    surfaceScalarField phikm(fvc::flux(UkmPtr_[phasei]));

    alphaPhiCorr =
        fvc::flux
        (
            phi_,
            phase,          // phase contains the volumetric concentration
            alphaScheme
        )
      + fvc::flux
        (
            phikm,
            phase,
            alpharScheme
        );              //  alpha*(Um+Ukm)

    // Ensure that the flux at inflow BCs is preserved
    fixedFluxOnPatches(alphaPhiCorr, phase);
}

////////////////////////////////////////////////////////////////////////////////////////////
void Foam::multiphaseDriftMixture::calculateAlphaPhis  //- The calculateAlphaPhis function is responsible for computing the alphaPhi values for all the phases
(
    PtrList<surfaceScalarField>& alphaPhiCorrs, //A list of surfaceScalarField pointers, which will store the calculated alphaPhi values for each phase.
    const UPtrList<phaseDrift>& phases // A list of phases in the simulation, each represented by a phaseDrift object.
)
{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    int phasei = 0;

/* This calculation is calcuated the second term in the fraction transport equation (\div(alpha_k (Um +Ukm))) where Um is used in calclulaing the Phi_ and Ukm is used to caluclate the phikm_
For each phase, it calculates the alphaPhi value using the alphaScheme and alpharScheme.
alphaScheme computes the alphaPhi using the "div(phi, alpha)" scheme, where phi_ is a flux field and alpha represents the phase's volume fraction.
The computed alphaPhi is stored as a surfaceScalarField with a name that includes the phase name followed by "Corr" (e.g., "phialph3Corr").
alpharScheme is used to calculate a second part of alphaPhi based on the phase's relative velocity.
Both parts are added to alphaPhiCorr.*/

    forAllConstIter(UPtrList<phaseDrift>, phases, iter)    //iteration for all phases: dispersed, continuous
    {
        const phaseDrift& alpha = iter();

        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField   
            (
                "phi" + alpha.name() + "Corr",
                fvc::flux
                (
                    phi_,
                    alpha,
                    alphaScheme
                )
            )
        );  

        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];   // alpha_k* Um

        surfaceScalarField phikm(fvc::flux(UkmPtr_[phasei]));       // alpha_k U_km

        alphaPhiCorr += fvc::flux //alphPhi = alphaPhi (Um) + alphaPhi(Ukm)
        (
            phikm,
            alpha,
            alpharScheme
        );                                 

        // Ensure that the flux at inflow BCs is preserved
        fixedFluxOnPatches(alphaPhiCorr, alpha);

        phasei++;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::multiphaseDriftMixture::massFraction(const phaseDrift& phasek) const
{
    return phasek.rho()*phasek/rho();
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
Foam::tmp<Foam::volVectorField>
Foam::multiphaseDriftMixture::sumCkUkr() const
{
    tmp<volVectorField> tCkUkr
    (
        new volVectorField
        (
            IOobject
            (
                "sumUkm",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("zero", dimVelocity, Zero)
        )
    );
    volVectorField& CkUkr = tCkUkr.ref();

    forAllConstIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
        CkUkr +=
            massFraction(iter())*relativeUModel_.Ukr(phasec_.first(), iter()); // sum(alpha*rho/rhom*Ukr)
    }



    return tCkUkr;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::volSymmTensorField> Foam::multiphaseDriftMixture::tauDm() const
{
    UPtrList<phaseDrift>::const_iterator iter = phases_.begin();

    label phaseId = 0;

    tmp<volSymmTensorField> tauDm(iter()*iter().rho()*sqr(UkmPtr_[phaseId]));

    for (++iter; iter != phases_.end(); ++iter)
    {
        phaseId++;
        tauDm.ref() += iter()*iter().rho()*sqr(UkmPtr_[phaseId]);
    }

    return tauDm;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::volScalarField> Foam::multiphaseDriftMixture::rho() const
{
    UPtrList<phaseDrift>::const_iterator iter = phases_.begin();

    tmp<volScalarField> trho = iter()*iter().rho();
    volScalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter()*iter().rho();
    }

    return trho;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::volScalarField> Foam::multiphaseDriftMixture::rhos() const
{
    PtrDictionary<phaseDrift>::const_iterator iter = phasesk_.begin();

    tmp<volScalarField> trho = iter()*iter().rho();
    volScalarField& rho = trho.ref();

    for (++iter; iter != phasesk_.end(); ++iter)
    {
        rho += iter()*iter().rho();
    }

    return trho;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::volScalarField> Foam::multiphaseDriftMixture::rhoc() const
{
    tmp<volScalarField> trhoc
    (
        new volScalarField
        (
            IOobject
            (
                "trhoc",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            phasec_.first().rho()
        )
    );

    return trhoc;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::scalarField>
Foam::multiphaseDriftMixture::rho(const label patchi) const
{
    UPtrList<phaseDrift>::const_iterator iter = phases_.begin();

    tmp<scalarField> trho = iter().boundaryField()[patchi]*iter().rho().value();
    scalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter().boundaryField()[patchi]*iter().rho().value();
    }

    return trho;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::volScalarField> Foam::multiphaseDriftMixture::mu() const
{
    return mu_;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::scalarField>
Foam::multiphaseDriftMixture::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::volScalarField>
Foam::multiphaseDriftMixture::nu() const
{
    return mu_/rho();
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

Foam::tmp<Foam::scalarField>
Foam::multiphaseDriftMixture::nu(const label patchi) const
{
    return mu_.boundaryField()[patchi]/rho(patchi);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

void Foam::multiphaseDriftMixture::solve()  
{
    //correct();

    const Time& runTime = mesh_.time();

    volScalarField& alpha = phasesk_.first();   //

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

    //

    //correctUkm();

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", rhoPhi_.dimensions(), 0.0)
        );     // initiialize

        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas();
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphas();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

void Foam::multiphaseDriftMixture::correct()
{
    // Correct continuos phase mu
    phasec_.first().correct();      // newtonian fluid: nu-> constant
    mu_ = muModel_->mu(phasec_.first().mu());  // mu=rho*nu  
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Foam::multiphaseDriftMixture::correctUkm()
{
    const volVectorField sumUkr(sumCkUkr());

    // Update Ukm for disperse phases
    for(label i = 0; i < phasesk_.size(); i++)
    {
        UkmPtr_[i] =
            relativeUModel_.Ukr(phasec_.first(), phases_[i]) - sumUkr; // Ukr-sum(Ck*Ukr)  which returns the diffusion velocity

    }
    // Ukm for continuos phase (Ufr = 0)
    UkmPtr_[phasesk_.size()] = -sumUkr; // continuous phase
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

void Foam::multiphaseDriftMixture::solveAlphas()
{
    static label nSolves=-1;
    nSolves++;

    const dictionary& alphaControls = mesh_.solverDict("alpha");  // mesh is an object that is declared and defined In meshobject.H and .C. SolverDict is to call the "alpha from the dictionary"
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr"))); // reading the nAlphaCorr from the dictonary. Which is by the way exists in system/fvSolution
    bool MULESCorr(readBool(alphaControls.lookup("MULESCorr")));// cheking if Mules is activated 
    bool limitAlphaPhi(readBool(alphaControls.lookup("limitAlphaPhi"))); // checking if limitAlphaPhi is activated 

    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size()); // declartion for alphaPhiCorrs as a PtrList that will store a number of surfaceScalarField objects, and it determines the size of the list based on the phases_.size() value.

    calculateAlphaPhis(alphaPhiCorrs, phases_);

    if (limitAlphaPhi)
    {
        label phasei = 0;
        forAllConstIters(phases_, iter)
        {
            MULES::limit
            (
                1.0/mesh_.time().deltaT().value(),
                geometricOneField(),
                iter(),
                phi_,
                alphaPhiCorrs[phasei],
                zeroField(),
                zeroField(),
               // UniformField<scalar>(alphaMax_),
                alphaMax_,
                zeroField(),
                true
            );
            phasei++;
        }

        // Limit the sum off all fluxes
        MULES::limitSum(alphaPhiCorrs);
    }

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("sumAlpha", dimless, 0.0)
    );


    // Reset rhoPhi
    rhoPhi_ = dimensionedScalar("rhoPhi", dimMass/dimTime, 0.0);



    label phasei = 0;

    // Calculate ShearRate
    CalcG();
    // Compute collision Frequency 
    ColFreq();
    // Compute collision effecenci 
    Coleff();
    // Comput Breakage Rate 
    BrkRate();
    // Initialize the summation of flocculation source terms
    forAll(SumSource_ ,cellI)
    {
        SumSource_[cellI] = 0;
    }
    for (label phasei = 0; phasei <= phasesk_.size()-1; phasei++)
    {
        // Calculate Flocculation Source Terms for Each Sediment Phase
        SrcList[phasei]=calSource(phasei);
        SumSource_ += SrcList[phasei]; 
    
    }
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
       // Flocculation source term for phase i
        volScalarField Sa = SrcList[phasei]; 
        phaseDrift& alpha = iter();

        surfaceScalarField alphaPhi
        (
            IOobject
            (
                "alphaPhi",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", phi_.dimensions(), 0.0)
        );

        if (MULESCorr) //Make sure that there is no any fraction value going below  zero (negative)
        {
        
             forAll(alpha,cellI)
            {   
                if(alpha[cellI]<0)
                {
                    alpha[cellI] = 0;
                }

            }
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            fvScalarMatrix alpha1Eqn
            (
                fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
              + fv::gaussConvectionScheme<scalar>
                (
                    mesh_,
                    phi_,
                    upwind<scalar>(mesh_, phi_)
                ).fvmDiv(phi_, alpha)  
            );
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            alpha1Eqn.solve();
            alphaPhi = alpha1Eqn.flux();
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
  }

        for (int aCorr = 0; aCorr<nAlphaCorr; aCorr++)
        {
            if (MULESCorr)
            {
                tmp<surfaceScalarField> talphaPhiCorr
                (
                    alphaPhiCorrs[phasei] - alphaPhi
                );

                volScalarField alpha10("alpha10", alpha);

                MULES::correct //This function will limit eac fraction
                (
                    geometricOneField(),
                    alpha,
		    alphaPhi,
                   // alphaPhiCorrs[phasei],  //F_H 	
                    talphaPhiCorr.ref(), // F_H - F_l
                    zeroField(),
                    zeroField(),
                   alphaMax_,   //Given
                    0
                );

                // Under-relax the correction for all but the 1st corrector
                if (aCorr == 0)
                {
                    alphaPhi += talphaPhiCorr();
                }
                else
                {
                    refCast<volScalarField>(alpha) = 0.5*alpha + 0.5*alpha10;
                    alphaPhi += 0.5*talphaPhiCorr();
                }
            }
            else
            {
                surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];

                MULES::explicitSolve
                (
                    alpha,
                    phi_,
                    alphaPhi,
                    alphaMax_,
                    0
                );
            }

            calculateAlphaPhi(alphaPhiCorrs[phasei], alpha, phasei);
        }

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // Apply the diffusion term separatelly for the limiter
       {
       
        //calculation of the alphaDiffusion 1
        forAll (alphas_, cellI)
        {

             if (alphas_[cellI] > alphaMax_)
                {
                alphaDiffusion_1_[cellI] = alphaDiffusion_;
                }

        }
        //
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha) - fvc::ddt(alpha)
              - fvm::laplacian(turbulencePtr_->nut() + alphaDiffusion_1_, alpha) -Sa 
            );
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            alpha1Eqn.solve(mesh_.solver("alpha1Diffusion"));

            alphaPhi += alpha1Eqn.flux();
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
  }

        rhoPhi_ += alphaPhi*alpha.rho();
        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;
        //

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Unit of concentration
        volScalarField unitalpha 
        (
            IOobject
            (
                "unitalpha",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar( "unitalpha", dimless, scalar(1.0))
        );

        //  Volume Concentration of phase i (saved as VolScalarField)
        if (phasei == 0)
        {   

                alphasPhase1_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        }
        if (phasei == 1)
        {   
        
                alphasPhase2_ = unitalpha*alpha.weightedAverage(mesh_.V()).value();

        }
        if (phasei == 2)
        {   

                alphasPhase3_= unitalpha*alpha.weightedAverage(mesh_.V()).value();
        }
        if (phasei == 3)
        {   

                alphasPhase4_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        }
        if (phasei == 4)
        {   

                alphasPhase5_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        }
        if (phasei == 5)
        {   

                alphasPhase6_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        }    
        if (phasei == 6)
        {   

                alphasPhase7_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        }    
        if (phasei == 7)
        {   

                alphasPhase8_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 8)
        {   

                alphasPhase9_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        }  
        if (phasei == 9)
        {   

                alphasPhase10_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 10)
        {   

                alphasPhase11_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 11)
        {   

                alphasPhase12_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 12)
        {   

                alphasPhase13_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 13)
        {   

                alphasPhase14_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 14)
        {   

                alphasPhase15_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 15)
        {   

                alphasPhase16_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 16)
        {   

                alphasPhase17_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 17)
        {   

                alphasPhase18_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 18)
        {   

                alphasPhase19_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 19)
        {   

                alphasPhase20_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 20)
        {   

                alphasPhase21_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 21)
        {   

                alphasPhase22_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 22)
        {   

                alphasPhase23_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 23)
        {   

                alphasPhase24_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 24)
        {   

                alphasPhase25_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 25)
        {   

                alphasPhase26_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
        if (phasei == 26)
        {   

                alphasPhase27_= unitalpha*alpha.weightedAverage(mesh_.V()).value();

        } 
    

        //
        sumAlpha += alpha;

        phasei++;
    }

sumAlpham_ = alphasPhase1_+alphasPhase2_+alphasPhase3_+alphasPhase4_+alphasPhase5_+alphasPhase6_+alphasPhase7_+alphasPhase8_+alphasPhase9_+alphasPhase10_+alphasPhase11_+alphasPhase12_+alphasPhase13_+alphasPhase14_+alphasPhase15_+alphasPhase16_+alphasPhase17_+alphasPhase18_+alphasPhase19_+alphasPhase20_+alphasPhase21_+alphasPhase22_+alphasPhase23_+alphasPhase24_+alphasPhase25_+alphasPhase26_+alphasPhase27_;    
   
 // ***************************Compute the vlomumetric concentration of continuous phase***************************



    //
     
     volScalarField& alphac= phasec_.first();
     alphac=1.0-sumAlpha;
     phaseDrift& alpha= phasec_.first();
     surfaceScalarField alphaPhi
        (
            IOobject
            (
                "alphaPhi",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", phi_.dimensions(), 0.0)
        );
     calculateAlphaPhi(alphaPhi, alpha, phasei);
     fixedFluxOnPatches(alphaPhi, alpha);
     rhoPhi_ += alphaPhi*alpha.rho();
     Info<< alpha.name() << " volume fraction, min, max = "
         << alpha.weightedAverage(mesh_.V()).value()
         << ' ' << min(alpha).value()
         << ' ' << max(alpha).value()
         << endl;     
     sumAlpha += alpha;

//****************************************************************************************************************
     calcAlphas();
     Info<< "Phase-sum volume fraction, min, max = "
         << sumAlpha.weightedAverage(mesh_.V()).value()
         << ' ' << min(sumAlpha).value()
         << ' ' << max(sumAlpha).value()
         << endl;


}
//****************************************************************************************************************
// Copy Diametre List
void Foam::multiphaseDriftMixture::copyDiameters()
{   
    label m = 0;
    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    diameters_[m] = iter().d().value();
    m++;
    }
}  
//****************************************************************************************************************
// calculate collision Diametre List
void Foam::multiphaseDriftMixture::calccollDiameters()
{   
    label x = 0;
    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    colldiameters_[x] = (diameters_[0]/2) * pow((x+1)/kc_,(1/df_));
    x++;
    }
}  
//****************************************************************************************************************
// Calculate the volume of the particles

void Foam::multiphaseDriftMixture::CalcVolumes()
{
    volumes_[0] = 1.0/6.0*Foam::constant::mathematical::pi*pow(diameters_[0],3.0);
    for (label j =1; j<=volumes_.size()-1;j++)
    {
        volumes_[j] = pow(2.0,j*1.0)*volumes_[0]; // Apply the division rule: v_{i+1} = 2* v_{i}
    }
    
}
//****************************************************************************************************************

// Calculate the radius of aggregates r=r_0 (V/V_0)^(1/2)

void Foam::multiphaseDriftMixture::Calcaggradius()
{
    aggradius_[0] = diameters_[0]/2;

    for (label xx =1; xx<=volumes_.size()-1;xx++)
    {
    aggradius_[xx] = aggradius_[0] * pow((volumes_[xx]/volumes_[0]),1/df_);

    }
    
}
//****************************************************************************************************************
// Calculate the porositiy of aggregates. phipor = 1-C(d_i/d_0)^(d_f-3)

void Foam::multiphaseDriftMixture::Calcaggporositiy()
{
    label xxx = 0;
    

    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    aggporositiy_[xxx] = 1 - (packc_* pow((diameters_[xxx]/diameters_[0]),df_-3));

    xxx++;
    }
    
}
//****************************************************************************************************************
// Calculate the permeability of aggregates. 

void Foam::multiphaseDriftMixture::Calcaggpermeability()
{
    label y = 0;
    

    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    aggpermeability_[y] =  (pow(diameters_[y],2)/72) * (3 + (3/(1-aggporositiy_[y])) - pow((8/(1-aggporositiy_[y]))-3, 0.33333333) ) ;

    y++;
    }
    
}
//****************************************************************************************************************
// Calculate zeta coef. 

void Foam::multiphaseDriftMixture::Calczetacoef()
{
    label yy = 0;
    

    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    zetacoef_[yy] =  aggradius_[yy] / pow(aggpermeability_[yy],0.5) ;


    
yy++;
    }
    
}
//****************************************************************************************************************
// Calculate J coef. 

void Foam::multiphaseDriftMixture::CalcJcoef()
{
   // label yyy = 0;
    

   // forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    for (label yyy =0; yyy<=volumes_.size()-1;yyy++)
    {
    Jcoef_[yyy] =  (2*pow(zetacoef_[yyy],2)) + 3 - (3 *tanh(zetacoef_[yyy]) / zetacoef_[yyy]) ;

    //yyy++;


    }
    
}
//****************************************************************************************************************
// Calculate c coef. 

void Foam::multiphaseDriftMixture::Calcccoef()
{
    label z = 0;
    

    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    ccoef_[z] =  (-1/Jcoef_[z]) * (pow(zetacoef_[z],5) + (6*pow(zetacoef_[z],3)) - ((tanh(zetacoef_[z])/zetacoef_[z])*((3*pow(zetacoef_[z],5))+(6*pow(zetacoef_[z],3))))) ;
   //ccoef_[z] = 1/Jcoef_[z] ;

    z++;
    }
    
}
//****************************************************************************************************************
// Calculate d coef. 

void Foam::multiphaseDriftMixture::Calcdcoef()
{
    label zz = 0;
    

    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    dcoef_[zz] =  (3/Jcoef_[zz]) * (pow(zetacoef_[zz],3)) * (1-(tanh(zetacoef_[zz])/zetacoef_[zz])) ;


    zz++;
    }
    
}

//****************************************************************************************************************
// Calculate omega coef. 

void Foam::multiphaseDriftMixture::Calcomegacoef()
{
    label l = 0;
    

    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    omegacoef_[l] =  (2*pow(zetacoef_[l],2)*(1-(tanh(zetacoef_[l])/zetacoef_[l])))/((2*pow(zetacoef_[l],2)) + (3 *(1-(tanh(zetacoef_[l])/zetacoef_[l])))) ;


    l++;
    }
    
}
//****************************************************************************************************************
// Calculate eta coef. 

void Foam::multiphaseDriftMixture::Calcetacoef()
{
    label nn = 0;
    

    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
    etacoef_[nn] =  1 - (dcoef_[nn]/zetacoef_[nn]) - (ccoef_[nn]/(pow(zetacoef_[nn],2)));



    nn++;
    }
    
}
//****************************************************************************************************************

// Calculate Shear Rate
void Foam::multiphaseDriftMixture::CalcG()
{  
     kk_ = turbulencePtr_->k()();
    volScalarField epsField = turbulencePtr_->epsilon()();
    volScalarField nuField = phasec_.first().nu()();
    volScalarField Gunit 
    (
        IOobject
        (
            "Gunit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Gunit", dimless/dimTime/dimTime, scalar(1.0))
    );
    G_ = pow(epsField/nuField/Gunit,1.0/2.0); //G = (epsilon/nu)^{1/2}

}
//****************************************************************************************************************
// Construct Collision Frequency Matrix
void Foam::multiphaseDriftMixture::ColFreq()
{
    scalarField diameterlist(phasesk_.size());

    label m = 0;
    forAllIter(PtrDictionary<phaseDrift>, phasesk_, iter)
    {
        diameterlist[m] = iter().d().value();
        m++;
    }
    volScalarField Vunit 
    (
        IOobject
        (
            "Vunit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Vunit", dimLength/dimTime, scalar(1.0))
    );
    volScalarField muunit 
    (
        IOobject
        (
            "muunit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("muunit", dimMass/dimLength/dimTime, scalar(1.0))
    );

    scalar Kb = 1.38064852e-23;
     forAll(Beta_,i)
    {
        word nameBeta("Beta"+ i);
       Beta_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    nameBeta,
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(nameBeta, dimless, scalar(0.0))
            )
        );
    }
    
    forAll (diameterlist, i)
    {
        forAll(diameterlist,j)
        {
        label I = j +i*phasesk_.size();
       // volScalarField Beta_shear = (1.0/6.0)*G_*pow((diameterlist[i]+diameterlist[j]),3.0);
        volScalarField Beta_shear = 1.294 * G_*pow((pow(etacoef_[i],0.5)*colldiameters_[i])+(pow(etacoef_[j],0.5)*colldiameters_[j]),3.0);
       //volScalarField Beta_settling = Foam::constant::mathematical::pi/4.0*pow((diameterlist[i]+diameterlist[j]),2.0)*(mag(mag(UkmPtr_[i]-UkmPtr_[j])/Vunit));
        volScalarField Beta_settling = Foam::constant::mathematical::pi *(pow((pow(etacoef_[i],0.5)*aggradius_[i])+(pow(etacoef_[j],0.5)*aggradius_[j]),2.0))*(mag(mag(UkmPtr_[i]-UkmPtr_[j])/Vunit));
        //volScalarField Beta_Brownian = 2.0/3.0*Kb*T_/(phasec_.first().mu()()/muunit)*pow((diameterlist[i]+diameterlist[j]),2.0)/(diameterlist[i]*diameterlist[j]);
        volScalarField Beta_Brownian = (2.0/3.0)*(Kb*T_/(phasec_.first().mu()()/muunit)) * ((1/(omegacoef_[i]*colldiameters_[i]))+((1/(omegacoef_[j]*colldiameters_[j])) )) * (colldiameters_[i]+colldiameters_[j]);
        Beta_[I] = Beta_shear+Beta_settling+Beta_Brownian; 
       // volScalarField Beta_direct = (G_) *(pow((colldiameters_[i]+colldiameters_[j]),3.0))   ;
       // Beta_[I] = C11_* Beta_direct; 

   
	}
    } 

}

//****************************************************************************************************************
// Construct Collision effeciency Matrix
void Foam::multiphaseDriftMixture::Coleff()
{
    scalarField diameterlist(phasesk_.size());

     forAll(Alphaa_,i)
    {
        word nameAlphaa("Alphaa"+ i);
        Alphaa_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    nameAlphaa,
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(nameAlphaa, dimless, scalar(0.0))
            )
        );
    }
    
    forAll (diameterlist, i)
    {
        forAll(diameterlist,j)
        {
        label I = j +i*phasesk_.size();

Alphaa_[I] = (exp(-0.1* pow((1-((i+1)/(j+1))),2) / (pow((i+1)*(j+1),0.1) ))) * 1;
   
	}
    } 

}
//****************************************************************************************************************
// Construct Breakage rate(
void Foam::multiphaseDriftMixture::BrkRate()
{

    volScalarField muunit 
    (
        IOobject
        (
            "muunit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("muunit", dimMass/dimLength/dimTime, scalar(1.0))
    );
    forAll (volumes_,i)
    {
    S_[i]= Eb_*G_*pow(((diameters_[i]-diameters_[0])/diameters_[0]),3.0-df_)*pow((phasec_.first().mu()()/muunit*G_)/(Fy_/pow(diameters_[i],2.0)),0.5);
    } 

}

//****************************************************************************************************************
// 
void Foam::multiphaseDriftMixture::testFun()
{
    ;
}
//****************************************************************************************************************
// // // Calculate source terms
Foam::tmp<Foam::volScalarField>
Foam::multiphaseDriftMixture::calSource(label i)  
{   
    
   //Declare source terms
    volScalarField Source
    (
        IOobject
        (
            "Source",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Source", dimless, scalar(0.0))
    );

    volScalarField Source1
    (
        IOobject
        (
            "Source1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source1", dimless, scalar(0.0))
    );
    volScalarField Source2
    (
        IOobject
        (
            "Source2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source2", dimless, scalar(0.0))
    );
    volScalarField Source3
    (
        IOobject
        (
            "Source3",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source3", dimless, scalar(0.0))
    );
    volScalarField Source4
    (
        IOobject
        (
            "Source4",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source4", dimless, scalar(0.0))
    );
    volScalarField Source5
    (
        IOobject
        (
            "Source5",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source5", dimless, scalar(0.0))
    );
    volScalarField Source6
    (
        IOobject
        (
            "Source6",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Source6", dimless, scalar(0.0))
    );
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //First Term: Unequal Size
    if (i>=2)
    {
        for (label j = 0; j <= i-2 ; j++)
        {   
		Source1+=pow(2.0,(j-i+1)*1.0)*Alpha_*Beta_[j +(i-1)*phasesk_.size()]*phasesk_[phaseskName_[i-1]]/volumes_[i-1]*phasesk_[phaseskName_[j]]/volumes_[j]*volumes_[i];
        }

    }
    else {
	    Source1 == 0;
 	 }
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	-----------------------
    //- Second Term: Equal Size Aggregation
    if (i>=1)
    {  
	Source2 = (0.5)*Alpha_*Beta_[i-1 +(i-1)*phasesk_.size()]*phasesk_[phaseskName_[i-1]]/volumes_[i-1]*phasesk_[phaseskName_[i-1]]/volumes_[i-1] *volumes_[i];
    } 
    else {
	    Source2 == 0;
 	 }
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     //- Third Term: Aggregation with smaller particles
    if(i>=1 && i<=phasesk_.size()-2)
   {
        for (label j = 0; j <= i-1 ; j++)
        
	{      
	Source3+= (-1.0)*pow(2.0,(j-i)*1.0) * Alpha_* Beta_[j+i*phasesk_.size()] * phasesk_[phaseskName_[j]]/volumes_[j] * phasesk_[phaseskName_[i]]/ volumes_[i] * volumes_[i] ;       
        } 
    }
    else {
	    Source3 == 0;
 	 }
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //- Fourth Term: Aggregation with euql or larger particles
  for (label j = i; j <= phasesk_.size()-2; j++)
    {   
         Source4+=(-1.0)*Alpha_*Beta_[j+i*phasesk_.size()]*phasesk_[phaseskName_[j]]/volumes_[j]*phasesk_[phaseskName_[i]]/volumes_[i] *volumes_[i];   
    }

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //- Fifth Term (self-breakage)
   if (i>=1) 
    {
        Source5=(-1.0)*S_[i]*phasesk_[phaseskName_[i]] / volumes_[i] * volumes_[i];
    }
    else {
	    Source5 == 0;
 	 }
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //- Sixth Term (due to breakage of i+1)
     if(i<=phasesk_.size()-2)
    {   
        //- Here set binary breakup by default; Gamma_ = 2, Gamma_/v{i+1}*v{i} = 2 *1/2 = 1
        Source6=S_[i+1]*phasesk_[phaseskName_[i+1]];
    }
      else {
	    Source6 == 0;
 	 }
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Source = Source1+Source2+Source3+Source4+Source5+Source6;

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Source1_ = Source1;
Source2_ = Source2;
Source3_ = Source3;
Source4_ = Source4;
Source5_ = Source5;
Source6_ = Source6;


dimensionedScalar dimCorr("dimCorr",dimless/dimTime,1.0);

return Source*dimCorr;

    
}   

//- Read 
bool Foam::multiphaseDriftMixture::read()
{
    if (transportModel::read())
    {
        // Not possible to chage phases on the run
        return true;
    }
    else
    {
        return false;
    }
}

// Calculate the Source Term


// ************************************************************************* //
