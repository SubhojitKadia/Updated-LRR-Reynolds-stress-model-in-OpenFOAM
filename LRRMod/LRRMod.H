/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::RASModels::LRRMod

Description
    Launder, Reece and Rodi Reynolds-stress turbulence model for
    incompressible and compressible flows.

    Reference:
    \verbatim
        Launder, B. E., Reece, G. J., & Rodi, W. (1975).
        Progress in the development of a Reynolds-stress turbulence closure.
        Journal of fluid mechanics, 68(03), 537-566.
		
		//////////////
		The pressure-strain term and nonlinear boundary functions are modified based on the recommendations in 
        Naot, D., & Rodi, W. (1982). Calculation of secondary currents in channel flow. 
		Journal of the Hydraulics Division, 108(HY8), 948–968.	
		Cokljat, D. P. (1993). Turbulence Models for Non-circular Ducts and Channels. 
		Doctoral Thesis, City University of London.
		Cokljat, D., & Younis, B. A. (1995). Second-Order Closure Study of Open-Channel Flows. 
		Journal of Hydraulic Engineering, 121(2), 94–107.
		//////////////
		
    \endverbatim

    Including the recommended generalized gradient diffusion model of
    Daly and Harlow:
    \verbatim
        Daly, B. J., & Harlow, F. H. (1970).
        Transport equations in turbulence.
        Physics of Fluids (1958-1988), 13(11), 2634-2649.
    \endverbatim

    Optional Gibson-Launder wall-reflection is also provided:
    \verbatim
        Gibson, M. M., & Launder, B. E. (1978).
        Ground effects on pressure fluctuations in the
        atmospheric boundary layer.
        Journal of Fluid Mechanics, 86(03), 491-511.
    \endverbatim

    The model coefficients are:// ref. Cokljat (1993), Cokljat & Younis (1995)
    \verbatim
        LRRModCoeffs
        {
            Cmu             0.09;
            C1              1.5;
            C2              0.4;
            Ceps1           1.45;
            Ceps2           1.9;
            Cs              0.22;
            Ceps            0.18;

            wallReflection  yes;
            kappa           0.41;
            Cref1           0.5;
            Cref2           0.1;

            couplingFactor  1.0;
        }
    \endverbatim

SourceFiles
    LRRMod.C

\*---------------------------------------------------------------------------*/

#ifndef LRRMod_H
#define LRRMod_H

#include "RASModel.H"
#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LRRMod Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class LRRMod
:
    public ReynoldsStress<RASModel<BasicMomentumTransportModel>>
{
protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;

            dimensionedScalar C1_;
            dimensionedScalar C2_;

            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Cs_;
            dimensionedScalar Ceps_;


        // Wall-refection coefficients

            Switch wallReflection_;
            dimensionedScalar kappa_;
            dimensionedScalar Cref1_;
            dimensionedScalar Cref2_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
			volVectorField CC_;//added


    // Protected Member Functions

        //- Update the eddy-viscosity
        virtual void correctNut();
		
		//- Source term for the epsilon equation
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("LRRMod");


    // Constructors

        //- Construct from components
        LRRMod
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        LRRMod(const LRRMod&) = delete;


    //- Destructor
    virtual ~LRRMod()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }
		
		//- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New
            (
                "omega",
                epsilon_/(Cmu_*k_),
                epsilon_.boundaryField().types()
            );
        }
		
		//- Return the cell centes//added
        virtual tmp<volVectorField> cellCentres() const
        {
            return CC_;
        }

        //- Return the effective diffusivity for R
        tmp<volSymmTensorField> DREff() const;
		//tmp<volScalarField> DREff() const;

        //- Return the effective diffusivity for epsilon
        tmp<volSymmTensorField> DepsilonEff() const;
		//tmp<volScalarField> DepsilonEff() const;

        //- Solve the turbulence equations and correct eddy-Viscosity and
        //  related properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LRRMod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LRRMod.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
