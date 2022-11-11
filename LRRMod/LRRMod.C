/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "LRRMod.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "wallDist.H"
//added
#include "bound.H"
#include "vectorIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void LRRMod<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = this->Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> LRRMod<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
LRRMod<BasicMomentumTransportModel>::LRRMod
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    ReynoldsStress<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
// existing coefficient were as per Launder et al. (1975), Gibson and Launder (1978), Launder (1992) (Wilcox 2006)
// revised coefficient are as per Cokljat (1993), Cokljat and Younis (1995)
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            //1.8
			1.5
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            //0.6
			0.4
        )
    ),
    Ceps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            this->coeffDict_,
            //1.44
			1.45
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            //1.92
			1.9
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            //0.25
			0.22
        )
    ),
    Ceps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps",
            this->coeffDict_,
            //0.15
			0.18
        )
    ),

    wallReflection_
    (
        Switch::lookupOrAddToDict
        (
            "wallReflection",
            this->coeffDict_,
            true
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    Cref1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cref1",
            this->coeffDict_,
            0.5
        )
    ),
    Cref2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cref2",
            this->coeffDict_,
            //0.3
			0.1
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
	CC_//added
    (
        IOobject
        (
            IOobject::groupName("cellCentres", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.0*this->mesh_.C()
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);

        this->boundNormalStress(this->R_);
        bound(epsilon_, this->epsilonMin_);
        k_ = 0.5*tr(this->R_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool LRRMod<BasicMomentumTransportModel>::read()
{
    if (ReynoldsStress<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        Ceps_.readIfPresent(this->coeffDict());

        wallReflection_.readIfPresent("wallReflection", this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        Cref1_.readIfPresent(this->coeffDict());
        Cref2_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> LRRMod<BasicMomentumTransportModel>::DREff() const
{
    return volSymmTensorField::New
    (
        "DREff",
        (Cs_*(this->k_/this->epsilon_))*this->R_ + I*this->nu()
    );
}


template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> LRRMod<BasicMomentumTransportModel>::DepsilonEff() const
{
    return volSymmTensorField::New
    (
        "DepsilonEff",
        (Ceps_*(this->k_/this->epsilon_))*this->R_ + I*this->nu()
    );
}


template<class BasicMomentumTransportModel>
void LRRMod<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    ReynoldsStress<RASModel<BasicMomentumTransportModel>>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));// for Pij
    const volTensorField& gradU = tgradU();// for Pij// read only
	
	tmp<volTensorField> ttgradU(fvc::grad(U)().T());//added for Dij, transpose of grad(U)
    const volTensorField& gradU1 = ttgradU();//added for Dij

    volSymmTensorField P(-twoSymm(R & gradU));// Production tensor Pij = -(Rim.dUj/dxm + Rjm.dUi/dxm)
	volSymmTensorField D(-twoSymm(R & gradU1));//added, Dij = -(Rim.dUm/dxj + Rjm.dUm/dxi)
	volSymmTensorField S(symm(gradU));//added, Sij = 1/2 (gradU + T(gradU))
    volScalarField G(this->GName(), 0.5*mag(tr(P)));// Production of TKE

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1_*alpha*rho*G*epsilon_/k_
      - fvm::Sp(Ceps2_*alpha*rho*epsilon_/k_, epsilon_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );
	
	epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Correct the trace of the tensorial production to be consistent
    // with the near-wall generation from the wall-functions
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label celli = curPatch.faceCells()[facei];
				P[celli] *= min
                (
                    G[celli]/(0.5*mag(tr(P[celli])) + small),
                    1.0
                );
            }
        }
    }
	
	// Reynolds stress equation//modified
    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(alpha, rho, R)
      + fvm::div(alphaRhoPhi, R)
      - fvm::laplacian(alpha*rho*DREff(), R)
      + fvm::Sp(C1_*alpha*rho*epsilon_/k_, R)
      ==
		alpha*rho*P
      - (2.0/3.0*(1 - C1_)*I)*alpha*rho*epsilon_
      - (((C2_ + 8.0)/11.0)*alpha*rho)*dev(P)
	  - (((8.0*C2_ - 2.0)/11.0)*alpha*rho)*(D - P + dev(P))
	  - (((30.0*C2_ - 2.0)/55.0)*alpha*rho*k_*2.0)*S
      + this->RSource()
      + fvModels.source(alpha, rho, R)
    );
    // Optionally add wall-refection term
	// reflect term is modified for the rapid pressure-strain part
	// "topSurface" boundary as symmetryPlane for free surface flows or as wall for duct flows
	// applicable for rectangular or square channels and ducts
    if (wallReflection_)
    {
        const volVectorField& n_(wallDist::New(this->mesh_).n());// modified later
        const volScalarField& y_(wallDist::New(this->mesh_).y());// modified later
		//// added
		// inputs
		const dimensionedScalar h = max(this->mesh_.Cf().component(1));//flow depth in m
		const dimensionedScalar b = 2.0*max(mag(this->mesh_.Cf().component(2))); // channel width in m. half channel width is simulated along the -z direction
		volScalarField yc("yc", y_);
		volScalarField zc("zc", y_);
		// initialization of the average distance from the boundary
		volScalarField y1("y1", y_);
		volScalarField y2("y2", y_);
		volScalarField z1("z1", y_);
		volScalarField z2("z2", y_);
		// initialization of the unit normal vector 
		volVectorField n1("n1", n_);
		volVectorField n2("n2", n_);
		volVectorField n3("n3", n_);
		volVectorField n4("n4", n_);
		
		// to get mean Rxy/k and Rxz/k near the walls
		const dimensionedScalar cell_Y = min(CC_().component(1));// center of 1st cell from bed
		const dimensionedScalar cell_Z = min(CC_().component(2));// center of 1st cell from left wall
		scalar sum1 = 0.0;
		scalar count1 = 0.0;
		scalar sum2 = 0.0;
		scalar count2 = 0.0;
		
		forAll(n_, ii)
		{
			n1[ii] = vector(0, 1, 0);// for top boundary
			n2[ii] = vector(0, -1, 0);// for bottom wall
			n3[ii] = vector(0, 0, -1);// for left wall
			n4[ii] = vector(0, 0, 1);// for right wall
			yc[ii] = CC_[ii].component(1); // y in the code = z in the manuscript
			zc[ii] = CC_[ii].component(2); // z in the code = y in the manuscript
			scalar y_1 = h.value() - yc[ii]; // cell center from the top boundary
			scalar z_1 = 0.5*b.value() - mag(zc[ii]); // cell center from the left wall
			scalar z_2 = b.value() - z_1; // cell center from the right wall
			y1[ii] = sqrt(3.1416*sqr(y_1)/(atan(z_2/y_1) + atan(z_1/y_1) 
					+ 1.0/(y_1/z_2 + z_2/y_1) + 1.0/(y_1/z_1 + z_1/y_1)));// for top boundary
			y2[ii] = sqrt(3.1416*sqr(yc[ii])/(atan(z_2/yc[ii]) + atan(z_1/yc[ii]) 
					+ 1.0/(yc[ii]/z_2 + z_2/yc[ii]) + 1.0/(yc[ii]/z_1 + z_1/yc[ii])));// for bottom wall
			z1[ii] = sqrt(3.1416*sqr(z_1)/(atan(yc[ii]/z_1) + atan(y_1/z_1) 
					+ 1.0/(z_1/yc[ii] + yc[ii]/z_1) + 1.0/(z_1/y_1 + y_1/z_1)));// for left wall
			z2[ii] = sqrt(3.1416*sqr(z_2)/(atan(y_1/z_2) + atan(yc[ii]/z_2) 
					+ 1.0/(z_2/y_1 + y_1/z_2) + 1.0/(z_2/yc[ii] + yc[ii]/z_2)));// for right wall
		
			// to get mean Rxy/k and Rxz/k near the walls
			// decomposepar not working due to this loop
			if (yc[ii] <= (cell_Y.value() + small) or yc[ii] <= (cell_Y.value() - small))// small = 1e-6
			{
				sum1 += mag(R[ii].component(1)/k_[ii]);//Rxy/k
				count1 += 1.0;
			}
			if (zc[ii] <= (cell_Z.value() + small) or zc[ii] <= (cell_Z.value() - small))
			{
				sum2 += mag(R[ii].component(2)/k_[ii]);//Rxz/k
				count2 += 1.0;
			}
		}		
		const volScalarField Lt1 = (pow(Cmu_, 0.75)/kappa_)*(pow(k_, 1.5)/epsilon_);// lurbulent length scale Lt1

		//Calculate Lt2 (bed, top wall - for ducts), Lt3 (side walls)
		scalar Lt2_factor = sum1/count1;// mean Rxy/k
		scalar Lt3_factor = sum2/count2;// mean Rxz/k
		
		const volScalarField Lt2 = (1.0/kappa_)*(pow(k_, 1.5)/epsilon_)*pow(Lt2_factor, 1.5);
		const volScalarField Lt3 = (1.0/kappa_)*(pow(k_, 1.5)/epsilon_)*pow(Lt3_factor, 1.5);
		
		const volScalarField fw2 = sqr(Lt2/y2);// nonlinear damping function for bottom wall
		const volScalarField fw3 = sqr(Lt3/z1);// nonlinear damping function for left wall 
		const volScalarField fw4 = sqr(Lt3/z2);// nonlinear damping function for right wall 
		
		//nonlinear damping function fw1 is applied depending on the topSurface patch condition: wall (duct flows) or symmetryPlane (free surface flows)
		
		const label& topBoundary = this->mesh_.boundary().findPatchID("topSurface");// top boundary should be named as "topSurface" in the blockMeshDict
		const fvPatch& topBoundaryPatch = this->mesh_.boundary()[topBoundary];
		volScalarField fw1("fw1", fw2);// initialization
		if (isA<wallFvPatch>(topBoundaryPatch))
		{
			fw1 = sqr(Lt2/y1);// nonlinear damping function for top wall - in case of ducts
		}
		else
		{
			fw1 = sqr(Lt1/(y1 + 0.16*Lt1));// nonlinear damping function for free surface
		}
		
        const volSymmTensorField reflect // modified
        (
			(Cref1_*(epsilon_/k_))*R - (Cref2_*((C2_ + 8.0)/11.0))*dev(P)
			- (Cref2_*((8.0*C2_ - 2.0)/11.0))*(D - P + dev(P))
			- ((Cref2_*((30.0*C2_ - 2.0)/55.0))*k_*2.0)*S
        );
		REqn.ref() +=
            (3.0*alpha*rho)
           *(fw1*(dev(symm((n1 & reflect)*n1))) + fw2*(dev(symm((n2 & reflect)*n2))) 
		   + fw3*(dev(symm((n3 & reflect)*n3))) + fw4*(dev(symm((n4 & reflect)*n4))));
    }

    REqn.ref().relax();
    fvConstraints.constrain(REqn.ref());
    solve(REqn);
    fvConstraints.constrain(R);
	
	this->boundNormalStress(R);
	
    k_ = 0.5*tr(R);
	
    correctNut();

    // Correct wall shear-stresses when applying wall-functions
    this->correctWallShearStress(R);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
