/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  dev
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0]; // turbulent viscosity in m^2/s

internalField   uniform 0;

boundaryField
{
    #includeEtc "caseDicts/setConstraintTypes"
	
    lowerWall
    {
        type            	nutUWallFunction;
        value           	uniform 0;
    }
	
    topSurface
    {
		type            codedFixedValue;//changing with kf and epsilonf
		patchType		symmetryPlane;
		value           uniform 0;
		name		    mynutFS;// Name of generated boundary condition, previously directType
		
		code // applicable for half width of a rectangular or square channel
		#{
			const fvPatch& boundaryPatch = this->patch();
			const vectorField& Cf = boundaryPatch.Cf();// face center coordinates for the patch
			const fvPatchSymmTensorField& Rp = boundaryPatch.lookupPatchField<volSymmTensorField, symmTensor>("R");// calling R for the patch
			const fvPatchScalarField& epsilonp = boundaryPatch.lookupPatchField<volScalarField, scalar>("epsilon");// calling epsilon for the patch
			fvPatchScalarField& nut = *this;
			
			forAll(Cf, facei)
			{
				scalar kf = 0.5*tr(Rp[facei]);// turbulent kinetic energy at the facei
				scalar epsilonf = epsilonp[facei];// dissipation rate of turbulent kinetic energy at the facei
				nut[facei] = 0.09*sqr(kf)/epsilonf;
			}
			// nut (free surface) = Cmu_*k(free surface)^2/epsilon(free surface);
		#};

		codeOptions
			#{
				-I$(LIB_SRC)/finiteVolume/lnInclude\
				-I$(LIB_SRC)/meshTools/lnInclude\
			#};
			
		codeInclude
			#{
				#include "fvCFD.H"
				#include <cmath>
				#include <iostream>
			#};
			
		codeLibs 
			#{ -lfiniteVolume -lmeshTools 
			#};
    }
	
    leftWall
    {
		type 				nutUWallFunction;
        value           	uniform 0;
    }
	
	rightSymmetry
    {
        type 				symmetryPlane;
    }
}
// ************************************************************************* //
