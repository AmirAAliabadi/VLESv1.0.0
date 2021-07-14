/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    powerLawTurbulentVelocityFvPatchVectorField

Description
    Boundary condition specifies a powerLawTurbulent velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the powerLawTurbulent coordinate y

SourceFiles
    powerLawTurbulentVelocityFvPatchVectorField.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef powerLawTurbulentVelocityFvPatchVectorField_H
#define powerLawTurbulentVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Random.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class powerLawTurbulentVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class powerLawTurbulentVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    
// Private data


        //- a*y^alpha = U
        scalar a_;
	
	//- Exponent
        scalar alpha_;

	//- zref
        scalar zref_;

	//- Vortex Method: No. of vortex points
        scalar N_; //Should initialize scalars before vectors - otherwise warning!
	//Random N_;

        //- Flow direction
        vector n_;

        //- Direction of the y-coordinate
        vector y_;

	//- Maximum turbulence intensity
	scalar IuMax_;	

	//- Multiplier controlling energy-containing eddy or vortex size
	scalar aSigma_;

	//- Multiplyer constant for characteristic life time of energy-containing eddies or vortices
	scalar aTau_;


public:

    //- Runtime type information
    TypeName("powerLawTurbulentVelocity");


    // Constructors

        //- Construct from patch and internal field
        powerLawTurbulentVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        powerLawTurbulentVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given powerLawTurbulentVelocityFvPatchVectorField
        //  onto a new patch
        powerLawTurbulentVelocityFvPatchVectorField
        (
            const powerLawTurbulentVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new powerLawTurbulentVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        powerLawTurbulentVelocityFvPatchVectorField
        (
            const powerLawTurbulentVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new powerLawTurbulentVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
        
        //- Return a
        scalar& a()
        {
            return a_;
        }

	//- Return alpha
        scalar& alpha()
        {
            return alpha_;
        }

	//- Return zref
        scalar& zref()
        {
            return zref_;
        }

	//- Return N
        scalar& N()
        {
            return N_;
        }

        //- Return flow direction
        vector& n()
        {
            return n_;
        }

        //- Return y direction
        vector& y()
        {
            return y_;
        }

	//- Return IuMax
        scalar& IuMax()
        {
            return IuMax_;
        }

	//- Return aSigma
        scalar& aSigma()
        {
            return aSigma_;
        }

	//- Return aTau
        scalar& aTau()
        {
            return aTau_;
        }
	
        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
