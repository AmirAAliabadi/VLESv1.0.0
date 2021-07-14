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

\*---------------------------------------------------------------------------*/

#include "powerLawTurbulentVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include <iostream>
#include <stdlib.h>
#include <time.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * *  Private Member Functions  *  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

powerLawTurbulentVelocityFvPatchVectorField::powerLawTurbulentVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    a_(0),
    alpha_(0),
    zref_(0),
    N_(0),
    n_(1, 0, 0),
    y_(0, 1, 0),
    IuMax_(0),
    aSigma_(0),
    aTau_(0)
{
}

powerLawTurbulentVelocityFvPatchVectorField::powerLawTurbulentVelocityFvPatchVectorField
(
    const powerLawTurbulentVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    a_(ptf.a_),
    alpha_(ptf.alpha_),
    zref_(ptf.zref_),
    N_(ptf.N_),
    n_(ptf.n_),
    y_(ptf.y_),
    IuMax_(ptf.IuMax_),
    aSigma_(ptf.aSigma_),
    aTau_(ptf.aTau_)
{}


powerLawTurbulentVelocityFvPatchVectorField::powerLawTurbulentVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    a_(readScalar(dict.lookup("a"))),
    alpha_(readScalar(dict.lookup("alpha"))),
    zref_(readScalar(dict.lookup("zref"))),
    N_(readScalar(dict.lookup("N"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y")),
    IuMax_(readScalar(dict.lookup("IuMax"))),
    aSigma_(readScalar(dict.lookup("aSigma"))),
    aTau_(readScalar(dict.lookup("aTau")))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("powerLawTurbulentVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
        //Info << "******IF******" << endl;
    }
    else
    {
        evaluate(Pstream::blocking);
        //Info << "******ELSE******" << endl;
    }
}


powerLawTurbulentVelocityFvPatchVectorField::powerLawTurbulentVelocityFvPatchVectorField
(
    const powerLawTurbulentVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    a_(fcvpvf.a_),
    alpha_(fcvpvf.alpha_),
    zref_(fcvpvf.zref_),
    N_(fcvpvf.N_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_),
    IuMax_(fcvpvf.IuMax_),
    aSigma_(fcvpvf.aSigma_),
    aTau_(fcvpvf.aTau_)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//------------------------------------------------------------------------------
// Mean velocity and kinetic energy in all faces on inlet patch
//------------------------------------------------------------------------------
void powerLawTurbulentVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

// Get range of inlet patch and face coordinates, c
boundBox bb(patch().patch().localPoints(), true);
const vectorField& c = patch().Cf();

// Global y coordinate at face centers of inlet patch
scalarField coord = (c & y_);

// Power law profile at inlet with 0 value for v and w components
vectorField Uinf = (n_*a_*pow(coord,alpha_));
scalarField Uinfx = Uinf.component(0); //get u component of vector field U

// Roughenss lenght z0 and turbulent kinetic energy k
scalar z0 = zref_/exp(1/alpha_);
scalarField Iu = 1/log(coord/z0);
// Should limit z>z0 or equivalently Iu.
forAll(coord,value)
{
	if (Iu[value] >= IuMax_)
	{
	Iu[value]=IuMax_;
	}
}
scalarField k = 1.5*(pow(Uinfx,2)*pow(Iu,2));

//------------------------------------------------------------------------------
// Tangential velocity fluctuations in all faces on inlet patch
//------------------------------------------------------------------------------


if (N_ >= c.size())
{
Info << "Error: N should be less than or equal to  " << c.size() << endl;
//Should be implemented so you could pick arbitrary number of vortex points and not only up to a limit equal to the face numbers at the inlet.
}


// Inlet patch area Ap
scalarField Af = patch().magSf(); //Vector containing area of each face in patch
scalar Ap = sum(Af); 

//Calculate turbulent length scale and sigma from patch geometry
vector bbMax = bb.max();
scalar lx = bbMax[1];
scalar ly = bbMax[2];
scalar L = (2*lx*ly)/(lx+ly);
scalar sigma = (aSigma_*L)/2; 
scalar sigma2 = 2*sqr(sigma);

//Calculate characteristic life time for the energy-containing eddies
//Use dynamic viscosity of air
scalar U0 = a_ * pow(zref_,alpha_);
scalar l0 = sigma;
scalar Rel0 = U0 * l0 / 0.000015;
scalar eta = l0 * pow(Rel0,-3.0/4.0);
scalar ueta = U0 * pow(Rel0,-1.0/4.0);
scalar epsilon = 0.000015 * pow(ueta,2.0) / pow(eta,2.0);
scalar tau0 = aTau_ * pow(l0,2.0/3.0) / pow(epsilon,1.0/3.0);

//Calculations necessary for sampling a new set of vortices
double dt = this->db().time().deltaTValue();
int NvortexSampleIter = floor(tau0 / dt);
if (NvortexSampleIter == 0)
{
	NvortexSampleIter = 1;
}
int timeIteration = floor(this->db().time().timeOutputValue()/this->db().time().deltaTValue());

//Sample a new set of vortices every fixed number of iterations
if (timeIteration % NvortexSampleIter == 0) 		
{
	vectorField ux = c*0; 	//Initialize vector of size = size c
	forAll(patch(),facei)	//Loop over all faces on inlet patch
	{
	vector x1 = c[facei]; 	//Pick coordinates of first face on inlet patch
	//Info << "x1" << x1 << endl;
	scalar x1x = x1[1]; 	//Pick x coordinate. How to make more general?
	scalar x1y = (x1 & y_); //Pick y coordinate		

	srand (time(NULL) ); 	//initialize the random seed
	vector ux_sum = vector(0, 0, 0);

	for (int iter = 1; iter <= N_; iter ++)
	{
		//Info << "iter = " << iter << endl;

		int RandIndex = rand() % c.size(); 	//Generate random number between 1 and 	size c
		//Info << "RandIndex = " << RandIndex << endl;
		vector xi = c[RandIndex]; 		//Picks random coordinates of inlet patch
		//Info << "xi" << xi << endl;
		// xi cannot be equal to x1
		if (xi == x1) 
		{
			scalar Lxi = sqrt(patch().magSf()[iter]); //sqrt of area of random face
			scalar addxi = 0.1*Lxi; //Move vortex point 10% of face Lxi in x and y (no contribution from vortex!)
			xi[1] = xi[1]+addxi;
			xi[2] = xi[2]+addxi;

			//Info << "new Randindex=" << RandIndex << endl;
			//Info << "new xi" << xi << endl;
		}

		scalar xix = xi[1]; 			//How to make more general
		scalar xiy = (xi & y_);

		scalar dist = sqr(x1x-xix)+sqr(x1y-xiy); //Square of distance between face 1 and random face 
		//Info << "dist = " << dist << endl;
		vector xx = xi-x1; //random face position - face 1 position
		//Info << "xx = " << xx << endl;

		scalar ks = k[RandIndex]; //Pick k value for the random face
		//Calculate the cirkulation, ghe, for the random face
		scalar gamma = 4*sqrt((Foam::constant::mathematical::pi*Ap*ks)/(3*N_*0.117783035656383));
		//Info << "gamma = " << ghe << endl;

		//Bound sigma to the local face size
		scalar delta = sqrt(patch().magSf()[RandIndex]); //Delta calculated from face area
		if (sigma <= delta)
		{
			sigma = delta; 
			sigma2 = 2*sqr(sigma);
		}

		//Calculate tangential fluctuations in one face for one vortex point
		//z1 is the streamwise direction of the flow, i.e. x in most cases
		vector z1 = vector(1,0,0);
		//Must randomize the sign for each vortex, 
		//Otherwise all vortices will have the same contribution to velocity near the walls
		//Note that the characteristic time for vortex life is not yet implemented  
		int sign=rand() % 100;
		if (sign > 50)
		{
			sign=1;
		}
		else
		{
			sign=-1;
		}
		//Calculate a local vortex size based on the 1/l(y)=1/lmax+1/(kappa(y+y0)) model
		//This model provides a length scale bounded by ky0 and lmax
		double sigma2facei=sigma2*0.41*(x1y+z0)/(0.41*(x1y+z0)+sigma2);
		vector ux1 = sign*(gamma*((xx ^ z1)/dist)*(1-exp(-dist/sigma2facei))*exp(-dist/sigma2facei));
		//Info << "ux1 = " << ux1 << endl;

		//Sum up the contribution to tangential fluctuations from all vortex points
		ux_sum += ux1;
		//Info << "ux_sum = " << ux_sum << endl;
		} 

		//Store tangential fluctuations for all faces:
		ux[facei] = (1/(2*Foam::constant::mathematical::pi))*ux_sum;
		//Info << "ux = " << ux << endl;
	}

	//vectorField::operator=(n_*a_*pow(coord,alpha_));
	vectorField::operator=(Uinf+ux);
	}
}

// Write
void powerLawTurbulentVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("a")
        << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha")
        << alpha_ << token::END_STATEMENT << nl;
    os.writeKeyword("zref")
        << zref_ << token::END_STATEMENT << nl;
    os.writeKeyword("N")
        << N_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    os.writeKeyword("IuMax")
        << IuMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("aSigma")
        << aSigma_ << token::END_STATEMENT << nl;
    os.writeKeyword("aTau")
        << aTau_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, powerLawTurbulentVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
