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

#include "inletFlowFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::testFvPatchVectorField::testFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    u_max_(0.0),
    t_interval_(1.0),
    r_pipe_(3.7),
    U_(p.size(), vector::zero)
{
    // const fvPatch& boundaryPatch = patch();
    // const vectorField& Cf = boundaryPatch.Cf();
    // vectorField& field = U_;

    // forAll(Cf, faceI)
    // {
    //     field[faceI] = vector(0, 0, 0);
    // }

    fvPatchField<vector>::operator=
    (
        U_
    ); 
}



Foam::testFvPatchVectorField::testFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    u_max_(dict.lookupOrDefault("u_max", 0.0)),
    t_interval_(dict.lookupOrDefault("time_interval", 1.0)),
    r_pipe_(dict.lookupOrDefault("r_pipe", 3.7)),
    U_(p.size(), vector::zero)
{
    // const fvPatch& boundaryPatch = patch();
    // const vectorField& Cf = boundaryPatch.Cf();
    // vectorField& field = U_;

    // forAll(Cf, faceI)
    // {
    //     field[faceI] = vector(0, 0, 0);
    // }

    fvPatchField<vector>::operator=
    (
        U_
    ); 
}



Foam::testFvPatchVectorField::testFvPatchVectorField
(
    const testFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    u_max_(ptf.u_max_),
    t_interval_(ptf.t_interval_),
    r_pipe_(ptf.r_pipe_),
    U_(ptf.U_)
{}



Foam::testFvPatchVectorField::testFvPatchVectorField
(
    const testFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    u_max_(ptf.u_max_),
    t_interval_(ptf.t_interval_),
    r_pipe_(ptf.r_pipe_),
    U_(ptf.U_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::testFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const Time& runTime = db().time();
    scalar time = runTime.value();

    if (cur_tstep_ == time)
    {
        return;
    }
    cur_tstep_ = time;
    step += 1;

    const scalar t = this->db().time().value();

    const volVectorField& U = db().lookupObject<volVectorField>("U");
    const label wall = U.mesh().boundaryMesh().findPatchID("wall");

    const fvPatchVectorField& uwall = U.boundaryField()[wall];
    
    // cast to terminalFvPatchVectorField
    const terminalFvPatchVectorField& U_wall = dynamic_cast<const terminalFvPatchVectorField&>(uwall);
    const scalar Uwall = U_wall.value().x();

    const scalar p1 = this->u_max_;
    const scalar p0 = 0.0; // p1*1e-2;
    const scalar t1 = this->t_interval_;
    const scalar t0 = 0.0;
    const scalar PI = 3.14159265358979323846;

    const scalar delta_p = p1 - p0;
    // const scalar delta_v = v1 - v0;
    const scalar delta_t = t1 - t0;

    const scalar rpipe = this->r_pipe_;

    const fvPatch& boundaryPatch = patch();
    const vectorField& Cf = boundaryPatch.Cf();
    vectorField& field = U_;

    scalar Umax = 0.0;

    if (t < t0) {
        Umax = p0;
    }
    else if (t < t1) {
        auto dt = (t - t0)/delta_t;
        auto out = p0 + delta_p * (1 - cos(dt*PI)) / 2.0;
        Umax = out;
    } else {
        Umax = p1;
    }

    std::cout << "Umax: " << Umax << std::endl;
    std::cout << "Uwall: " << Uwall << std::endl;

    forAll(Cf, faceI)
    {
        const scalar y = Cf[faceI].y();
        const scalar z = Cf[faceI].z();
        const scalar r = sqrt(z*z + y*y);
        const scalar U = Umax * (1 - powf(r/rpipe, 2.0)) + Uwall;
        field[faceI] = vector(U, 0, 0);
    }

    // Field<vector>::operator=(value_);
    vectorField::operator=(field);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::testFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    Info << endl << "Evaluating terminalRheoFvPatchVectorField" << endl;

    fixedValueFvPatchVectorField::evaluate(commsType);
}


void Foam::testFvPatchVectorField::write(Ostream& os) const
{
    // fixedValueFvPatchVectorField::write(os);
    writeEntry(os, "type", "testInletFlow");
    writeEntry(os, "value", *this);

    writeEntry(os, "u_max", u_max_);
    writeEntry(os, "time_interval", t_interval_);
    writeEntry(os, "r_pipe", r_pipe_);
    // writeEntry(os, "value", *this);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        testFvPatchVectorField
    );
}


// ************************************************************************* //
