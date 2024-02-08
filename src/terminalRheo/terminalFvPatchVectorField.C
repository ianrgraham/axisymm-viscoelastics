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

#include "terminalFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::terminalFvPatchVectorField::terminalFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    value_(vector(0.0, 0.0, 0.0)),
    mass_(0.75e-6),
    dt_(1e-4),
    force_name_("__force__"),
    force_dict_(Foam::dictionary("forceDict")),
    tolerance_(1e-10),
    reverse_sign_(true)
{
    Field<vector>::operator=(value_);
}



Foam::terminalFvPatchVectorField::terminalFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    value_(dict.lookupOrDefault("init_value", vector(0.0, 0.0, 0.0))),
    mass_(dict.lookupOrDefault("mass", 0.75e-6)),
    dt_(dict.lookupOrDefault("dt", 1e-4)),
    force_name_("__force__"),
    force_dict_(dict.subDict("forceDict")),
    tolerance_(dict.lookupOrDefault("tolerance", 1e-10)),
    reverse_sign_(dict.lookupOrDefault("reverse_sign", true))
{}



Foam::terminalFvPatchVectorField::terminalFvPatchVectorField
(
    const terminalFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    value_(ptf.value_),
    mass_(ptf.mass_),
    dt_(ptf.dt_),
    force_name_(ptf.force_name_),
    force_dict_(ptf.force_dict_),
    tolerance_(ptf.tolerance_),
    reverse_sign_(ptf.reverse_sign_)
{}



Foam::terminalFvPatchVectorField::terminalFvPatchVectorField
(
    const terminalFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    value_(ptf.value_),
    mass_(ptf.mass_),
    dt_(ptf.dt_),
    force_name_(ptf.force_name_),
    force_dict_(ptf.force_dict_),
    tolerance_(ptf.tolerance_),
    reverse_sign_(ptf.reverse_sign_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::terminalFvPatchVectorField::updateCoeffs()
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

    const word& name = this->force_name_;
    const dictionary& dict = this->force_dict_;
    Foam::functionObjects::rheoForcesv2 force(name, runTime, dict);

    force.calcForcesMoment();
    vector eff_force = force.forceEff();

    eff_force.y() = 0.0;
    eff_force.z() = 0.0;

    if (fabs(eff_force.x()) < tolerance_)
    {
        eff_force.x() = 0.0;
    }

    if (reverse_sign_)
    {
        eff_force.x() *= -1.0;
    }

    // integrate velocity to third order implicit multistep method (Adams-Moulton)
    this->value_ = 
        (dt_ / mass_) * (
            (5.0 / 12.0) * eff_force
            + (8.0 / 12.0) * force0_
            - (1.0 / 12.0) * force00_
        )
        + vel0_;

    force00_ = force0_;
    force0_ = eff_force;
    vel0_ = this->value_;

    Field<vector>::operator=(this->value_);

    Info << "value = " << value_ << endl;

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::terminalFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    Info << endl << "Evaluating terminalRheoFvPatchVectorField" << endl;

    fixedValueFvPatchVectorField::evaluate(commsType);
}


void Foam::terminalFvPatchVectorField::write(Ostream& os) const
{
    // fixedValueFvPatchVectorField::write(os);
    writeEntry(os, "type", "terminalRheo");
    writeEntry(os, "value", *this);

    writeEntry(os, "mass", mass_);
    writeEntry(os, "dt", dt_);
    writeEntry(os, "tolerance", tolerance_);
    writeEntry(os, "reverse_sign", reverse_sign_);
    writeEntry(os, "force_name", force_name_);
    writeEntry(os, "forceDict", force_dict_);

    // writeEntry(os, "value", *this);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        terminalFvPatchVectorField
    );
}


// ************************************************************************* //
