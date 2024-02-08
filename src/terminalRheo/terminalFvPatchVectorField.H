/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::terminalFvPatchVectorField

Description
    This boundary condition supplies a fixed value constraint, and is the base
    class for a number of other boundary conditions.

Usage
    \table
        Property     | Description             | Required    | Default value
        value        | Patch face values       | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            terminal;
        value           uniform 0;  // Example for scalar field usage
    }
    \endverbatim

SourceFiles
    terminalFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef terminalFvPatchVectorField_H
#define terminalFvPatchVectorField_H

#include <deque>
#include <optional>
// #include "lapack.h"
#include "fixedValueFvPatchFields.H"
#include "scalar.H"
#include "vectorField.H"
#include "rheoForcesv2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class terminalFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class terminalFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Current value
        vector value_;

        scalar mass_;

        scalar dt_;

        //- Force name
        const word force_name_;

        //- Force object dictionary
        const dictionary force_dict_;

        scalar tolerance_;

        bool reverse_sign_;

        scalar cur_tstep_ = -1.0;

        int step = 0;

        vector vel0_ = vector(0.0, 0.0, 0.0);

        vector force0_ = vector(0.0, 0.0, 0.0);

        vector force00_ = vector(0.0, 0.0, 0.0);

public:

    //- Runtime type information
    TypeName("terminalRheo");


    // Constructors

        //- Construct from patch and internal field
        terminalFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        terminalFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given terminalFvPatchVectorField
        //  onto a new patch
        terminalFvPatchVectorField
        (
            const terminalFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        terminalFvPatchVectorField(const terminalFvPatchVectorField&) = delete;

        //- Copy constructor setting internal field reference
        terminalFvPatchVectorField
        (
            const terminalFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new terminalFvPatchVectorField(*this, iF)
            );
        }

    // Member Functions

        // Evaluation functions

            virtual void updateCoeffs();

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
            );

        //- Write
        virtual void write(Ostream&) const;

        //- value getter
        vector value() const
        {
            return value_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
