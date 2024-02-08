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
    Foam::inletFlowFvPatchVectorField

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
    inletFlowFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef inletFlowFvPatchVectorField_H
#define inletFlowFvPatchVectorField_H

#include <deque>
#include <optional>
// #include "lapack.h"
#include "fixedValueFvPatchFields.H"
#include "scalar.H"
#include "vectorField.H"
#include "rheoForcesv2.H"
#include "terminalFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class inletFlowFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class testFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    // Private Data
        scalar u_max_;

        scalar t_interval_;

        scalar r_pipe_;

        vectorField U_;

        scalar cur_tstep_ = -1.0;

        int step = 0;

public:
    TypeName("testInletFlow");
    // Constructors

        //- Construct from patch and internal field
        testFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        testFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given inletFlowFvPatchVectorField
        //  onto a new patch
        testFvPatchVectorField
        (
            const testFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        testFvPatchVectorField(const testFvPatchVectorField&) = delete;

        //- Copy constructor setting internal field reference
        testFvPatchVectorField
        (
            const testFvPatchVectorField&,
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
                new testFvPatchVectorField(*this, iF)
            );
        }

    // Member Functions

        // Evaluation functions

            virtual void updateCoeffs();

            // - Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
            );

        //- Write
        virtual void write(Ostream&) const;

        //- value getter
        vector value() const
        {
            return vector(0, 0, 0);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
