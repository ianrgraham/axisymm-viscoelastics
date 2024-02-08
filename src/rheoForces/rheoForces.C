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

#include "rheoForces.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "kinematicMomentumTransportModel.H"
#include "dynamicMomentumTransportModel.H"
#include "phaseKinematicMomentumTransportModel.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rheoForces, 0);

    addToRunTimeSelectionTable(functionObject, rheoForces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::rheoForces::createFileNames
(
    const dictionary& dict
) const
{
    DynamicList<word> names(1);

    const word forceType(dict.lookup("type"));

    // Name for file(fileID::mainFile=0)
    names.append(forceType);

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));
        const label nb = binDict.lookup<label>("nBin");
        if (nb > 0)
        {
            // Name for file(fileID::binsFile=1)
            names.append(forceType + "_bins");
        }
    }

    return move(names);
}


void Foam::functionObjects::rheoForces::writeFileHeader(const label i)
{
    const word forceTypes
    (
        porosity_
      ? "(pressure viscous porous polymer)"
      : "(pressure viscous polymer)"
    );

    switch (fileID(i))
    {
        case fileID::mainFile:
        {
            // force data

            writeHeader(file(i), "Forces");
            writeHeaderValue(file(i), "CofR", coordSys_.origin());
            writeCommented(file(i), "Time");

            file(i)
                << "forces" << forceTypes << tab
                << "moments" << forceTypes;

            if (localSystem_)
            {
                file(i)
                    << tab
                    << "localForces" << forceTypes << tab
                    << "localMoments" << forceTypes;
            }

            break;
        }
        case fileID::binsFile:
        {
            // bin data

            writeHeader(file(i), "Force bins");
            writeHeaderValue(file(i), "bins", nBin_);
            writeHeaderValue(file(i), "start", binMin_);
            writeHeaderValue(file(i), "delta", binDx_);
            writeHeaderValue(file(i), "direction", binDir_);

            vectorField binPoints(nBin_);
            writeCommented(file(i), "x co-ords  :");
            forAll(binPoints, pointi)
            {
                binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
                file(i) << tab << binPoints[pointi].x();
            }
            file(i) << nl;

            writeCommented(file(i), "y co-ords  :");
            forAll(binPoints, pointi)
            {
                file(i) << tab << binPoints[pointi].y();
            }
            file(i) << nl;

            writeCommented(file(i), "z co-ords  :");
            forAll(binPoints, pointi)
            {
                file(i) << tab << binPoints[pointi].z();
            }
            file(i) << nl;

            writeCommented(file(i), "Time");

            for (label j = 0; j < nBin_; j++)
            {
                const word jn('(' + Foam::name(j) + ')');
                const word f("rheoForces" + jn + forceTypes);
                const word m("moments" + jn + forceTypes);

                file(i)<< tab << f << tab << m;
            }
            if (localSystem_)
            {
                for (label j = 0; j < nBin_; j++)
                {
                    const word jn('(' + Foam::name(j) + ')');
                    const word f("localForces" + jn + forceTypes);
                    const word m("localMoments" + jn + forceTypes);

                    file(i)<< tab << f << tab << m;
                }
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }

    file(i)<< endl;
}


void Foam::functionObjects::rheoForces::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!obr_.foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database."
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !obr_.foundObject<volVectorField>(UName_)
         || !obr_.foundObject<volScalarField>(pName_)

        )
        {
            FatalErrorInFunction
                << "Could not find " << UName_ << ", " << pName_
                << exit(FatalError);
        }

        if
        (
            rhoName_ != "rhoInf"
         && !obr_.foundObject<volScalarField>(rhoName_)
        )
        {
            FatalErrorInFunction
                << "Could not find " << rhoName_
                << exit(FatalError);
        }
    }

    initialised_ = true;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::rheoForces::devTau() const
{
    typedef incompressible::momentumTransportModel icoModel;
    typedef compressible::momentumTransportModel cmpModel;
    typedef phaseIncompressible::momentumTransportModel phaseIcoModel;
    typedef phaseCompressible::momentumTransportModel phaseCmpModel;

    const word& modelName = momentumTransportModel::typeName;
    const word phaseModelName =
        phaseName_ == word::null
      ? word::null
      : IOobject::groupName(momentumTransportModel::typeName, phaseName_);

    if (obr_.foundObject<icoModel>(modelName))
    {
        const incompressible::momentumTransportModel& model =
            obr_.lookupObject<icoModel>(modelName);

        return alpha()*rho()*model.devSigma();
    }
    else if (obr_.foundObject<cmpModel>(modelName))
    {
        const cmpModel& model =
            obr_.lookupObject<cmpModel>(modelName);

        return alpha()*model.devTau();
    }
    else if (obr_.foundObject<phaseIcoModel>(phaseModelName))
    {
        const phaseIcoModel& model =
            obr_.lookupObject<phaseIcoModel>(phaseModelName);

        return rho()*model.devSigma();
    }
    else if (obr_.foundObject<phaseCmpModel>(phaseModelName))
    {
        const phaseCmpModel& model =
            obr_.lookupObject<phaseCmpModel>(phaseModelName);

        return model.devTau();
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        // Legacy support for icoFoam

        const dictionary& transportProperties =
             obr_.lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::rheoForces::mu() const
{
    typedef incompressible::momentumTransportModel icoModel;
    typedef compressible::momentumTransportModel cmpModel;
    typedef phaseIncompressible::momentumTransportModel phaseIcoModel;
    typedef phaseCompressible::momentumTransportModel phaseCmpModel;

    const word& modelName = momentumTransportModel::typeName;
    const word phaseModelName =
        phaseName_ == word::null
      ? word::null
      : IOobject::groupName(momentumTransportModel::typeName, phaseName_);

    if (obr_.foundObject<icoModel>(modelName))
    {
        const incompressible::momentumTransportModel& model =
            obr_.lookupObject<icoModel>(modelName);

        return rho()*model.transport().nu();
    }
    else if (obr_.foundObject<cmpModel>(modelName))
    {
        const cmpModel& model =
            obr_.lookupObject<cmpModel>(modelName);

        return model.transport().mu();
    }
    else if (obr_.foundObject<phaseIcoModel>(phaseModelName))
    {
        const phaseIcoModel& model =
            obr_.lookupObject<phaseIcoModel>(phaseModelName);

        return rho()*model.transport().nu();
    }
    else if (obr_.foundObject<phaseCmpModel>(phaseModelName))
    {
        const phaseCmpModel& model =
            obr_.lookupObject<phaseCmpModel>(phaseModelName);

        return model.transport().mu();
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        // Legacy support for icoFoam

        const dictionary& transportProperties =
             obr_.lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::rheoForces::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return volScalarField::New
        (
            "rho",
            mesh_,
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }
    else
    {
        return(obr_.lookupObject<volScalarField>(rhoName_));
    }
}


Foam::scalar Foam::functionObjects::rheoForces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorInFunction
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::rheoForces::alpha() const
{
    if (phaseName_ == word::null)
    {
        return volScalarField::New
        (
            "alpha",
            mesh_,
            dimensionedScalar(dimless, 1)
        );
    }
    else
    {
        return obr_.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phaseName_)
        );
    }
}


Foam::tmp<Foam::scalarField> Foam::functionObjects::rheoForces::alpha
(
    const label patchi
) const
{
    if (phaseName_ == word::null)
    {
        return tmp<scalarField>
        (
            new scalarField(mesh_.boundary()[patchi].size(), 1)
        );
    }
    else
    {
        return obr_.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phaseName_)
        ).boundaryField()[patchi];
    }
}


void Foam::functionObjects::rheoForces::applyBins
(
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fPor,
    const vectorField& fPol,
    const vectorField& d
)
{
    if (nBin_ == 1)
    {
        // Info << "0" << endl;
        force_[0][0] += sum(fN); // pressure
        // Info << "fN" << endl;
        force_[1][0] += sum(fT); // viscous
        // Info << "fT" << endl;
        force_[2][0] += sum(fPor); // porous
        // Info << "fPor" << endl;
        force_[3][0] += sum(fPol);
        // Info << "fPol" << endl;
        // Info << "1" << endl;
        moment_[0][0] += sum(Md^fN);
        // Info << "mN" << endl;
        moment_[1][0] += sum(Md^fT);
        // Info << "mT" << endl;
        moment_[2][0] += sum(Md^fPor);
        moment_[3][0] += sum(Md^fPol);
    }
    else
    {
        scalarField dd((d & binDir_) - binMin_);

        forAll(dd, i)
        {
            label bini = min(max(floor(dd[i]/binDx_), 0), force_[0].size() - 1);

            force_[0][bini] += fN[i];
            force_[1][bini] += fT[i];
            force_[2][bini] += fPor[i];
            force_[3][bini] += fPol[i];
            moment_[0][bini] += Md[i]^fN[i];
            moment_[1][bini] += Md[i]^fT[i];
            moment_[2][bini] += Md[i]^fPor[i];
            moment_[3][bini] += Md[i]^fPol[i];
        }
    }
}


void Foam::functionObjects::rheoForces::writeForces()
{
    Log << type() << " " << name() << " write:" << nl
        << "    sum of rheoForces:" << nl
        << "        pressure : " << sum(force_[0]) << nl
        << "        viscous  : " << sum(force_[1]) << nl
        << "        porous   : " << sum(force_[2]) << nl
        << "        polymer  : " << sum(force_[3]) << nl
        << "    sum of moments:" << nl
        << "        pressure : " << sum(moment_[0]) << nl
        << "        viscous  : " << sum(moment_[1]) << nl
        << "        porous   : " << sum(moment_[2]) << nl
        << "        polymer  : " << sum(moment_[3]) << nl
        << endl;

    writeTime(file(fileID::mainFile));

    if (porosity_)
    {
        file(fileID::mainFile) << tab << setw(1) << '('
            << sum(force_[0]) << setw(1) << ' '
            << sum(force_[1]) << setw(1) << ' '
            << sum(force_[2]) << setw(1) << ' '
            << sum(force_[3]) << setw(3) << ") ("
            << sum(moment_[0]) << setw(1) << ' '
            << sum(moment_[1]) << setw(1) << ' '
            << sum(moment_[2]) << setw(1) << ' '
            << sum(moment_[3]) << setw(1) << ')';
    }
    else
    {
        file(fileID::mainFile) << tab << setw(1) << '('
            << sum(force_[0]) << setw(1) << ' '
            << sum(force_[1]) << setw(1) << ' '
            << sum(force_[3]) << setw(3) << ") ("
            << sum(moment_[0]) << setw(1) << ' '
            << sum(moment_[1]) << setw(1) << ' '
            << sum(moment_[3]) << setw(1) << ')';
    }

    if (localSystem_)
    {
        if (porosity_)
        {
            file(fileID::mainFile) << tab << setw(1) << '('
                << sum(coordSys_.localVector(force_[0])) << setw(1) << ' '
                << sum(coordSys_.localVector(force_[1])) << setw(1) << ' '
                << sum(coordSys_.localVector(force_[2])) << setw(1) << ' '
                << sum(coordSys_.localVector(force_[3])) << setw(3) << ") ("
                << sum(coordSys_.localVector(moment_[0])) << setw(1) << ' '
                << sum(coordSys_.localVector(moment_[1])) << setw(1) << ' '
                << sum(coordSys_.localVector(moment_[2])) << setw(1) << ' '
                << sum(coordSys_.localVector(moment_[3])) << setw(1) << ')';
        }
        else
        {
            file(fileID::mainFile) << tab << setw(1) << '('
                << sum(coordSys_.localVector(force_[0])) << setw(1) << ' '
                << sum(coordSys_.localVector(force_[1])) << setw(1) << ' '
                << sum(coordSys_.localVector(force_[3])) << setw(3) << ") ("
                << sum(coordSys_.localVector(moment_[0])) << setw(1) << ' '
                << sum(coordSys_.localVector(moment_[1])) << setw(1) << ' '
                << sum(coordSys_.localVector(moment_[3])) << setw(1) << ')';
        }
    }

    file(fileID::mainFile) << endl;
}


void Foam::functionObjects::rheoForces::writeBins()
{
    if (nBin_ == 1)
    {
        return;
    }

    List<Field<vector>> f(force_);
    List<Field<vector>> m(moment_);

    if (binCumulative_)
    {
        for (label i = 1; i < f[0].size(); i++)
        {
            f[0][i] += f[0][i-1];
            f[1][i] += f[1][i-1];
            f[2][i] += f[2][i-1];
            f[3][i] += f[3][i-1];

            m[0][i] += m[0][i-1];
            m[1][i] += m[1][i-1];
            m[2][i] += m[2][i-1];
            m[3][i] += m[3][i-1];
        }
    }

    writeTime(file(fileID::binsFile));


    forAll(f[0], i)
    {
        if (porosity_)
        {
            file(fileID::binsFile)
                << tab << setw(1) << '('
                << f[0][i] << setw(1) << ' '
                << f[1][i] << setw(1) << ' '
                << f[2][i] << setw(1) << ' '
                << f[3][i] << setw(3) << ") ("
                << m[0][i] << setw(1) << ' '
                << m[1][i] << setw(1) << ' '
                << m[2][i] << setw(1) << ' '
                << m[3][i] << setw(1) << ')';
        }
        else
        {
            file(fileID::binsFile)
                << tab << setw(1) << '('
                << f[0][i] << setw(1) << ' '
                << f[1][i] << setw(1) << " "
                << f[3][i] << setw(3) << ") ("
                << m[0][i] << setw(1) << ' '
                << m[1][i] << setw(1) << ' '
                << m[3][i] << setw(1) << ')';
        }
    }

    if (localSystem_)
    {
        List<Field<vector>> lf(4);
        List<Field<vector>> lm(4);
        lf[0] = coordSys_.localVector(force_[0]);
        lf[1] = coordSys_.localVector(force_[1]);
        lf[2] = coordSys_.localVector(force_[2]);
        lf[3] = coordSys_.localVector(force_[3]);
        lm[0] = coordSys_.localVector(moment_[0]);
        lm[1] = coordSys_.localVector(moment_[1]);
        lm[2] = coordSys_.localVector(moment_[2]);
        lm[3] = coordSys_.localVector(moment_[3]);

        if (binCumulative_)
        {
            for (label i = 1; i < lf[0].size(); i++)
            {
                lf[0][i] += lf[0][i-1];
                lf[1][i] += lf[1][i-1];
                lf[2][i] += lf[2][i-1];
                lf[3][i] += lf[3][i-1];
                lm[0][i] += lm[0][i-1];
                lm[1][i] += lm[1][i-1];
                lm[2][i] += lm[2][i-1];
                lm[3][i] += lm[3][i-1];
            }
        }

        forAll(lf[0], i)
        {
            if (porosity_)
            {
                file(fileID::binsFile)
                    << tab << setw(1) << '('
                    << lf[0][i] << setw(1) << ' '
                    << lf[1][i] << setw(1) << ' '
                    << lf[2][i] << setw(1) << ' '
                    << lf[3][i] << setw(3) << ") ("
                    << lm[0][i] << setw(1) << ' '
                    << lm[1][i] << setw(1) << ' '
                    << lm[2][i] << setw(1) << ' '
                    << lm[3][i] << setw(1) << ')';
            }
            else
            {
                file(fileID::binsFile)
                    << tab << setw(1) << '('
                    << lf[0][i] << setw(1) << ' '
                    << lf[1][i] << setw(1) << " "
                    << lf[3][i] << setw(3) << ") ("
                    << lm[0][i] << setw(1) << ' '
                    << lm[1][i] << setw(1) << ' '
                    << lm[3][i] << setw(1) << ')';
            }
        }
    }

    file(fileID::binsFile) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rheoForces::rheoForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    force_(4),
    moment_(4),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    phaseName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(vGreat),
    pRef_(0),
    coordSys_("coordSys", vector::zero),
    localSystem_(false),
    porosity_(false),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(great),
    binPoints_(),
    binCumulative_(true),
    initialised_(false)
{
    read(dict);
}


Foam::functionObjects::rheoForces::rheoForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, obr, dict),
    logFiles(obr_, name),
    force_(4),
    moment_(4),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    phaseName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(vGreat),
    pRef_(0),
    coordSys_("coordSys", vector::zero),
    localSystem_(false),
    porosity_(false),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(great),
    binPoints_(),
    binCumulative_(true),
    initialised_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::rheoForces::~rheoForces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::rheoForces::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    initialised_ = false;

    Log << type() << " " << name() << ":" << nl;

    directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

    if (directForceDensity_)
    {
        // Optional entry for fDName
        fDName_ = dict.lookupOrDefault<word>("fD", "fD");
    }
    else
    {
        // Optional phase entry
        phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

        // Optional U, p and rho entries
        pName_ =
            dict.lookupOrDefault<word>
            (
                "p",
                IOobject::groupName("p", phaseName_)
            );
        UName_ =
            dict.lookupOrDefault<word>
            (
                "U",
                IOobject::groupName("U", phaseName_)
            );
        rhoName_ =
            dict.lookupOrDefault<word>
            (
                "rho",
                IOobject::groupName("rho", phaseName_)
            );

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            dict.lookup("rhoInf") >> rhoRef_;
        }

        // Reference pressure, 0 by default
        pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);
    }

    // Centre of rotation for moment calculations
    // specified directly, from coordinate system, or implicitly (0 0 0)
    if (dict.found("CofR"))
    {
        coordSys_ = coordinateSystem("coordSys", vector(dict.lookup("CofR")));
        localSystem_ = false;
    }
    else
    {
        coordSys_ = coordinateSystem("coordSys", dict);
        localSystem_ = true;
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Log << "    Including porosity effects" << endl;
    }
    else
    {
        Log << "    Not including porosity effects" << endl;
    }

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));
        binDict.lookup("nBin") >> nBin_;

        if (nBin_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Number of bins (nBin) must be zero or greater"
                << exit(FatalIOError);
        }
        else if ((nBin_ == 0) || (nBin_ == 1))
        {
            nBin_ = 1;
            forAll(force_, i)
            {
                force_[i].setSize(1);
                moment_[i].setSize(1);
            }
        }

        if (nBin_ > 1)
        {
            binDict.lookup("direction") >> binDir_;
            binDir_ /= mag(binDir_);

            binMin_ = great;
            scalar binMax = -great;
            forAllConstIter(labelHashSet, patchSet_, iter)
            {
                const label patchi = iter.key();
                const polyPatch& pp = pbm[patchi];
                const scalarField d(pp.faceCentres() & binDir_);
                binMin_ = min(min(d), binMin_);
                binMax = max(max(d), binMax);
            }
            reduce(binMin_, minOp<scalar>());
            reduce(binMax, maxOp<scalar>());

            // slightly boost binMax so that region of interest is fully
            // within bounds
            binMax = 1.0001*(binMax - binMin_) + binMin_;

            binDx_ = (binMax - binMin_)/scalar(nBin_);

            // create the bin points used for writing
            binPoints_.setSize(nBin_);
            forAll(binPoints_, i)
            {
                binPoints_[i] = (i + 0.5)*binDir_*binDx_;
            }

            binDict.lookup("cumulative") >> binCumulative_;

            // allocate storage for rheoForces and moments
            forAll(force_, i)
            {
                force_[i].setSize(nBin_);
                moment_[i].setSize(nBin_);
            }
        }
    }

    if (nBin_ == 1)
    {
        // allocate storage for rheoForces and moments
        force_[0].setSize(1);
        force_[1].setSize(1);
        force_[2].setSize(1);
        force_[3].setSize(1);
        moment_[0].setSize(1);
        moment_[1].setSize(1);
        moment_[2].setSize(1);
        moment_[3].setSize(1);
    }

    resetNames(createFileNames(dict));

    return true;
}


void Foam::functionObjects::rheoForces::calcForcesMoment()
{
    initialise();

    force_[0] = Zero;
    force_[1] = Zero;
    force_[2] = Zero;
    force_[3] = Zero;

    moment_[0] = Zero;
    moment_[1] = Zero;
    moment_[2] = Zero;
    moment_[3] = Zero;

    if (directForceDensity_)
    {
        const volVectorField& fD = obr_.lookupObject<volVectorField>(fDName_);

        const surfaceVectorField::Boundary& Sfb =
            mesh_.Sf().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            const label patchi = iter.key();

            const vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - coordSys_.origin()
            );

            const scalarField sA(mag(Sfb[patchi]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            const vectorField fN
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );

            // Tangential force (total force minus normal fN)
            const vectorField fT(sA*fD.boundaryField()[patchi] - fN);

            //- Porous and polymer force
            const vectorField fDummy(Md.size(), Zero);

            applyBins(Md, fN, fT, fDummy, fDummy, mesh_.C().boundaryField()[patchi]);
        }
    }
    else if (mesh_.foundObject<constitutiveModel>("constitutiveProperties"))
    {

        // Info << "Using constitutiveProperties" << endl;
        // get the pressure field
        // Info << "Getting pressure" << endl;
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        // Face areas of the boundary
        // Info << "Getting surface vector field" << endl;
        const surfaceVectorField::Boundary& Sfb =
            mesh_.Sf().boundaryField();

        // Scale pRef by density for incompressible simulations
        const scalar pRef = pRef_/rho(p);
        
        constitutiveModel &constEq_ = const_cast<constitutiveModel &>(
              mesh_.lookupObject<constitutiveModel>("constitutiveProperties"));

        // get the boundary fields for the polymer stress and total stress
        auto tauP = constEq_.tau();
        auto tauTot = constEq_.tauTotal();
        auto devTau = tauTot - tauP;
        // Info << "Getting polymer stress" << endl;
        const volSymmTensorField::Boundary& tauPb = tauP().boundaryField();
        // Info << "Getting total stress" << endl;
        const volSymmTensorField::Boundary& devTaub = devTau().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            // get the patch index
            const label patchi = iter.key();

            const vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - coordSys_.origin()
            );

            // // compute deviatoric shear stress on the boundary patch (minus any polymer contribution)
            // const symmTensorField devTaub
            // (
            //     dev(tauTotb[patchi] - tauPb[patchi])
            // );

            // pressure force
            const vectorField fN
            (
                alpha(patchi)
               *rho(p)
               *Sfb[patchi]
               *(p.boundaryField()[patchi] - pRef)
            );

            // polymer force
            const vectorField fPol(Sfb[patchi] & tauPb[patchi]);

            // viscous force
            // tensor inner product (&)
            const vectorField fT(Sfb[patchi] & devTaub[patchi]);

            const vectorField fDummy(Md.size(), Zero);

            applyBins(Md, fN, fT, fDummy, fPol, mesh_.C().boundaryField()[patchi]);
        }
    }
    else
    {
        // get the pressure field
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        // Face areas of the boundary
        const surfaceVectorField::Boundary& Sfb =
            mesh_.Sf().boundaryField();

        // Viscous stress of the boundary (maybe deviatoric?)
        tmp<volSymmTensorField> tdevTau = devTau();
        const volSymmTensorField::Boundary& devTaub =
            tdevTau().boundaryField();

        // Scale pRef by density for incompressible simulations
        const scalar pRef = pRef_/rho(p);

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            // get the patch index
            const label patchi = iter.key();

            const vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - coordSys_.origin()
            );

            // Info << "rheoForces: rho = " << rho(p) << endl;

            // pressure force
            const vectorField fN
            (
                alpha(patchi)
               *rho(p)
               *Sfb[patchi]
               *(p.boundaryField()[patchi] - pRef)
            );

            // Info << "rheoForces: fN = " << sum(fN) << endl;

            // viscous force
            // tensor inner product
            const vectorField fT(Sfb[patchi] & devTaub[patchi]);

            // Info << "rheoForces: fT = " << sum(fT) << endl;

            const vectorField fDummy(Md.size(), Zero);

            applyBins(Md, fN, fT, fDummy, fDummy, mesh_.C().boundaryField()[patchi]);
        }
    }

    if (porosity_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            const porosityModel& pm = *iter();

            const vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                const label zoneI = cellZoneIDs[i];
                const cellZone& cZone = mesh_.cellZones()[zoneI];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - coordSys_.origin());

                const vectorField fDummy(Md.size(), Zero);

                applyBins(Md, fDummy, fDummy, fP, fDummy, d);
            }
        }
    }

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineScatter(moment_);
}


Foam::vector Foam::functionObjects::rheoForces::forceEff() const
{
    return sum(force_[0]) + sum(force_[1]) + sum(force_[2]) + sum(force_[3]);
}

Foam::vector Foam::functionObjects::rheoForces::forcePress() const
{
    return sum(force_[0]);
}

Foam::vector Foam::functionObjects::rheoForces::forceVisc() const
{
    return sum(force_[1]);
}

Foam::vector Foam::functionObjects::rheoForces::forcePor() const
{
    return sum(force_[2]);
}

Foam::vector Foam::functionObjects::rheoForces::forcePoly() const
{
    return sum(force_[3]);
}

Foam::vector Foam::functionObjects::rheoForces::momentEff() const
{
    return sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2]) + sum(moment_[3]);
}


bool Foam::functionObjects::rheoForces::execute()
{
    return true;
}


bool Foam::functionObjects::rheoForces::write()
{
    calcForcesMoment();

    if (Pstream::master())
    {
        logFiles::write();

        writeForces();

        writeBins();

        Log << endl;
    }

    return true;
}


// ************************************************************************* //
