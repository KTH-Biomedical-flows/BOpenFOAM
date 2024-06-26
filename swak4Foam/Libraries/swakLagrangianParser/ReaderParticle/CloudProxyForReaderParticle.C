/*---------------------------------------------------------------------------*\
|                       _    _  _     ___                       | The         |
|     _____      ____ _| | _| || |   / __\__   __ _ _ __ ___    | Swiss       |
|    / __\ \ /\ / / _` | |/ / || |_ / _\/ _ \ / _` | '_ ` _ \   | Army        |
|    \__ \\ V  V / (_| |   <|__   _/ / | (_) | (_| | | | | | |  | Knife       |
|    |___/ \_/\_/ \__,_|_|\_\  |_| \/   \___/ \__,_|_| |_| |_|  | For         |
|                                                               | OpenFOAM    |
-------------------------------------------------------------------------------
License
    This file is part of swak4Foam.

    swak4Foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    swak4Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with swak4Foam.  If not, see <http://www.gnu.org/licenses/>.

Contributors/Copyright:
    2012-2013, 2016-2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "CloudProxyForReaderParticle.H"

#include "DebugOStream.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(CloudProxyForReaderParticle,0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CloudProxyForReaderParticle::CloudProxyForReaderParticle
(
    const cloud& c
)
:
    CloudProxyForParticle<ReaderParticleCloud>(c),
    readerCloud_(
        dynamicCast<const ReaderParticleCloud>(c)
    )
{
    internalAddFields<scalar>(readerCloud_.getScalarFieldNames());
    internalAddFields<scalar>(readerCloud_.getLabelFieldNames());
    internalAddFields<vector>(readerCloud_.getVectorFieldNames());
    internalAddFields<tensor>(readerCloud_.getTensorFieldNames());
    internalAddFields<symmTensor>(readerCloud_.getSymmTensorFieldNames());
    internalAddFields<sphericalTensor>(readerCloud_.getSphericalTensorFieldNames());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CloudProxyForReaderParticle::~CloudProxyForReaderParticle()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<scalar> > CloudProxyForReaderParticle::getScalarField(
    const word &name
) const
{
    if(readerCloud_.hasScalar(name)) {
        return tmp<Field<scalar> >(
            new Field<scalar>(readerCloud_.getScalarField(name))
        );
    }
    if(readerCloud_.hasLabel(name)) {
        Field<label> orig(readerCloud_.getLabelField(name));
        tmp<Field<scalar> > result(new Field<scalar>(orig.size()));
        forAll(orig,i)
        {
            const_cast<scalar&>(result()[i])=orig[i];
        }
        return result;
    }

    return CloudProxyForParticle<ReaderParticleCloud>::getScalarField(name);
}

tmp<Field<vector> > CloudProxyForReaderParticle::getVectorField(
    const word &name
) const
{
    if(readerCloud_.hasVector(name)) {
        return tmp<Field<vector> >(
            new Field<vector>(readerCloud_.getVectorField(name))
        );
    }

    return CloudProxyForParticle<ReaderParticleCloud>::getVectorField(name);
}

tmp<Field<tensor> > CloudProxyForReaderParticle::getTensorField(
    const word &name
) const
{
    if(readerCloud_.hasTensor(name)) {
        return tmp<Field<tensor> >(
            new Field<tensor>(readerCloud_.getTensorField(name))
        );
    }

    return CloudProxyForParticle<ReaderParticleCloud>::getTensorField(name);
}

tmp<Field<symmTensor> > CloudProxyForReaderParticle::getSymmTensorField(
    const word &name
) const
{
    if(readerCloud_.hasSymmTensor(name)) {
        return tmp<Field<symmTensor> >(
            new Field<symmTensor>(readerCloud_.getSymmTensorField(name))
        );
    }

    return CloudProxyForParticle<ReaderParticleCloud>::getSymmTensorField(name);
}

tmp<Field<sphericalTensor> > CloudProxyForReaderParticle::getSphericalTensorField(
    const word &name
) const
{
    if(readerCloud_.hasSphericalTensor(name)) {
        return tmp<Field<sphericalTensor> >(
            new Field<sphericalTensor>(readerCloud_.getSphericalTensorField(name))
        );
    }

    return CloudProxyForParticle<ReaderParticleCloud>::getSphericalTensorField(name);
}

CloudProxy::addcloudConstructorToTable< CloudProxyForReaderParticle > add_lookup_CloudProxyForReaderParticlecloudConstructorToCloudProxyTable_("Cloud<ReaderParticle>");

} // namespace end

// ************************************************************************* //
