template<class ParcelType>
template<class TrackCloudType>
inline Foam::PlateletParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud
)
:
    ParcelType::trackingData(cloud),
    tauInterp_
    (
        interpolation<symmTensor>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.tau()
        )
    )
{}