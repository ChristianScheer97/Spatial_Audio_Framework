/*
 * Copyright 2017-2018 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/**
 *@addtogroup SOFA_Reader
 *@{
 * @file saf_sofa_reader.h 
 * @brief Main header for the sofa reader module (#SAF_SOFA_READER_MODULE)
 *
 * @note This SOFA reader may optionally use netcdf if "SAF_ENABLE_NETCDF" is
 *       defined. Otherwise, the reader will use libmysofa [1] in conjunction
 *       with zlib, which is included in framework/resources/zlib; i.e. no
 *       external libraries need to be linked by default.
 *
 * @see [1] https://github.com/hoene/libmysofa (BSD-3-Clause license)
 *
 * @author Leo McCormack
 * @date 21.11.2017
 * @license ISC
 */

#ifndef __SAF_SOFA_READER_H_INCLUDED__
#define __SAF_SOFA_READER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef SAF_ENABLE_SOFA_READER_MODULE

#include "libmysofa/mysofa.h"

#ifdef SAF_ENABLE_NETCDF
#include "netcdf.h"
#endif /* SAF_ENABLE_NETCDF */

/** SOFA file reader options */
typedef enum{
    /** The default option is #SAF_SOFA_READER_OPTION_LIBMYSOFA */
    SAF_SOFA_READER_OPTION_DEFAULT,

    /** This option uses the libmysofa library to load SOFA files, which is
     *  adopted from: https://github.com/hoene/libmysofa (BSD-3-Clause license)
     *
     *  The benefits of this option is that it only depends on zlib, which is
     *  included in SAF. While the downsides of this option, is that zlib has
     *  file size limits for each chunk (<4GB) and it is quite slow at
     *  decompressing large files. */
    SAF_SOFA_READER_OPTION_LIBMYSOFA,

    /** If SAF_ENABLE_NETCDF is defined, then an alternative SOFA reader may be
     *  used. This version requires netcdf to be linked to SAF, along with its
     *  dependencies. The netcdf loader gets around the file size limits of
     *  the libmysofa loader and is also approximately 3 times faster.
     *  Therefore, if you intend to load many large SOFA files
     *  (especially microphone arrays or Ambisonic IRs), then this alternative
     *  SOFA reader is either required (to get around the file size limit) or
     *  may be preferred due to the shorter loading times. The downsides of
     *  using the netcdf option is that it is NOT thread-safe! and requires
     *  these additional external libraries to be linked to SAF. */
    SAF_SOFA_READER_OPTION_NETCDF

} SAF_SOFA_READER_OPTIONS;

/** SOFA file use cases (HRTF, BRIR, SRIR etc.) */
typedef enum {
    /** The default use case is HRTF */
    SAF_SOFA_READER_USECASE_DEFAULT,

    SAF_SOFA_READER_USECASE_HRIR,

    SAF_SOFA_READER_USECASE_BRIR
} SAF_SOFA_READER_USECASE;


/* ========================================================================== */
/*                          Public Structures/Enums                           */
/* ========================================================================== */

/**
 * SOFA container struct comprising all possible data that can be extracted
 * from SOFA 1.0 and 2.1 files; as laid down in the GeneralFIR, SimpleFreeFieldHRIR, 
 * MultiSpeakerBRIR and SingleRoomMIMOSRIR specifications:
 *    https://www.sofaconventions.org/mediawiki/index.php/GeneralFIR
 *    https://www.sofaconventions.org/mediawiki/index.php/SimpleFreeFieldHRIR
 *    https://www.sofaconventions.org/mediawiki/index.php/MultiSpeakerBRIR
 *    https://www.sofaconventions.org/mediawiki/index.php/SingleRoomMIMOSRIR
 */
typedef struct _saf_sofa_container{
    /* All possible SOFA variables (defaults={-1|NULL}) */
    int nSources;                 /**< Number of source/measurement positions */
    int nReceivers;               /**< Number of ears/number of mics etc. */
    int DataLengthIR;             /**< Length of the IRs, in samples */
    float* DataIR;                /**< The impulse response (IR) Data;
                                   *   FLAT:nSources x nReceivers x DataLengthIR*/
    float DataSamplingRate;       /**< Sampling rate used to measure the IRs */
    float* DataDelay;             /**< Delay in samples; nReceivers x 1 */
    float* SourcePosition;        /**< Source positions (refer to
                                   *   SourcePositionType & SourcePositionUnits
                                   *   for the convention and units);
                                   *   FLAT: nSources x 3 */
    float* ReceiverPosition;      /**< Receiver positions (refer to
                                   *   ReceiverPositionType &
                                   *   ReceiverPositionUnits for the convention
                                   *   and units);
                                   *   FLAT: nReceivers x 3 */
    int nListeners;               /**< Number of listener positions */
    int nEmitters;                /**< Number of emitter positions */
    float* ListenerPosition;      /**< Listener position (The object
                                   *   incorporating all receivers; refer to
                                   *   ListenerPositionType &
                                   *   ListenerPositionUnits for the convention
                                   *   and units); FLAT: nListeners x 3  */
    float* ListenerUp;            /**< Vector pointing upwards from the listener
                                   *   position (Cartesian); 1 x 3 or FLAT: nListeners x 3  */
    float* ListenerView;          /**< Vector pointing forwards from the
                                   *   listener position (Cartesian); 3 x 1 */
    float* EmitterPosition;       /**< Positions of acoustic excitation used for
                                   *   the measurement (refer to
                                   *   EmitterPositionType &
                                   *   EmitterPositionUnits for the convention
                                   *   and units); FLAT: nEmitters x 3  or nEmitters x 3 x mMeasurements*/
    float* EmitterUp;             /**< Vector pointing upwards from the emitter
                                   *   position (Cartesian); [E C I] or [E C M]; nEmitters x3 or nEmitter x 3 x mMeasurements */ 
    float* EmitterView;           /**< Vector pointing forwards from the
                                   *   emitter position (Cartesian) [E C I] or [E C M]; nEmitters x3 or nEmitter x 3 x mMeasurements  */
    float* RoomTemperature;       /**< Temperature during measurements, given in Kelvin (if not differently set in RoomTemperatureUnits). [I] or [M]
                                   *   1 x 1 or FLAT: mMeasurements x 1 */
    float* RoomVolume;            /**< Volume of the room [I] or [M I]
                                   *   1 x 1 or FLAT: mMeasurements x 1 */
    float* RoomCornerA;           /**< Cartesian yxordinate of edge A [I C] or [M C]
                                   *   1 x 3 or FLAT: mMeasurements x 3*/
    float* RoomCornerB;           /**< Cartesian yxordinate of edge B [I C] or [M C]
                                   *   1 x 3 or FLAT: mMeasurements x 3*/
    int RoomCorners;              /**< The value of this attribute is to be ignored. It only exist to for RoomCorners:Type and RoomCorners:Units  [II]*/
    float* ReceiverView;          /**< View vector for the orientation. [R C I] or [R C M]
                                   *   rReceivers x 3  or  rReceivers x 3 x rMeasurements */
    float* ReceiverUp;            /**< Up vector for the orientation. [R C I] or [R C M] 
                                   *   rReceivers x 3 x 1 or rReceivers x 3 x mMeasurements */
    float* SourceView;            /**< Vector pointing forwards from the
                                   *   source position (Cartesian) [I C] or [M C]; 3 x 1 or 3  */
    float* SourceUp;              /**< /**< Vector pointing upwards from the listener [I C] or [M C]
                                   *   position (Cartesian); 1 x 3 or FLAT: mMeasurements x 3  */ 
    
    int* MeasurementDate;         /**< Optional M-dependent date and time of the measurement. 
                                   *   8 or 8 * mMeasurements */


    /* All possible SOFA strings (defaults=NULL) */
    char* EmitterDescriptions;    /**< E-dependent version of the attribute EmitterDescription. [E S] or [E S M]
                                   *   When more than one string is considered in a file, S shall represent the size of the character array with the longest string dimension.
                                   *   Entries shorter than S shall be padded with null characters up to the length of S.
                                       eEmitters x sCharacters or eEmitters x sCharacters x mMEasurements */
    char* ReceiverDescriptions;   /**< R-dependent version of the attribute ReceiverDescription. [R S] or [R S M]
                                   *   When more than one string is considered in a file, S shall represent the size of the character array with the longest string dimension.
                                   *   Entries shorter than S shall be padded with null characters up to the length of S.
                                   *   rReceivers x sCharacters or sReceivers x sCharacters x mMeasurements */
    
                                   
                                   /* All possible SOFA variable attributes (defaults=NULL) */
    char* ListenerPositionType;   /**< {'cartesian'|'spherical'} */
    char* ListenerPositionUnits;  /**< {'degree, degree, metre'|'metre'} */
    char* ListenerViewType;       /**< {'cartesian'|'spherical'} */
    char* ListenerViewUnits;      /**< {'degree, degree, metre'|'metre'} */
    char* ReceiverPositionType;   /**< {'cartesian'|'spherical'} */
    char* ReceiverPositionUnits;  /**< {'degree, degree, metre'|'metre'} */
    char* ReceriverViewType;      /**< {'cartesian'} */
    char* ReceriverViewUnits;     /**< {'metre'} */
    char* RoomCornersType;        /**< {'cartesian'} */
    char* RoomCornersUnits;       /**< {'metre'} */
    char* RoomTemperaturUnits;    /**< {'kelvin' | 'celsius' | 'fahrenheit'} */
    char* RoomVolumeUnits;        /**< {'cubic metre' */
    char* SourcePositionType;     /**< {'cartesian'|'spherical'} */
    char* SourcePositionUnits;    /**< {'degree, degree, metre'|'metre'} */
    char* SourceViewType;         /**< {'cartesian'} */
    char* SourceViewUnits;        /**< {'metre'} */
    char* EmitterPositionType;    /**< {'cartesian'|'spherical'} */
    char* EmitterPositionUnits;   /**< {'degree, degree, metre'|'metre'} */
    char* EmitterViewType;        /**< {'cartesian'} */
    char* EmitterViewUnits;        /**< {'metre'} */
    char* DataSamplingRateUnits;  /**< {'hertz'} */

    /* All possible SOFA global attributes (defaults=NULL) */
    char* Conventions;            /**< {'SOFA'} */
    char* Version;                /**< Version number */
    char* SOFAConventions;        /**< {'GeneralFIR'|
                                        'GeneralTF'|
                                   *    'SimpleFreeFieldHRIR'|
                                        'SingleRoomMIMOSRIR'|
                                        'MultiSpeakerBRIR'} */
    char* SOFAConventionsVersion; /**< SOFA convention number */
    char* APIName;                /**< API name */
    char* APIVersion;             /**< API version */
    char* ApplicationName;        /**< Name of Application that created file */
    char* ApplicationVersion;     /**< Ver. of Application that created file */
    char* AuthorContact;          /**< Contact information */
    char* Comment;                /**< File comments */
    char* DataType;               /**< {'FIR'|'TF'|'FIR-E'} */
    char* History;                /**< History information */
    char* License;                /**< License under which file is provided */
    char* Organisation;           /**< Organisation reponsible for the file */
    char* References;             /**< References */
    char* RoomDescription;        /**< Informal verbal description of the room */
    char* RoomGeometry;           /**< URI to a file describing the room geometry */
    char* RoomLocation;           /**< Location of the room */
    char* RoomShortName;          /**< Short name of the Room */
    char* RoomType;               /**< Room type (free field, shoebox, dae etc.) */
    char* Origin;                 /**< Where this file came from */
    char* Organization;           /**< Organization that provided the data */
    char* DateCreated;            /**< Date file was created */
    char* DateModified;           /**< Date file was modified */
    char* Title;                  /**< Title of file */
    char* DatabaseName;           /**< Name of the database. Used for classification of the data. */
    char* ListenerShortName;      /**< Name of the listener/dummyhead/mic etc.*/
    char* ListenerDescription;    /**< Description of the listener */
    char* ReceiverShortName;      /**< Short name of the receiver */
    char* ReceiverDescription;    /**< Description of the receiver */
    char* SourceShortName;        /**< Short name of the source */
    char* SourceDescription;      /**< Description of the source */
    char* EmitterShortName;       /**< Short name of the emitter */
    char* EmitterDescription;     /**< Description of the emitter */

    /* libmysofa handle, which is used if SAF_ENABLE_NETCDF is not defined */
    void* hLMSOFA;                /**< libmysofa handle */

}saf_sofa_container;

/** SOFA loader error codes */
typedef enum{
    /** None of the error checks failed */
    SAF_SOFA_OK,
    /** Not a SOFA file, or no such file was found in the specified location */
    SAF_SOFA_ERROR_INVALID_FILE_OR_FILE_PATH,
    /** Dimensions of the SOFA data were not as expected */
    SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED,
    /** The data-type of the SOFA data was not as expected */
    SAF_SOFA_ERROR_FORMAT_UNEXPECTED,
    /** NetCDF is not thread safe! */
    SAF_SOFA_ERROR_NETCDF_IN_USE,
    /** The wrong reader option was chosen. This error accurs, 
    if you try to load a BRIR with the SAF_SOFA_READER_USECASE_BRIR flag 
    without the reader option SAF_SOFA_READER_OPTION_NETCDF*/
    SAF_SOFA_ERROR_INVALID_READER_OPTION

} SAF_SOFA_ERROR_CODES;


/* ========================================================================== */
/*                              Main Functions                                */
/* ========================================================================== */

/**
 * Fills a 'sofa_container' with data found in a SOFA file (GeneralFIR,
 * SimpleFreeFieldHRIR, SingleRoomMIMOSRIR or MultiSpeakerBRIR), as detailed 
 * in the SOFA 1.0 and 2.1 standard [1,2,3,4,5]
 *
 * @warning This loader currently does not support TF SOFA files!
 * @note If you encounter a SOFA file that this SOFA loader cannot load, (or it
 *       misses some of the data) then please send it to the developers :-)
 * @test test__saf_sofa_open(), test__mysofa_load(), test__sofa_comparison()
 *
 * @param[in] hSOFA         The sofa_container
 * @param[in] sofa_filepath SOFA file path (including .sofa extension)
 * @param[in] option        See #SAF_SOFA_READER_OPTIONS
 * @returns An error code (see #SAF_SOFA_ERROR_CODES)
 *
 * @see [1] Majdak, P., Iwaya, Y., Carpentier, T., Nicol, R., Parmentier, M.,
 *          Roginska, A., Suzuki, Y., Watanabe, K., Wierstorf, H., Ziegelwanger,
 *          H. and Noisternig, M., 2013, May. Spatially oriented format for
 *          acoustics: A data exchange format representing head-related transfer
 *          functions. In Audio Engineering Society Convention 134. Audio
 *          Engineering Society.
 * @see [2] https://www.sofaconventions.org/mediawiki/index.php/GeneralFIR
 * @see [3] https://www.sofaconventions.org/mediawiki/index.php/SimpleFreeFieldHRIR
 * @see [4] https://www.sofaconventions.org/mediawiki/index.php/SingleRoomMIMOSRIR
 * @see [5] https://www.sofaconventions.org/mediawiki/index.php/MultiSpeakerBRIR
 */
SAF_SOFA_ERROR_CODES saf_sofa_open_universal(saf_sofa_container* hSOFA,
                                   char* sofa_filepath,
                                   SAF_SOFA_READER_OPTIONS option,
                                    SAF_SOFA_READER_USECASE usecase);

SAF_SOFA_ERROR_CODES saf_sofa_open(saf_sofa_container* hSOFA,
                                   char* sofa_filepath,
                                   SAF_SOFA_READER_OPTIONS option);

/**
 * Frees all SOFA data in a sofa_container
 *
 * @param[in] hSOFA The sofa_container
 */
void saf_sofa_close(saf_sofa_container* hSOFA);

#endif /* SAF_ENABLE_SOFA_READER_MODULE */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_SOFA_READER_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup SOFA_Reader */
