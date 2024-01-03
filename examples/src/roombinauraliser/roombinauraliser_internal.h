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
 * @file: roombinauraliser_internal.h
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain.
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @author Leo McCormack
 * @date 25.09.2017
 * @license ISC
 */

#ifndef __roombinauraliser_INTERNAL_H_INCLUDED__
#define __roombinauraliser_INTERNAL_H_INCLUDED__

#include "roombinauraliser.h"  /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */
#include "netcdf.h"        /* Include NetCDF*/
#include "../../../framework/resources/afSTFT/afSTFTlib.h"
//#include " ../../../modules/saf_sofa_reader/libmysofa/internal/mysofa_internal.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(roombinauraliser_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define roombinauraliser_FRAME_SIZE ( FRAME_SIZE )          /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define roombinauraliser_FRAME_SIZE ( 128 )                 /**< Framesize, in time-domain samples */
# endif
#endif
#define HOP_SIZE ( 128 )                                  /**< STFT hop size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                     /**< Number of frequency bands */
#define TIME_SLOTS ( roombinauraliser_FRAME_SIZE / HOP_SIZE ) /**< Number of STFT timeslots */

/* Checks: */
#if (roombinauraliser_FRAME_SIZE % HOP_SIZE != 0)
# error "roombinauraliser_FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for roombinauraliser. Contains variables for audio buffers,
 * afSTFT, HRTFs, internal variables, flags, user parameters.
 * Note: if this is modified, identically modify _roombinauraliserNF struct.
 */
typedef struct _roombinauraliser
{
    /* audio buffers */
    float** inputFrameTD;            /**< time-domain input frame; #MAX_NUM_INPUTS x #roombinauraliser_FRAME_SIZE */
    float** outframeTD;              /**< time-domain output frame; #NUM_EARS x #roombinauraliser_FRAME_SIZE */
    float_complex*** inputframeTF;   /**< time-frequency domain input frame; #HYBRID_BANDS x #MAX_NUM_INPUTS x #TIME_SLOTS */
    float_complex*** outputframeTF;  /**< time-frequency domain output frame; #HYBRID_BANDS x #NUM_EARS x #TIME_SLOTS */
    int fs;                          /**< Host sampling rate, in Hz */
    float freqVector[HYBRID_BANDS];  /**< Frequency vector (filterbank centre frequencies) */
    void* hSTFT;                     /**< afSTFT handle */
    
    /* sofa file info */
    char* sofa_filepath;             /**< absolute/relevative file path for a sofa file */
    float** hrirs;                   /**< time domain HRIRs; FLAT: N_hrir_dirs x #NUM_EARS x hrir_len */
    float* hrir_dirs_deg;            /**< directions of the HRIRs in degrees [azi elev]; FLAT: N_hrir_dirs x 2 */
    int N_hrir_dirs;                 /**< number of HRIR directions in the current sofa file */
    int hrir_loaded_len;             /**< length of the loaded HRIRs, in samples */
    int hrir_runtime_len;            /**< length of the HRIRs being used for processing (after any resampling), in samples */
    int hrir_loaded_fs;              /**< sampling rate of the loaded HRIRs  */
    int hrir_runtime_fs;             /**< sampling rate of the HRIRs being used for processing (after any resampling) */
    float* weights;                  /**< Integration weights for the HRIR measurement grid */

    /* Diffuse field eq data */
    float* fabian_cir;
    float_complex* ctf_fb;         /**< hrtf filterbank coefficients; nBands x nCH x N_hrirs */

    int N_samples_fabian_cir;

    /* vbap gain table */
    int hrtf_vbapTableRes[2];        /**< [0] azimuth, and [1] elevation grid resolution, in degrees */
    int N_hrtf_vbap_gtable;          /**< Number of interpolation weights/directions */
    int* hrtf_vbap_gtableIdx;        /**< N_hrtf_vbap_gtable x 3 */
    float* hrtf_vbap_gtableComp;     /**< N_hrtf_vbap_gtable x 3 */
    
    /* hrir filterbank coefficients */
    float** itds_s;                  /**< interaural-time differences for each HRIR (in seconds); nSources x nBands x 1 */
    double_complex** hrtf_fb;         /**< hrtf filterbank coefficients; nBands x nCH x N_hrirs */
    float** hrtf_fb_mag;             /**< magnitudes of the hrtf filterbank coefficients; nBands x nCH x N_hrirs */
    float_complex hrtf_interp[MAX_NUM_INPUTS][HYBRID_BANDS][NUM_EARS]; /**< Interpolated HRTFs */
    
    /* flags/status */
    CODEC_STATUS codecStatus;               /**< see #CODEC_STATUS */
    float progressBar0_1;                   /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;                  /**< Current (re)initialisation step, string */
    char* progressBarTooltip;               /**< Tooltip for current (re)initialisation step, string */
    PROC_STATUS procStatus;                 /**< see #PROC_STATUS */
    int recalc_hrtf_interpFLAG;             /**< 1: re-calculate/interpolate the HRTF, 0: do not */
    REINIT_MODES reInitHRTFsAndGainTables;  /**< 1: reinitialise the HRTFs and interpolation tables, 2: only resample HRTFs, 3: only apply diffuse field EQ, 0: do nothing */
    int recalc_M_rotFLAG;                   /**< 1: re-calculate the rotation matrix, 0: do not */
    int VBAP_3d_FLAG;                       /**< 1: VBAP in 3 Dimensions, 0: VBAP in 2 Dimensions */

    /* misc. */
    float rot_deg[2];                       /**< Intermediate rotated reference frame, in degrees */
    float rot_xyz[3];                       /**< Intermediate rotated reference frame, as unit-length Cartesian coordinates */
    float src_dirs_xyz[MAX_NUM_INPUTS][3];  /**< Intermediate source directions, as unit-length Cartesian coordinates  */
    int nTriangles;                         /**< Number of triangles in the convex hull of the spherical arrangement of HRIR directions/points */

    /* user parameters */
    int nSources;                           /**< Current number of input/source signals */
    float src_dirs_deg[MAX_NUM_INPUTS][2];  /**< Current source/panning directions, in degrees */
    INTERP_MODES interpMode;                /**< see #INTERP_MODES */
    DIFF_EQ_MODES diffEqMode;
    int useDefaultHRIRsFLAG;                /**< 1: use default HRIRs in database, 0: use those from SOFA file */
    int enableBRIRsDiffuseEQ;               /**< flag to diffuse-field equalisation to the currently loaded BRIRs */
    int enableRotation;                     /**< 1: enable rotation, 0: disable */
    int enablePartConv;                     /**< 1 enable partitioned convolution, 0: disable */
    float yaw;                              /**< yaw (Euler) rotation angle, in degrees */
    float roll;                             /**< roll (Euler) rotation angle, in degrees */
    float pitch;                            /**< pitch (Euler) rotation angle, in degrees */
    int bFlipYaw;                           /**< flag to flip the sign of the yaw rotation angle */
    int bFlipPitch;                         /**< flag to flip the sign of the pitch rotation angle */
    int bFlipRoll;                          /**< flag to flip the sign of the roll rotation angle */
    int useRollPitchYawFlag;                /**< rotation order flag, 1: r-p-y, 0: y-p-r */
    float src_gains[MAX_NUM_INPUTS];        /**< Gains applied per source */

} roombinauraliser_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void roombinauraliser_setCodecStatus(void* const hBin,
                                 CODEC_STATUS newStatus);

/**
 * Interpolates between (up to) 3 HRTFs via amplitude-normalised VBAP gains.
 *
 * The HRTF magnitude responses and HRIR ITDs are interpolated seperately before
 * re-introducing the phase.
 *
 * @param[in]  hBin          roombinauraliser handle
 * @param[in]  mode          see #INTERP_MODES 
 * @param[in]  azimuth_deg   Source azimuth in DEGREES
 * @param[in]  elevation_deg Source elevation in DEGREES
 * @param[out] h_intrp       Interpolated HRTF
 */
void roombinauraliser_interpHRTFs(void* const hBin,
                              INTERP_MODES mode,
                              float azimuth_deg,
                              float elevation_deg,
                              float_complex h_intrp[MAX_NUM_INPUTS][HYBRID_BANDS][NUM_EARS]);

/**
 * Initialise the HRTFs: either loading the default set or loading from a SOFA
 * file; and then generate a VBAP gain table for interpolation.
 *
 * @note Call roombinauraliser_initTFT() (if needed) before calling this function
 */
void roombinauraliser_initHRTFsAndGainTables(void* const hBin);

/**
 * Initialise the filterbank used by roombinauraliser.
 *
 * @note Call this function before roombinauraliser_initHRTFsAndGainTables()
 */
void roombinauraliser_initTFT(void* const hBin);

static double fabian_ir = {6.347623757418916490e-01,-3.229786055630429198e-01,1.802470897633837166e-01,-4.306831691412307123e-02,2.021296727078609190e-01,-1.278432979319796617e-03,9.211553921811316270e-02,1.112025697219754666e-02,7.102338822541492958e-02,2.116948839261354430e-03,3.533534474459527897e-02,-4.782989278813806755e-03,1.961155907668466619e-02,-8.324748182822951093e-03,1.890523538003034687e-02,-5.151430063348769184e-03,1.720650523149344635e-02,1.534134679661368350e-03,1.382324977924552500e-02,1.913833833957444969e-03,1.158019821704660106e-02,5.637062863037708982e-03,1.140537984280109302e-02,3.263822223500500649e-03,8.795929652983784305e-03,1.387702335816867142e-03,7.939434305533257102e-03,2.960683105060799665e-04,6.152246061429689564e-03,-1.074279428851980598e-04,4.396600362471919307e-03,1.461175100134327691e-04,3.282734279044176480e-03,5.852363377029396602e-04,3.019699471470465788e-03,6.767663603708469651e-04,3.362895951400109596e-03,6.903585806225464271e-04,3.830261417749242919e-03,7.804607314488692113e-04,4.022008517589376400e-03,1.188374153931281638e-03,3.886273952225215329e-03,1.763702340519725286e-03,3.659871236683768375e-03,2.204471970341330238e-03,3.490221418746291063e-03,2.231831089350578244e-03,3.459951380202113021e-03,1.921618628821775837e-03,3.385438356845829090e-03,1.435524263665358949e-03,3.069915761480931145e-03,1.020093874171674846e-03,2.491780464480404522e-03,7.384330485773573424e-04,1.780032758710380515e-03,5.001776305672350157e-04,1.181742902257324791e-03,2.089139677159372198e-04,8.342981001236377068e-04,-1.436960705956030543e-04,7.332742609899506932e-04,-4.467872706527811419e-04,7.045842202551601003e-04,-5.819032337736750483e-04,6.271958237695415004e-04,-5.311760281885025961e-04,4.657348183797568484e-04,-3.975940147125833941e-04,2.999485271200753671e-04,-3.226001460338417384e-04,2.031517603165563168e-04,-4.043234225591980794e-04,1.937316056310017212e-04,-6.016456786100988227e-04,1.712022696134814616e-04,-8.126918653247878870e-04,5.642620238534838033e-05,-9.491824799252341228e-04,-1.794864082295016035e-04,-1.004400902821920433e-03,-4.529972035075383084e-04,-1.054187599141759909e-03,-6.583482570512188582e-04,-1.149442939304545254e-03,-7.305723239770893756e-04,-1.291074970891296306e-03,-6.998708779752012277e-04,-1.398947012674592707e-03,-6.441176647479178842e-04,-1.401926412700444735e-03,-6.296651772440758172e-04,-1.279979034212476484e-03,-6.440878022616795080e-04,-1.105199989069475261e-03,-6.326425617335096341e-04,-9.625115130869978976e-04,-5.436314907687025185e-04,-9.006902624184477549e-04,-3.925641247500917443e-04,-8.932135272764027845e-04,-2.433453056069400993e-04,-8.688018558849731291e-04,-1.677868946779114366e-04,-7.866683299209764136e-04,-1.704707425258825083e-04,-6.564725226885076686e-04,-1.993057750803129097e-04,-5.436230146743372890e-04,-1.821098818488891462e-04,-4.932153731812146356e-04,-8.811930370803075789e-05,-4.982755378507687603e-04,4.047903540198115454e-05,-5.087682612819327673e-04,1.367125029648291431e-04,-4.608065974128649654e-04,1.659864469949497005e-04,-3.472715591582585569e-04,1.572848156935205345e-04,-2.113831248794097834e-04,1.804512530579007357e-04,-1.090367411148658242e-04,2.806144489495618734e-04,-9.502422741275887521e-05,4.297870178485137336e-04,-1.934093319360942161e-04,3.295705182845652258e-04,-1.011353792534521563e-04,9.054033130876370519e-05,-3.264446062156329973e-05,5.004780685438044499e-05,-2.827877429436320044e-05,1.018938067101397614e-05,-3.748103605267806413e-05,-7.563044336005482154e-06,-3.955124109493061502e-05,-1.398747623301070225e-05,-3.543049052078638134e-05,-1.030693463865157080e-05,-2.779770031220376745e-05,-5.537350974872968198e-06,-2.013198950549518145e-05,-3.731379920921948326e-06,-1.549381312407833230e-05,-3.401703883545242497e-06,-1.303655797884452877e-05,-3.378006420237262070e-06,-1.178074444674698332e-05,-2.873806752926856326e-06,-1.096920838037915763e-05,-1.997708322409193002e-06,-9.551962520890339809e-06,-7.752050529048797014e-07,-7.491975695235659856e-06,3.070530341089298475e-07,-5.174365714767041586e-06,1.402136940959674404e-06,-2.910591281285936466e-06,2.711339979440113721e-06,-9.658212816898394514e-07,4.243900570135445956e-06,6.322999923591176586e-07,5.791084468178990245e-06,2.036470383781585818e-06,7.097078122154883612e-06,3.371567284493380548e-06,7.995602272418335946e-06,4.608377486226298897e-06,8.501099102866963779e-06,5.562949740222741015e-06,8.723820820193048918e-06,6.061910819202179885e-06,8.764889910813323327e-06,6.082213854643997076e-06,8.638659299445634022e-06,5.773608093803410439e-06,8.304193645789331973e-06,5.348522183495517043e-06,7.753198727925229213e-06,4.956483907377884335e-06,7.065953604858477663e-06,4.619603460185547701e-06,6.384686530622318592e-06,4.281295980858042064e-06,5.828461921494376299e-06,3.897186108047870638e-06,5.418819140904903168e-06,3.482799162167997189e-06,5.076829754286610089e-06,3.095952656865587077e-06,4.697511895659460863e-06,2.770437798516925444e-06,4.229279227976893718e-06,2.481095274344526692e-06,3.696067641788190734e-06,2.159476721941442368e-06,3.161283659297449417e-06,1.751146222286991929e-06,2.668529513605654564e-06,1.263349163993694559e-06,2.209349018888094886e-06,7.628735285650453436e-07,1.750232237157596976e-06,3.238855673936957624e-07,1.278573749216673027e-06,-1.676131758171181321e-08,8.267254411261153952e-07,-2.746158187333242040e-07,4.547930803729659316e-07,-4.837585802628400635e-07,2.041634794770004007e-07,-6.597905759872652094e-07,6.949969103445809214e-08,-7.863818491988552931e-07,4.663265469659930386e-09,-8.368478079453864608e-07,-3.919615703916220715e-08,-8.071175481201342491e-07,-8.178677795202159725e-08,-7.293683696737789078e-07,-1.106800165369704630e-07,-6.513132202943487213e-07,-1.078565918042252748e-07,-6.034958604758106348e-07,-7.444622294744570050e-08,-5.817858854067767986e-07,-3.486836780761190102e-08,-5.580810271100726710e-07,-1.567323052567202658e-08,-5.089673824829190696e-07,-2.121177819223922314e-08,-4.357118129847898169e-07,-3.146131648863904114e-08,-3.597972810950774524e-07,-2.085375840036693904e-08,-3.012759447756772103e-07,1.810938964710525451e-08,-2.608179656399772839e-07,6.712920498991806555e-08,-2.225000084547728792e-07,1.015123886347437471e-07,-1.685384311074530895e-07,1.120519164072019998e-07,-9.806510439177268810e-08,1.112391649570293542e-07,-3.039341863756408877e-08,1.162228222096052485e-07,1.226189978661771360e-08,1.272632411674036148e-07,1.991670334215750343e-08,1.242260081137231127e-07,-2.947287521939403714e-09,7.869989944466298288e-08,-2.011494816388753958e-08};

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __roombinauraliser_INTERNAL_H_INCLUDED__ */
