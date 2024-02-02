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
#include "netcdf.h"        /* Include NetCDF */
#include "../../../../Partitioned-Convolution/lib/Convolution.h"   /* Include Partitioned Convolution */
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
#define LATENCY 256

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
    
    /* partitioned convolution */
    
    Convolution* hPart_current_left;
    Convolution* hPart_current_right;
    Convolution* hPart_new_left;
    Convolution* hPart_new_right;    /**< partitioned convolution handles (2x2) for each ear and to be able to seamlessly crossfade between new and old audio */
    
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
    float_complex* ctf_fb;           /**< ctf filterbank coefficients; nBands */
    int N_samples_fabian_cir;

    /* vbap gain table */
    int hrtf_vbapTableRes[2];        /**< [0] azimuth, and [1] elevation grid resolution, in degrees */
    int N_hrtf_vbap_gtable;          /**< Number of interpolation weights/directions */
    int* hrtf_vbap_gtableIdx;        /**< N_hrtf_vbap_gtable x 3 */
    float* hrtf_vbap_gtableComp;     /**< N_hrtf_vbap_gtable x 3 */
    
    /* hrir filterbank coefficients */
    float** itds_s;                  /**< interaural-time differences for each HRIR (in seconds); nSources x nBands x 1 */
    float_complex** hrtf_fb;         /**< hrtf filterbank coefficients; nBands x nCH x N_hrirs */
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
    DIFF_EQ_MODES diffEqMode;               /**< see #DIFF_EQ_MODES */
    int useDefaultHRIRsFLAG;                /**< 1: use default HRIRs in database, 0: use those from SOFA file */
    int enableHRIRsDiffuseEQ;               /**< flag to diffuse-field equalisation to the currently loaded HRTFs */
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
 * outputs the processed audio from partitioned convolution into output buffer
 *
 *
 */
void roombinauraliser_convoutput(void* SELF, dft_sample_t* output, int num_samples, void* const hBin);


/**
 * Initialise the filterbank used by roombinauraliser.
 *
 * @note Call this function before roombinauraliser_initHRTFsAndGainTables()
 */
void roombinauraliser_initTFT(void* const hBin);
static float fabian_ir[] = {6.347623467445373535e-01,-3.229786157608032227e-01,1.802470833063125610e-01,-4.306831583380699158e-02,2.021296769380569458e-01,-1.278433017432689667e-03,9.211553633213043213e-02,1.112025696784257889e-02,7.102338969707489014e-02,2.116948831826448441e-03,3.533534333109855652e-02,-4.782989155501127243e-03,1.961155980825424194e-02,-8.324747905135154724e-03,1.890523545444011688e-02,-5.151430144906044006e-03,1.720650494098663330e-02,1.534134731628000736e-03,1.382324937731027603e-02,1.913833781145513058e-03,1.158019807189702988e-02,5.637062713503837585e-03,1.140537951141595840e-02,3.263822291046380997e-03,8.795930072665214539e-03,1.387702301144599915e-03,7.939434610307216644e-03,2.960683195851743221e-04,6.152246147394180298e-03,-1.074279425665736198e-04,4.396600183099508286e-03,1.461175124859437346e-04,3.282734192907810211e-03,5.852363537997007370e-04,3.019699361175298691e-03,6.767663871869444847e-04,3.362895920872688293e-03,6.903585745021700859e-04,3.830261528491973877e-03,7.804607157595455647e-04,4.022008739411830902e-03,1.188374124467372894e-03,3.886274062097072601e-03,1.763702370226383209e-03,3.659871174022555351e-03,2.204471966251730919e-03,3.490221453830599785e-03,2.231831196695566177e-03,3.459951374679803848e-03,1.921618590131402016e-03,3.385438350960612297e-03,1.435524318367242813e-03,3.069915808737277985e-03,1.020093914121389389e-03,2.491780556738376617e-03,7.384330383501946926e-04,1.780032762326300144e-03,5.001776153221726418e-04,1.181742874905467033e-03,2.089139743475243449e-04,8.342980872839689255e-04,-1.436960737919434905e-04,7.332742679864168167e-04,-4.467872786335647106e-04,7.045842357911169529e-04,-5.819032085128128529e-04,6.271958118304610252e-04,-5.311760469339787960e-04,4.657348326873034239e-04,-3.975940053351223469e-04,2.999485295731574297e-04,-3.226001572329550982e-04,2.031517651630565524e-04,-4.043234221171587706e-04,1.937316119438037276e-04,-6.016456754878163338e-04,1.712022640276700258e-04,-8.126918692141771317e-04,5.642620089929550886e-05,-9.491824894212186337e-04,-1.794864074327051640e-04,-1.004400895908474922e-03,-4.529971920419484377e-04,-1.054187654517591000e-03,-6.583482609130442142e-04,-1.149442978203296661e-03,-7.305723265744745731e-04,-1.291075022891163826e-03,-6.998708704486489296e-04,-1.398946973495185375e-03,-6.441176519729197025e-04,-1.401926390826702118e-03,-6.296651554293930531e-04,-1.279979012906551361e-03,-6.440877914428710938e-04,-1.105200033634901047e-03,-6.326425354927778244e-04,-9.625115199014544487e-04,-5.436314968392252922e-04,-9.006902691908180714e-04,-3.925641358364373446e-04,-8.932135533541440964e-04,-2.433453046251088381e-04,-8.688018424436450005e-04,-1.677869004197418690e-04,-7.866683299653232098e-04,-1.704707392491400242e-04,-6.564725190401077271e-04,-1.993057812796905637e-04,-5.436229985207319260e-04,-1.821098849177360535e-04,-4.932153970003128052e-04,-8.811930456431582570e-05,-4.982755635865032673e-04,4.047903712489642203e-05,-5.087682511657476425e-04,1.367125078104436398e-04,-4.608065937645733356e-04,1.659864501561969519e-04,-3.472715616226196289e-04,1.572848123032599688e-04,-2.113831287715584040e-04,1.804512576200067997e-04,-1.090367441065609455e-04,2.806144475471228361e-04,-9.502422471996396780e-05,4.297870036680251360e-04,-1.934093306772410870e-04,3.295705246273428202e-04,-1.011353815556503832e-04,9.054032852873206139e-05,-3.264446058892644942e-05,5.004780541639775038e-05,-2.827877506206277758e-05,1.018938110064482316e-05,-3.748103699763305485e-05,-7.563044164271559566e-06,-3.955123975174501538e-05,-1.398747644998366013e-05,-3.543049024301581085e-05,-1.030693420034367591e-05,-2.779769965854939073e-05,-5.537351171369664371e-06,-2.013198900385759771e-05,-3.731379820237634704e-06,-1.549381340737454593e-05,-3.401703907002229244e-06,-1.303655790252378210e-05,-3.378006340426509269e-06,-1.178074489871505648e-05,-2.873806806746870279e-06,-1.096920823329128325e-05,-1.997708295675693080e-06,-9.551962648401968181e-06,-7.752050805720500648e-07,-7.491975793527672067e-06,3.070530283366679214e-07,-5.174365924176527187e-06,1.402136945216625463e-06,-2.910591319960076362e-06,2.711339902816689573e-06,-9.658213002694537863e-07,4.243900548317469656e-06,6.322999865915335249e-07,5.791084277007030323e-06,2.036470277744228952e-06,7.097078196238726377e-06,3.371567345311632380e-06,7.995602572918869555e-06,4.608377366821514443e-06,8.501098818669561297e-06,5.562949809245765209e-06,8.723820428713224828e-06,6.061910880816867575e-06,8.764889571466483176e-06,6.082213985791895539e-06,8.638658982818014920e-06,5.773607881565112621e-06,8.304193215735722333e-06,5.348522336134919897e-06,7.753198588034138083e-06,4.956483735440997407e-06,7.065953468554653227e-06,4.619603259925497696e-06,6.384686457749921829e-06,4.281295787222916260e-06,5.828461780765792355e-06,3.897186161339050159e-06,5.418819000624353066e-06,3.482799229459487833e-06,5.076829893368994817e-06,3.095952706644311547e-06,4.697511940321419388e-06,2.770437731669517234e-06,4.229279056744417176e-06,2.481095179973635823e-06,3.696067551572923549e-06,2.159476707674912177e-06,3.161283757435739972e-06,1.751146214701293502e-06,2.668529532456886955e-06,1.263349190594453830e-06,2.209349077020306140e-06,7.628735261278052349e-07,1.750232286212849431e-06,3.238855583731492516e-07,1.278573790841619484e-06,-1.676131766714661353e-08,8.267254543170565739e-07,-2.746158145328081446e-07,4.547930814169376390e-07,-4.837585834138735663e-07,2.041634843408246525e-07,-6.597906008209974971e-07,6.949969133529521059e-08,-7.863818609621375799e-07,4.663265329440946516e-09,-8.368477892872761004e-07,-3.919615565450840222e-08,-8.071175443546962924e-07,-8.178678001513617346e-08,-7.293683665920980275e-07,-1.106800198158452986e-07,-6.513132007057720330e-07,-1.078565929901742493e-07,-6.034958346390340012e-07,-7.444621985541743925e-08,-5.817859118906199001e-07,-3.486836774868606881e-08,-5.580810125138668809e-07,-1.567322982509722351e-08,-5.089673891234269831e-07,-2.121177899994108884e-08,-4.357118257303227438e-07,-3.146131533071638842e-08,-3.597972693114570575e-07,-2.085375783167364716e-08,-3.012759464127157116e-07,1.810938954349694541e-08,-2.608179556773393415e-07,6.712920708196179476e-08,-2.225000059752346715e-07,1.015123913816751156e-07,-1.685384347638319014e-07,1.120519144137688272e-07,-9.806510092857934069e-08,1.112391672108969942e-07,-3.039341933686046104e-08,1.162228215889626881e-07,1.226189993985826732e-08,1.272632346172031248e-07,1.991670295353742404e-08,1.242260054823418614e-07,-2.947287525145725340e-09,7.869989815389999421e-08,-2.011494792952817079e-08};

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __roombinauraliser_INTERNAL_H_INCLUDED__ */
