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
 * @example roombinauraliser.h
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain.
 *
 * ### Files
 * roombinauraliser.h (include), roombinauraliser_internal.h, roombinauraliser.c,
 * roombinauraliser_internal.c
 * ### Include Header
 */

/**
 * @file: roombinauraliser.h
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

#ifndef __roombinauraliser_H_INCLUDED__
#define __roombinauraliser_H_INCLUDED__

#include "_common.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the roombinauraliser
 *
 * @param[in] phBin (&) address of roombinauraliser handle
 */
void roombinauraliser_create(void** const phBin);

/**
 * Destroys an instance of the roombinauraliser
 *
 * @param[in] phBin (&) address of roombinauraliser handle
 */
void roombinauraliser_destroy(void** const phBin);

/**
 * Initialises an instance of roombinauraliser with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hBin       roombinauraliser handle
 * @param[in] samplerate Host samplerate.
 */
void roombinauraliser_init(void* const hBin,
                       int samplerate);

/**
 * Intialises the codec variables, based on current global/user parameters
 *
 * @note This function is fully threadsafe. It can even be called periodically
 *       via a timer on one thread, while calling _process() on another thread.
 *       Since, if a set function is called (that warrants a re-init), then a
 *       flag is triggered internally and the next time this function is called,
 *       it will wait until the current process() function has completed before
 *       reinitialising the relevant parameters. If the _initCodec() takes
 *       longer than the time it takes for process() to be called again, then
 *       process() is simply bypassed until the codec is ready.
 * @note This function does nothing if no re-initialisations are required.
 *
 * @param[in] hBin roombinauraliser handle
 */
void roombinauraliser_initCodec(void* const hBin);

/**
 * Binauralises the input signals at the user specified directions
 *
 * @param[in] hBin      roombinauraliser handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 */
void roombinauraliser_process(void* const hBin,
                          const float *const * inputs,
                          float** const outputs,
                          int nInputs,
                          int nOutputs,
                          int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as roombinauraliser is currently configured, at next available opportunity.
 */
void roombinauraliser_refreshSettings(void* const hBin);
    
/** Sets the panning azimuth for a specific channel index, in DEGREES */
void roombinauraliser_setSourceAzi_deg(void* const hBin,
                                   int index,
                                   float newAzi_deg);

/** Sets the panning elevation for a specific channel index, in DEGREES */
void roombinauraliser_setSourceElev_deg(void* const hBin,
                                    int index,
                                    float newElev_deg);

/** Sets the number of input channels/sources to binauralise */
void roombinauraliser_setNumSources(void* const hBin, int new_nSources);

/**
 * Sets flag to dictate whether the default HRIRs in the Spatial_Audio_Framework
 * should be used (1), or a custom HRIR set loaded via a SOFA file (0).
 *
 * @note If the custom set fails to load correctly, roombinauraliser will revert to
 *       the defualt set. Use roombinauraliser_getUseDefaultHRIRsflag() to check if
 *       loading was successful.
 */
void roombinauraliser_setUseDefaultHRIRsflag(void* const hBin, int newState);

/**
 * Sets the file path for a .sofa file, in order to employ a custom HRIR set for
 * the decoding.
 *
 * @note If the custom set fails to load correctly, roombinauraliser will revert to
 *       the defualt set. Use roombinauraliser_getUseDefaultHRIRsflag() to check if
 *       loading was successful.
 *
 * @param[in] hBin roombinauraliser handle
 * @param[in] path File path to .sofa file (WITH file extension)
 */
void roombinauraliser_setSofaFilePath(void* const hBin, const char* path);

/** Enable (1) or disable (0) the diffuse-field EQ applied to the HRTFs */
void roombinauraliser_setEnableHRIRsDiffuseEQ(void* const hBin, int newState);

/** Loads an input preset (see #SOURCE_CONFIG_PRESETS enum) */
void roombinauraliser_setInputConfigPreset(void* const hBin, int newPresetID);

/** Sets the flag to enable/disable (1 or 0) rotation */
void roombinauraliser_setEnableRotation(void* const hBin, int newState);

/** Enables/disables partitioned convolution */
void roombinauraliser_setEnablePartConv(void* const hBin, int newState);

/** Sets the 'yaw' rotation angle, in DEGREES */
void roombinauraliser_setYaw(void* const hBin, float newYaw);

/** Sets the 'pitch' rotation angle, in DEGREES */
void roombinauraliser_setPitch(void* const hBin, float newPitch);

/** Sets the 'roll' rotation angle, in DEGREES */
void roombinauraliser_setRoll(void* const hBin, float newRoll);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'yaw' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void roombinauraliser_setFlipYaw(void* const hBin, int newState);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'pitch' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void roombinauraliser_setFlipPitch(void* const hBin, int newState);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'roll' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void roombinauraliser_setFlipRoll(void* const hBin, int newState);

/**
 * Sets a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw" (1)
 * rotation order.
 */
void roombinauraliser_setRPYflag(void* const hBin, int newState);

/** NOT IMPLEMENTED YET */
void roombinauraliser_setInterpMode(void* const hBin, int newMode);

/**
 * Sets gain factor for an input source.
 */
void roombinauraliser_setSourceGain(void* const hBin, int srcIdx, float newGain);

/**
 * Set a source to solo.
 */

void roombinauraliser_muteSource(void* const hBin, int srcIdx, bool muted);

/**
 * Mutes or unmutes a source.
 */
void roombinauraliser_setSourceSolo(void* const hBin, int srcIdx);

/**
 * Unsolo / unmute all sources.
 */
void roombinauraliser_setUnSolo(void* const hBin);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int roombinauraliser_getFrameSize(void);

/** Returns current codec status codec status (see #CODEC_STATUS enum) */
CODEC_STATUS roombinauraliser_getCodecStatus(void* const hBin);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * - 0: intialisation/processing has started
 * - 1: intialisation/processing has ended
 */
float roombinauraliser_getProgressBar0_1(void* const hBin);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void roombinauraliser_getProgressBarText(void* const hBin, char* text, char* tooltip);

/** Returns the source azimuth for a given index, in DEGREES */
float roombinauraliser_getSourceAzi_deg(void* const hBin, int index);

/** Returns the source elevation for a given index, in DEGREES */
float roombinauraliser_getSourceElev_deg(void* const hBin, int index);

/** Returns the number of inputs/sources in the current layout */
int roombinauraliser_getNumSources(void* const hBin);

/** Returns the maximum number of input sources supported by roombinauraliser */
int roombinauraliser_getMaxNumSources(void);

/** Returns the number of ears possessed by the average homo sapien */
int roombinauraliser_getNumEars(void);

/** Returns the number of directions in the currently used HRIR set */
int roombinauraliser_getNDirs(void* const hBin);

/**
 * Returns the number of triangular groupings (faces) returned by the Convex
 * Hull
 */
int roombinauraliser_getNTriangles(void* const hBin);

/** Returns the HRIR/HRTF azimuth for a given index, in DEGREES */
float roombinauraliser_getHRIRAzi_deg(void* const hBin, int index);

/** Returns the HRIR/HRTF elevation for a given index, in DEGREES */
float roombinauraliser_getHRIRElev_deg(void* const hBin, int index);

/** Returns the length of HRIRs in time-domain samples */
int roombinauraliser_getHRIRlength(void* const hBin);

/** Returns the HRIR sample rate */
int roombinauraliser_getHRIRsamplerate(void* const hBin);

/**
 * Returns the value of a flag used to dictate whether the default HRIRs in the
 * Spatial_Audio_Framework should be used (1), or a custom HRIR set loaded via a
 * SOFA file (0).
 *
 * @note If the custom set fails to load correctly, roombinauraliser will revert to
 *       the defualt set, so this will be '1'
 */
int roombinauraliser_getUseDefaultHRIRsflag(void* const hBin);

/**
 * Returns the file path for a .sofa file.
 *
 * @note If the custom set fails to load correctly, roombinauraliser will revert to
 *       the defualt set. Use 'roombinauraliser_getUseDefaultHRIRsflag()' to check
 *       if loading was successful.
 *
 * @param[in] hBin roombinauraliser handle
 * @returns        File path to .sofa file (WITH file extension)
 */
char* roombinauraliser_getSofaFilePath(void* const hBin);

/**
 * Returns the flag indicating whether the diffuse-field EQ applied to the HRTFs
 * is enabled (1) or disabled (0).
 */
int roombinauraliser_getEnableHRIRsDiffuseEQ(void* const hBin);

/** Returns the DAW/Host sample rate */
int roombinauraliser_getDAWsamplerate(void* const hBin);

/**
 * Returns the flag value which dictates whether to enable/disable sound-field
 * rotation (0: disabled, 1: enabled)
 */
int roombinauraliser_getEnableRotation(void* const hBin);

/** Returns the 'yaw' rotation angle, in DEGREES */
float roombinauraliser_getYaw(void* const hBin);

/** Returns the 'pitch' rotation angle, in DEGREES */
float roombinauraliser_getPitch(void* const hBin);

/** Returns the 'roll' rotation angle, in DEGREES */
float roombinauraliser_getRoll(void* const hBin);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int roombinauraliser_getFlipYaw(void* const hBin);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int roombinauraliser_getFlipPitch(void* const hBin);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int roombinauraliser_getFlipRoll(void* const hBin);

/**
 * Returns a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw"
 * (1) rotation order.
 */
int roombinauraliser_getRPYflag(void* const hBin);


/**
 * Returns the processing delay in samples (may be used for delay compensation
 * purposes)
 */
int roombinauraliser_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __roombinauraliser_H_INCLUDED__ */
