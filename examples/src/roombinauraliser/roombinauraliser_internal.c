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
 * @file: roombinauraliser_internal.c
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

#include "roombinauraliser_internal.h"

void roombinauraliser_setCodecStatus(void* const hBin, CODEC_STATUS newStatus)
{
    roombinauraliser_data *pData = (roombinauraliser_data*)(hBin);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}

void roombinauraliser_interpHRTFs
(
    void* const hBin,
    float azimuth_deg,
    float elevation_deg,
    float_complex h_intrp[MAX_NUM_INPUTS][HYBRID_BANDS][NUM_EARS]
)
{
    roombinauraliser_data *pData = (roombinauraliser_data*)(hBin);
    int i, band;
    int aziIndex, elevIndex, N_azi, idx3d;
    float_complex weights_cmplx[3], hrtf_fb3[MAX_NUM_INPUTS][NUM_EARS][3];
    float aziRes, elevRes, weights[3];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
     
    /* find closest pre-computed VBAP direction */
    aziRes = (float)pData->hrtf_vbapTableRes[0];
    elevRes = (float)pData->hrtf_vbapTableRes[1];
    N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
    aziIndex = (int)(matlab_fmodf(azimuth_deg + 180.0f, 360.0f) / aziRes + 0.5f);
    elevIndex = (int)((elevation_deg + 90.0f) / elevRes + 0.5f);
    if (!pData->VBAP_3d_FLAG)
        elevIndex = 0;
    idx3d = elevIndex * N_azi + aziIndex;
    for (i = 0; i < 3; i++)
        weights[i] = pData->hrtf_vbap_gtableComp[idx3d*3 + i];
    
    /* Interpolate */
    for (i = 0; i < 3; i++)
        weights_cmplx[i] = cmplxf(weights[i], 0.0f);
    for (int source=0; source < pData->nSources; source++){
        for (band = 0; band < HYBRID_BANDS; band++) {
            printf(" azi \t\t|  elev \t\t|  idx \n");
            for (i = 0; i < 3; i++){
                printf("%.3f \t\t| %.3f \t\t| %i \n", azimuth_deg, elevation_deg, pData->hrtf_vbap_gtableIdx[idx3d*3+i]);
                hrtf_fb3[source][0][i] = pData->hrtf_fb[band*NUM_EARS*(source+1)*(pData->N_hrir_dirs) + 0*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i] ];
                hrtf_fb3[source][1][i] = pData->hrtf_fb[band*NUM_EARS*(source+1)*(pData->N_hrir_dirs) + 1*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i] ];
            }
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, 1, 3, &calpha,
                        (float_complex*)hrtf_fb3[source], 3,
                        (float_complex*)weights_cmplx, 1, &cbeta,
                        (float_complex*)h_intrp[source][band], 1);
        }
    }
}

void roombinauraliser_initHRTFsAndGainTables(void* const hBin)
{
    roombinauraliser_data *pData = (roombinauraliser_data*)(hBin);
    int i, new_len;
    float* hrtf_vbap_gtable, *hrirs_resampled;//, *hrir_dirs_rad;
#ifdef SAF_ENABLE_SOFA_READER_MODULE
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
#endif
    
    strcpy(pData->progressBarText,"Loading HRIRs");
    pData->progressBar0_1 = 0.2f;
    
    /* load sofa file or load default hrir data */
#ifdef SAF_ENABLE_SOFA_READER_MODULE
    if(!pData->useDefaultHRIRsFLAG && pData->sofa_filepath!=NULL){
        /* Load SOFA file */
        //error = saf_sofa_open(&sofa, pData->sofa_filepath, SAF_SOFA_READER_OPTION_DEFAULT);
        error = saf_sofa_open_universal(&sofa, pData->sofa_filepath, SAF_SOFA_READER_OPTION_NETCDF, SAF_SOFA_READER_USECASE_BRIR);


        /* Load defaults instead */
        if(error!=SAF_SOFA_OK || sofa.nReceivers!=NUM_EARS){
            pData->useDefaultHRIRsFLAG = 1;
            saf_print_warning("Unable to load the specified SOFA file, or it contained something other than 2 channels. Using default HRIR data instead.");
        }
        else{
            /* Copy SOFA data */
            pData->hrir_loaded_fs = (int)sofa.DataSamplingRate;
            pData->hrir_loaded_len = sofa.DataLengthIR;
            pData->N_hrir_dirs = sofa.nSources;
            pData->nEmitters = sofa.nEmitters;
            pData->new_nSources = pData->nSources = sofa.nEmitters;
            pData->hrirs = realloc1d(pData->hrirs, pData->N_hrir_dirs*NUM_EARS*pData->nSources*(pData->hrir_loaded_len)*sizeof(float));
            memcpy(pData->hrirs, sofa.DataIR, pData->N_hrir_dirs*NUM_EARS*pData->nSources*(pData->hrir_loaded_len)*sizeof(float));
            pData->hrir_dirs_deg = realloc1d(pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));

            cblas_scopy(pData->N_hrir_dirs, sofa.ListenerView, 3, pData->hrir_dirs_deg, 2); /* azi */
            cblas_scopy(pData->N_hrir_dirs, &sofa.ListenerView[1], 3, &pData->hrir_dirs_deg[1], 2); /* elev */
            
            
            
            /* Set Emitters to Points specified in BRIR */
            if (!strcmp(sofa.EmitterPositionUnits, "metre")) { /* strcmp returns 0 if both strings are equal, therefor !strcmp = true */
                /* cartesian coordinates */
                for (int i=0; i<sofa.nEmitters; i++){
                    pData->src_dirs_xyz[i][0] = sofa.EmitterPosition[3*i+0];
                    pData->src_dirs_xyz[i][1] = sofa.EmitterPosition[3*i+1];
                    pData->src_dirs_xyz[i][2] = sofa.EmitterPosition[3*i+2];
                    
                    float temp_sph[3];
                    cart2sph(&sofa.EmitterPosition[3*i], 1, 1, temp_sph);
                    pData->src_dirs_deg[i][0] = temp_sph[0];
                    pData->src_dirs_deg[i][1] = temp_sph[1];
                }
            }
            else {
                /* spherical coordinates */
                for (int i=0; i<sofa.nEmitters; i++){
                    if (sofa.EmitterPosition[3*i+0] > 180)
                        pData->src_dirs_deg[i][0] = sofa.EmitterPosition[3*i+0] - 360.0;
                    else if (sofa.EmitterPosition[3*i+0] < -180)
                        pData->src_dirs_deg[i][0] = sofa.EmitterPosition[3*i+0] + 360.0;
                    else
                        pData->src_dirs_deg[i][0] = sofa.EmitterPosition[3*i+0];
                    
                    if (sofa.EmitterPosition[3 * i + 1] > 90)
                        pData->src_dirs_deg[i][1] = sofa.EmitterPosition[3 * i + 1] - 180;
                    else if (sofa.EmitterPosition[3 * i + 1] < -90)
                        pData->src_dirs_deg[i][1] = sofa.EmitterPosition[3 * i + 1] + 180.0;
                    else
                        pData->src_dirs_deg[i][1] = sofa.EmitterPosition[3*i+1];
                    
                    float temp_sph[3];
                    sph2cart(&sofa.EmitterPosition[3*i], 1, 1, temp_sph);
                    pData->src_dirs_xyz[i][0] = temp_sph[0];
                    pData->src_dirs_xyz[i][1] = temp_sph[1];
                    pData->src_dirs_xyz[i][2] = temp_sph[2];
                }
            }
        }

        /* Clean-up */
        saf_sofa_close(&sofa);
    }
#else
    pData->useDefaultHRIRsFLAG = 1; /* Can only load the default HRIR data */
#endif
    if(pData->useDefaultHRIRsFLAG){
        /* Copy default HRIR data */
        pData->hrir_loaded_fs = __default_hrir_fs;
        pData->hrir_loaded_len = __default_hrir_len;
        pData->N_hrir_dirs = __default_N_hrir_dirs;
        pData->hrirs = realloc1d(pData->hrirs, pData->N_hrir_dirs*NUM_EARS*(pData->hrir_loaded_len)*sizeof(float));
        memcpy(pData->hrirs, (float*)__default_hrirs, pData->N_hrir_dirs*NUM_EARS*(pData->hrir_loaded_len)*sizeof(float));
        pData->hrir_dirs_deg = realloc1d(pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
        memcpy(pData->hrir_dirs_deg, (float*)__default_hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
    }

    /* Convert from the 0..360 convention, to -180..180 */
    convert_0_360To_m180_180(pData->hrir_dirs_deg, pData->N_hrir_dirs);

    /* estimate the ITDs for each HRIR */
    strcpy(pData->progressBarText,"Estimating ITDs");
    pData->progressBar0_1 = 0.4f;
    pData->itds_s = realloc1d(pData->itds_s, pData->N_hrir_dirs*sizeof(float));
    // estimateITDs(pData->hrirs, pData->N_hrir_dirs, pData->hrir_loaded_len, pData->hrir_loaded_fs, pData->itds_s);

    /* Resample the HRIRs if needed */
    if(pData->hrir_loaded_fs!=pData->fs){
        strcpy(pData->progressBarText,"Resampling the HRIRs");
        pData->progressBar0_1 = 0.5f;
        hrirs_resampled = NULL;
        resampleHRIRs(pData->hrirs, pData->N_hrir_dirs, pData->hrir_loaded_len, pData->hrir_loaded_fs, pData->fs, 1, &hrirs_resampled, &new_len);
        pData->hrirs = realloc1d(pData->hrirs, pData->N_hrir_dirs*NUM_EARS*new_len*sizeof(float));
        cblas_scopy(pData->N_hrir_dirs*NUM_EARS*new_len, hrirs_resampled, 1, pData->hrirs, 1);
        free(hrirs_resampled);
        pData->hrir_runtime_fs = pData->fs;
        pData->hrir_runtime_len = new_len;
    }
    else{
        pData->hrir_runtime_fs = pData->hrir_loaded_fs;
        pData->hrir_runtime_len = pData->hrir_loaded_len;
    }
    
    /* generate VBAP gain table */
    strcpy(pData->progressBarText,"Generating interpolation table");
    pData->progressBar0_1 = 0.6f;
    hrtf_vbap_gtable = NULL;
    pData->hrtf_vbapTableRes[0] = 2;
    pData->hrtf_vbapTableRes[1] = 5;
    
    pData->VBAP_3d_FLAG = 1;

    float elevation_max = -90.0f;
    float elevation_min = +90.0f;
    for (int i = 1; i < pData->N_hrir_dirs; i += 2) {
        elevation_max = SAF_MAX(elevation_max, pData->hrir_dirs_deg[i]);
        elevation_min = SAF_MIN(elevation_min, pData->hrir_dirs_deg[i]);
    }

    if (fabsf((elevation_min + 90) / (180) - (elevation_max + 90) / (180)) < 1e-6) /* dont criticize the normalize*/
        pData->VBAP_3d_FLAG = 0;

    if (pData->VBAP_3d_FLAG) {
        generateVBAPgainTable3D(pData->hrir_dirs_deg, pData->N_hrir_dirs, pData->hrtf_vbapTableRes[0], pData->hrtf_vbapTableRes[1], 1, 0, 0.0f,
            &hrtf_vbap_gtable, &(pData->N_hrtf_vbap_gtable), &(pData->nTriangles));
    }
    else {
        generateVBAPgainTable2D(pData->hrir_dirs_deg, pData->N_hrir_dirs, pData->hrtf_vbapTableRes[0],
            &hrtf_vbap_gtable, &(pData->N_hrtf_vbap_gtable), &(pData->nTriangles));
    }

    if(hrtf_vbap_gtable==NULL){
        /* if generating vbap gain tabled failed, re-calculate with default HRIR set */
        pData->useDefaultHRIRsFLAG = 1;
        roombinauraliser_initHRTFsAndGainTables(hBin);
    }
    
    /* compress VBAP table (i.e. remove the zero elements) */
    
    pData->hrtf_vbap_gtableComp = realloc1d(pData->hrtf_vbap_gtableComp, pData->N_hrtf_vbap_gtable * 3 * sizeof(float));
    pData->hrtf_vbap_gtableIdx = realloc1d(pData->hrtf_vbap_gtableIdx, pData->N_hrtf_vbap_gtable * 3 * sizeof(int));
    compressVBAPgainTable3D(hrtf_vbap_gtable, pData->N_hrtf_vbap_gtable, pData->N_hrir_dirs, pData->hrtf_vbap_gtableComp, pData->hrtf_vbap_gtableIdx);
    /* --> 3D GainTable compression also works in 2D */

    
    /* convert hrirs to filterbank coefficients */
    pData->progressBar0_1 = 0.6f;
    pData->hrtf_fb = realloc1d(pData->hrtf_fb, HYBRID_BANDS * NUM_EARS * pData->nSources * (pData->N_hrir_dirs)*sizeof(float_complex));
    for (int source = 0; source < pData->nSources; source++)
        HRIRs2HRTFs_afSTFT(pData->hrirs+source*pData->N_hrir_dirs*HYBRID_BANDS * NUM_EARS, pData->N_hrir_dirs, pData->hrir_runtime_len, HOP_SIZE, 0, 1, pData->hrtf_fb+source*(HYBRID_BANDS * NUM_EARS*(pData->N_hrir_dirs)));

    /* HRIR pre-processing */
    if(pData->enableHRIRsDiffuseEQ){
        /* get integration weights */
        strcpy(pData->progressBarText,"Applying HRIR diffuse-field EQ");
        pData->progressBar0_1 = 0.9f;
        if(pData->N_hrir_dirs<=1000){//3600
//            hrir_dirs_rad = malloc1d(pData->N_hrir_dirs*2*sizeof(float));
//            memcpy(hrir_dirs_rad, pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
//            cblas_sscal(pData->N_hrir_dirs*2, SAF_PI/180.0f, hrir_dirs_rad, 1);
            pData->weights = realloc1d(pData->weights, pData->N_hrir_dirs*sizeof(float));
//            calculateGridWeights(hrir_dirs_rad, pData->N_hrir_dirs, -1, pData->weights);
//            free(hrir_dirs_rad);
            getVoronoiWeights(pData->hrir_dirs_deg, pData->N_hrir_dirs, 0, pData->weights);
        }
        else{
            pData->weights = realloc1d(pData->weights, pData->N_hrir_dirs*sizeof(float));
            for(int idx=0; idx < pData->N_hrir_dirs; idx++)
                pData->weights[idx] = 4.f*SAF_PI / (float)pData->N_hrir_dirs;
        }
        diffuseFieldEqualiseHRTFs(pData->N_hrir_dirs, pData->itds_s, pData->freqVector, HYBRID_BANDS, pData->weights, 1, 0, pData->hrtf_fb);
    }

    /* calculate magnitude responses */
    pData->hrtf_fb_mag = realloc1d(pData->hrtf_fb_mag, HYBRID_BANDS*NUM_EARS*(pData->N_hrir_dirs)*sizeof(float)); 
    for(i=0; i<HYBRID_BANDS*NUM_EARS* (pData->N_hrir_dirs); i++)
        pData->hrtf_fb_mag[i] = cabsf(pData->hrtf_fb[i]);

    /* The HRTFs should be re-interpolated */
    pData->recalc_hrtf_interpFLAG = 1;
    
    /* clean-up */
    free(hrtf_vbap_gtable);
    
    /* multiConv filter copy*/
    pData->nfilters = pData->nSources*NUM_EARS;
    int nSamples = pData->filter_length = pData->hrir_runtime_len;
    pData->filters = malloc1d(pData->nfilters*nSamples*sizeof(float));
    //cblas_sscal(nSamples*10*NUM_EARS*MAX_NUM_CHANNELS, 1.0f/1000.0f, pData->hrirs, 1);
    int dir = 1;
    for (int ear=0;ear<NUM_EARS; ear++)
        for (int src=0;src<pData->nSources; src++)
            memcpy(&(pData->filters[ear*src*nSamples]), &pData->hrirs[dir*ear*src*nSamples], nSamples*sizeof(float));
    
    if ((pData->reInitFilters == 1) && (pData->filters !=NULL)) {
        pData->reInitFilters = 2;
        if (pData->hMultiConv != NULL)
            saf_multiConv_destroy(&(pData->hMultiConv));
        saf_multiConv_create(&(pData->hMultiConv),
                             roombinauraliser_FRAME_SIZE,
                             pData->filters,
                             pData->hrir_runtime_len,
                             pData->nfilters,
                             pData->enablePartitionedConv);

        /* Resize buffers */
        pData->inputFrameTD  = (float**)realloc2d((void**)pData->inputFrameTD, MAX_NUM_CHANNELS, roombinauraliser_FRAME_SIZE, sizeof(float));
        pData->outframeTD = (float**)realloc2d((void**)pData->outframeTD, MAX_NUM_CHANNELS, roombinauraliser_FRAME_SIZE, sizeof(float));
        memset(FLATTEN2D(pData->inputFrameTD), 0, MAX_NUM_CHANNELS*roombinauraliser_FRAME_SIZE*sizeof(float));

        /* reset FIFO buffers */
        pData->FIFO_idx = 0;
        memset(pData->inFIFO, 0, MAX_NUM_CHANNELS*roombinauraliser_FRAME_SIZE*sizeof(float));
        memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*roombinauraliser_FRAME_SIZE*sizeof(float));

        pData->reInitFilters = 0;
    }
    
}

void roombinauraliser_initTFT
(
    void* const hBin
)
{
    roombinauraliser_data *pData = (roombinauraliser_data*)(hBin);
 
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), pData->new_nSources, NUM_EARS, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(pData->new_nSources!=pData->nSources){
        afSTFT_channelChange(pData->hSTFT, pData->new_nSources, NUM_EARS);
        afSTFT_clearBuffers(pData->hSTFT);
    }
    pData->nSources = pData->new_nSources;
}

void roombinauraliser_loadPreset
(
    SOURCE_CONFIG_PRESETS preset,
    float dirs_deg[MAX_NUM_INPUTS][2],
    int* newNCH,
    int* nDims
)
{
    float sum_elev;
    int ch, i, nCH;
    
    switch(preset){
        default:
        case SOURCE_CONFIG_PRESET_DEFAULT:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = 0.0f;
            break;
        case SOURCE_CONFIG_PRESET_MONO:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __mono_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_STEREO:
            nCH = 2;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __stereo_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_5PX:
            nCH = 5;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __5pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_7PX:
            nCH = 7;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __7pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_8PX:
            nCH = 8;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __8pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_9PX:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_10PX:
            nCH = 10;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __10pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_11PX:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_11PX_7_4:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_7_4_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_13PX:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __13pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_22PX:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __22pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_22P2_9_10_3:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9_10_3p2_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_MCC:
            nCH = 45;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCC_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_MCC_SUBSET:
            nCH = 37;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCCsubset_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_APAJA:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja_dirs_deg[ch][i];
            break; 
        case SOURCE_CONFIG_PRESET_AALTO_LR:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_LR_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_DTU_AVIL:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __DTU_AVIL_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_ZYLIA_LAB:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Zylia_Lab_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_4:
            nCH = 4;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_2_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_12:
            nCH = 12;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_4_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_24:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_6_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_36:
            nCH = 36;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_8_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_48:
            nCH = 48;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_9_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_60:
            nCH = 60;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_10_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_9:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_9_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_16:
            nCH = 16;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_16_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_25:
            nCH = 25;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_25_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_49:
            nCH = 49;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_49_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_64:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_64_dirs_deg[ch][i];
            break;
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_INPUTS; ch++){
        for(i=0; i<2; i++){
            dirs_deg[ch][i] = __default_LScoords64_rad[ch][i]* (180.0f/SAF_PI);
        }
    }
    
    /* For dynamically changing the number of TFT channels */
    (*newNCH) = nCH;
    
    /* estimate number of dimensions. (Obviously fails if using 2D setups thare are on an angle.
       However, in these cases, triangulation should fail and revert to 2D anyway) */
    sum_elev = 0.0f;
    for(i=0; i<nCH; i++){
        sum_elev += dirs_deg[i][1];
    }
    if(sum_elev < 0.01f)
        (*nDims) = 2;
    else
        (*nDims) = 3;
}
 










 
