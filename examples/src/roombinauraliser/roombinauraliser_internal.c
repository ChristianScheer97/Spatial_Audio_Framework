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
    INTERP_MODES mode,
    float azimuth_deg,
    float elevation_deg,
    float_complex h_intrp[MAX_NUM_INPUTS][HYBRID_BANDS][NUM_EARS]
)
{
    roombinauraliser_data *pData = (roombinauraliser_data*)(hBin);
    int i, band, source;
    int aziIndex, elevIndex, N_azi, idx3d;
    float_complex weights_cmplx[3], hrtf_fb3[MAX_NUM_INPUTS][NUM_EARS][3], ipd[MAX_NUM_INPUTS];
    float aziRes, elevRes, weights[3], itds3[MAX_NUM_INPUTS][3], itdInterp;
    float magnitudes3[MAX_NUM_INPUTS][HYBRID_BANDS][3][NUM_EARS], magInterp[MAX_NUM_INPUTS][HYBRID_BANDS][NUM_EARS];
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
    
    switch (mode) {
        case INTERP_TRI:
            /* Interpolate */
            for (i = 0; i < 3; i++)
                weights_cmplx[i] = cmplxf(weights[i], 0.0f);
            for (source=0; source < pData->nSources; source++){
                for (band = 0; band < HYBRID_BANDS; band++) {
                    for (i = 0; i < 3; i++){
                        hrtf_fb3[source][0][i] = pData->hrtf_fb[source][band*NUM_EARS*(pData->N_hrir_dirs) + 0*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                        hrtf_fb3[source][1][i] = pData->hrtf_fb[source][band*NUM_EARS*(pData->N_hrir_dirs) + 1*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                    }
                    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, 1, 3, &calpha,
                                (float_complex*)hrtf_fb3[source], 3,
                                (float_complex*)weights_cmplx, 1, &cbeta,
                                (float_complex*)h_intrp[source][band], 1);
                }
            }
            break;
        case INTERP_TRI_PS:
            /* retrieve the 3 itds and hrtf magnitudes */
            for (source=0; source<pData->nSources; source++){
                for (i = 0; i < 3; i++) {
                    itds3[source][i] = pData->itds_s[source][pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                    for (band = 0; band < HYBRID_BANDS; band++) {
                        magnitudes3[source][band][i][0] = pData->hrtf_fb_mag[source][band*NUM_EARS*(pData->N_hrir_dirs) + 0*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                        magnitudes3[source][band][i][1] = pData->hrtf_fb_mag[source][band*NUM_EARS*(pData->N_hrir_dirs) + 1*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                    }
                }
              
                /* interpolate hrtf magnitudes and itd */
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, 1, 3, 1.0f,
                            (float*)weights, 3,
                            (float*)itds3[source], 1, 0.0f,
                            &itdInterp, 1);
                for (band = 0; band < HYBRID_BANDS; band++)
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, 2, 3, 1.0f,
                                (float*)weights, 3,
                                (float*)magnitudes3[source][band], 2, 0.0f,
                                (float*)magInterp[source][band], 2);
                
                /* introduce interaural phase difference */
                for (band = 0; band < HYBRID_BANDS; band++) {
                    if(pData->freqVector[band]<1.5e3f)
                        ipd[source] = cmplxf(0.0f, (matlab_fmodf(2.0f*SAF_PI*(pData->freqVector[band]) * itdInterp + SAF_PI, 2.0f*SAF_PI) - SAF_PI)/2.0f);
                    else
                        ipd[source] = cmplxf(0.0f, 0.0f);
                    h_intrp[source][band][0] = crmulf(cexpf(ipd[source]),        magInterp[source][band][0]);
                    h_intrp[source][band][1] = crmulf(conjf(cexpf(ipd[source])), magInterp[source][band][1]);
                }
            }
            break;
    }
}

void roombinauraliser_initHRTFsAndGainTables(void* const hBin)
{
    roombinauraliser_data *pData = (roombinauraliser_data*)(hBin);
    int new_len;
    float* hrtf_vbap_gtable, **hrirs_resampled, *resample;//, *hrir_dirs_rad
    if (pData->reInitHRTFsAndGainTables == REINIT_NONE)
        return;
#ifdef SAF_ENABLE_SOFA_READER_MODULE
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
#endif
    
    strcpy(pData->progressBarText,"Loading BRIRs");
    strcpy(pData->progressBarTooltip,"Reading impulse response data from the specified SOFA file");
    pData->progressBar0_1 = 0.1f;
    
    /* load sofa file or load default hrir data */
#ifdef SAF_ENABLE_SOFA_READER_MODULE
    if(pData->reInitHRTFsAndGainTables == REINIT_FULL)
        if (!pData->useDefaultHRIRsFLAG && pData->sofa_filepath!=NULL){
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
            pData->nSources = sofa.nEmitters;
            pData->hrirs = (float**)malloc2d(pData->nSources, pData->N_hrir_dirs*NUM_EARS*(pData->hrir_loaded_len), sizeof(float));
            
            for(int cycle=0; cycle<pData->N_hrir_dirs*NUM_EARS; cycle++)
                for(int source=0; source<pData->nSources; source++)
                    memcpy(pData->hrirs[source]+pData->hrir_loaded_len*cycle,
                           &sofa.DataIR[(cycle*pData->nSources+source)*pData->hrir_loaded_len],
                           pData->hrir_loaded_len*sizeof(float));
                
            pData->hrir_dirs_deg = realloc1d(pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
            cblas_scopy(pData->N_hrir_dirs, sofa.ListenerView, 3, pData->hrir_dirs_deg, 2); /* azi */
            cblas_scopy(pData->N_hrir_dirs, &sofa.ListenerView[1], 3, &pData->hrir_dirs_deg[1], 2); /* elev */
            
            /* Set Emitters to Points specified in BRIR */
            if (!strcmp(sofa.EmitterPositionUnits, "metre")) { /* strcmp returns 0 if both strings are equal, therefor !strcmp = true */
                
                /* cartesian coordinates */
                for (int i=0; i<pData->nSources; i++){
                    pData->src_dirs_xyz[i][0] = sofa.EmitterPosition[3*i+0];
                    pData->src_dirs_xyz[i][1] = sofa.EmitterPosition[3*i+1];
                    pData->src_dirs_xyz[i][2] = sofa.EmitterPosition[3*i+2];
                    
                    /* also convert to spherical */
                    float temp_sph[3];
                    cart2sph(&sofa.EmitterPosition[3*i], 1, 1, temp_sph);
                    pData->src_dirs_deg[i][0] = temp_sph[0];
                    pData->src_dirs_deg[i][1] = temp_sph[1];
                }
            }
            else {
                /* spherical coordinates */
                for (int i=0; i<pData->nSources; i++) {
                    
                    /* azimuth */
                    if (sofa.EmitterPosition[3*i+0] > 180)
                        pData->src_dirs_deg[i][0] = sofa.EmitterPosition[3*i+0] - 360.0;
                    else if (sofa.EmitterPosition[3*i+0] < -180)
                        pData->src_dirs_deg[i][0] = sofa.EmitterPosition[3*i+0] + 360.0;
                    else
                        pData->src_dirs_deg[i][0] = sofa.EmitterPosition[3*i+0];
                    
                    /* elevation */
                    if (sofa.EmitterPosition[3 * i + 1] > 90)
                        pData->src_dirs_deg[i][1] = sofa.EmitterPosition[3 * i + 1] - 180;
                    else if (sofa.EmitterPosition[3 * i + 1] < -90)
                        pData->src_dirs_deg[i][1] = sofa.EmitterPosition[3 * i + 1] + 180.0;
                    else
                        pData->src_dirs_deg[i][1] = sofa.EmitterPosition[3*i+1];
                    
                    /* also convert to cartesian */
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
    if(pData->reInitHRTFsAndGainTables == REINIT_FULL && pData->useDefaultHRIRsFLAG){
        /* Build default BRIR from Binauraliser HRIR data */
        pData->hrir_loaded_fs = __default_hrir_fs;
        pData->hrir_loaded_len = __default_hrir_len;
        pData->N_hrir_dirs = __default_N_hrir_dirs;

        pData->hrirs = (float**)realloc2d((void**)pData->hrirs, 2, pData->N_hrir_dirs *2* (pData->hrir_loaded_len), sizeof(float));
        cblas_scopy(pData->N_hrir_dirs * NUM_EARS * pData->hrir_loaded_len, (float*)__default_hrirs, 1, pData->hrirs[0], 1);
        cblas_scopy(pData->N_hrir_dirs * NUM_EARS * pData->hrir_loaded_len, (float*)__default_hrirs, 1, pData->hrirs[1], 1);
        pData->nSources = 2;
        pData->src_dirs_xyz[0][0] = pData->src_dirs_xyz[1][0] = 2;
        pData->src_dirs_xyz[0][1] = pData->src_dirs_xyz[1][1] = 2;
        pData->src_dirs_xyz[0][2] = pData->src_dirs_xyz[1][2] = 0;
        pData->src_dirs_deg[0][0] = 33;
        pData->src_dirs_deg[1][0] = -33;
        pData->src_dirs_deg[0][1] = 0;
        pData->src_dirs_deg[1][1] = 0;
        
        pData->hrir_dirs_deg = realloc1d(pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
        memcpy(pData->hrir_dirs_deg, (float*)__default_hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));


    }
    if(pData->reInitHRTFsAndGainTables == REINIT_FULL) {
        /* Convert from the 0..360 convention, to -180..180 */
        convert_0_360To_m180_180(pData->hrir_dirs_deg, pData->N_hrir_dirs);
        
        /* estimate the ITDs for each HRIR */
        strcpy(pData->progressBarText,"Estimating ITDs");
        strcpy(pData->progressBarTooltip,"Calculating time difference between both ears for all sources and directions");
        pData->progressBar0_1 = 0.3f;
        
        pData->itds_s = (float**)realloc2d((void**)pData->itds_s, pData->nSources, pData->N_hrir_dirs, sizeof(float));
        
        /* truncate hrirs for faster cross-correlation processing */
        int truncated_len = 1000;
        if (pData->hrir_loaded_len > truncated_len) {
            float* hrirs_truncated = (float*)malloc1d(NUM_EARS*pData->N_hrir_dirs*truncated_len*sizeof(float));
            for (int source=0; source<pData->nSources; source++) {
                for (int dir=0; dir<pData->N_hrir_dirs; dir++)
                    for (int ear=0; ear<NUM_EARS; ear++)
                        memcpy(hrirs_truncated+(NUM_EARS*dir+ear)*truncated_len, pData->hrirs[source]+(NUM_EARS*dir+ear), truncated_len*sizeof(float));
                estimateITDs(hrirs_truncated, pData->N_hrir_dirs, truncated_len, pData->hrir_loaded_fs, pData->itds_s[source]);
            }
            free(hrirs_truncated);
        }
        else
            for (int source=0; source>pData->nSources; source++)
                estimateITDs(pData->hrirs[source], pData->N_hrir_dirs, pData->hrir_loaded_len, pData->hrir_loaded_fs, pData->itds_s[source]);
    }
    if(pData->reInitHRTFsAndGainTables == REINIT_FULL || pData->reInitHRTFsAndGainTables == REINIT_RESAMPLE) {
        /* Resample the HRIRs if needed */
        if(pData->hrir_loaded_fs!=pData->fs){
            hrirs_resampled = NULL;
            resample = NULL;
            for(int source=0; source<pData->nSources; source++){
                char buffer[PROGRESSBARTEXT_CHAR_LENGTH];
                snprintf(buffer, sizeof(buffer), "Resampling BRIRs (Source %d/%d)",source+1, pData->nSources);
                strcpy(pData->progressBarText, buffer);
                strcpy(pData->progressBarTooltip, "Resampling the impulse responses to match the DAW's sampling rate. This may take some time...");
                pData->progressBar0_1 = 0.5f+0.2f*(source)/pData->nSources;
                
                resampleHRIRs(pData->hrirs[source], pData->N_hrir_dirs, pData->hrir_loaded_len, pData->hrir_loaded_fs, pData->fs, 1, &resample, &new_len);
                if (hrirs_resampled == NULL)
                    hrirs_resampled = (float**)realloc2d((void**)hrirs_resampled, pData->nSources, pData->N_hrir_dirs*NUM_EARS*new_len,sizeof(float));
                hrirs_resampled[source] = resample;
            }
            pData->hrirs = (float**)realloc2d((void**)pData->hrirs, pData->nSources, pData->N_hrir_dirs*NUM_EARS*new_len, sizeof(float));
            for (int source=0; source<pData->nSources; source++)
                memcpy(pData->hrirs[source], hrirs_resampled[source], pData->N_hrir_dirs*NUM_EARS*new_len*sizeof(float));
            
            free(resample);
            free(hrirs_resampled);
            
            pData->hrir_runtime_fs = pData->fs;
            pData->hrir_loaded_fs = pData->fs; /* needed to enable sample rate switching without reloading sofa */
            pData->hrir_runtime_len = new_len;
        }
        else{
            pData->hrir_runtime_fs = pData->hrir_loaded_fs;
            pData->hrir_runtime_len = pData->hrir_loaded_len;
        }
    }
    
    if(pData->reInitHRTFsAndGainTables == REINIT_FULL) {
        /* generate VBAP gain table */
        strcpy(pData->progressBarText,"Generating interpolation table");
        strcpy(pData->progressBarTooltip,"Calculating VBAP weights and filterbank coefficients");
        pData->progressBar0_1 = 0.7f;
        
        hrtf_vbap_gtable = NULL;
        pData->hrtf_vbapTableRes[0] = 2;
        pData->hrtf_vbapTableRes[1] = 5;
        pData->VBAP_3d_FLAG = 1;
        
        float elevation_max = -90.0f;
        float elevation_min = +90.0f;
        for (int i = 1; i < pData->N_hrir_dirs; i += 2) { /* only compare elevation data, skip azimuth */
            elevation_max = SAF_MAX(elevation_max, pData->hrir_dirs_deg[i]);
            elevation_min = SAF_MIN(elevation_min, pData->hrir_dirs_deg[i]);
        }
        
        /* Differentiate between 3D and 2D VBAP */
        if /* 2D */ (fabsf((elevation_min + 90) / (180) - (elevation_max + 90) / (180)) < 1e-6) { /* dont criticize the normalize */
            pData->VBAP_3d_FLAG = 0;
            generateVBAPgainTable2D(pData->hrir_dirs_deg, pData->N_hrir_dirs, pData->hrtf_vbapTableRes[0],
                                    &hrtf_vbap_gtable, &(pData->N_hrtf_vbap_gtable), &(pData->nTriangles));
        }
        else /* 3D */
            generateVBAPgainTable3D(pData->hrir_dirs_deg, pData->N_hrir_dirs, pData->hrtf_vbapTableRes[0],
                                    pData->hrtf_vbapTableRes[1], 1, 0, 0.0f,
                                    &hrtf_vbap_gtable, &(pData->N_hrtf_vbap_gtable), &(pData->nTriangles));
        
        
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
        
        /* clean-up */
        free(hrtf_vbap_gtable);
        hrtf_vbap_gtable = NULL;
    }
    
    /* convert hrirs to filterbank coefficients */
    pData->progressBar0_1 = 0.85f;
    pData->hrtf_fb = (float_complex**)realloc2d((void**)pData->hrtf_fb, pData->nSources, HYBRID_BANDS * NUM_EARS * (pData->N_hrir_dirs), sizeof(float_complex));
    for (int source = 0; source < pData->nSources; source++)
        HRIRs2HRTFs_afSTFT(pData->hrirs[source], pData->N_hrir_dirs, pData->hrir_runtime_len, HOP_SIZE, 0, 1, pData->hrtf_fb[source]);

    /* HRIR pre-processing */

    /* Apply diffuse field equalisation */
    if (pData->enableHRIRsDiffuseEQ) {

        /* dummy head (FABIAN) diffuse field equalisation */
        if (pData->diffEqMode == DIFF_EQ_FABIAN_CTF) {
            strcpy(pData->progressBarText, "Applying dummy head diffuse-field EQ");
            strcpy(pData->progressBarTooltip, "Applying dummy head diffuse-field EQ");
            pData->progressBar0_1 = 0.95f;
            pData->N_samples_fabian_cir = 256;

            pData->fabian_cir = realloc1d(pData->fabian_cir, pData->N_samples_fabian_cir * sizeof(float));
            memcpy(pData->fabian_cir, &fabian_ir, pData->N_samples_fabian_cir * sizeof(float));
            pData->ctf_fb = (float_complex*)realloc1d(pData->ctf_fb , HYBRID_BANDS * sizeof(float_complex));

            // print all CIR samples
            for (int i = 0; i < pData->N_samples_fabian_cir; i++)
                printf("CIR[%d] = %.3f\n", i, pData->fabian_cir[i]);

            /* convert FABIAN dummy head cir to filter bank coefficients */
            //afAnalyse(pData->fabian_cir, pData->N_samples_fabian_cir, 1, HOP_SIZE, 0, 1, pData->ctf_fb);
            afSTFT_FIRtoFilterbankCoeffs((float*)pData->fabian_cir, 1, 1, pData->N_samples_fabian_cir, HOP_SIZE, 0, 1, pData->ctf_fb);

            /* perform equalisation */
            int source, band, ear, nd;
            for (source = 0; source < pData->nSources; source++)
                for (band = 0; band < HYBRID_BANDS; band++)
                    for (ear = 0; ear < NUM_EARS; ear++ )
                        for (nd = 0; nd < pData->N_hrir_dirs; nd++) {
                            pData->hrtf_fb[source][band*NUM_EARS*pData->N_hrir_dirs + ear*pData->N_hrir_dirs + nd] = ccmulf(
                                    pData->ctf_fb[band],
                                    pData->hrtf_fb[source][band*NUM_EARS*pData->N_hrir_dirs + ear*pData->N_hrir_dirs + nd]); }

            // print all HRTF coefficients of source 0 and direction 0 after equalisation for the right ear
            for (int i = 0; i < HYBRID_BANDS; i++)
                printf("HRTF[%d] = %.3f + %.3f i\n", i, crealf(pData->hrtf_fb[0][i*NUM_EARS*pData->N_hrir_dirs + 1]), cimagf(pData->hrtf_fb[0][i*NUM_EARS*pData->N_hrir_dirs + 1]));
        }

        /* equalise diffuse field with loaded BRIR data */
        else if (pData->diffEqMode == DIFF_EQ_BRIR_CTF) {
            strcpy(pData->progressBarText,"Applying BRIR diffuse-field EQ");
            strcpy(pData->progressBarTooltip,"Applying BRIR diffuse-field EQ");
            pData->progressBar0_1 = 0.95f;

            if(pData->N_hrir_dirs<=3600){
                pData->weights = realloc1d(pData->weights, pData->N_hrir_dirs*sizeof(float));
                float* hrir_dirs_rad = (float*) malloc1d(pData->N_hrir_dirs*2*sizeof(float));
                memcpy(hrir_dirs_rad, pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
                cblas_sscal(pData->N_hrir_dirs*2, SAF_PI/180.f, hrir_dirs_rad, 1);
                sphElev2incl(hrir_dirs_rad, pData->N_hrir_dirs, 0, hrir_dirs_rad);
                //for (int i=0; i<pData->N_hrir_dirs; i++)
                //    printf("%.3f, %.3f \n", hrir_dirs_rad[2*i], hrir_dirs_rad[2*i+1]);
                int supOrder = calculateGridWeights(hrir_dirs_rad, pData->N_hrir_dirs, -1, pData->weights);
                if(supOrder < 1){
                    if(pData->VBAP_3d_FLAG) {
                        saf_print_warning("Could not calculate grid weights");
                        free(pData->weights);
                        pData->weights = NULL;
                    }
                    else {
                        /* DEQ weight calculation for 2D BRIRs go brrr */
                        saf_print_warning("Could not calculate grid weights");
                        free(pData->weights);
                        pData->weights = NULL;
                    }
                }
            }
            else{
                saf_print_warning("Too many grid points to calculate grid weights. i.e., we're not assuming that the HRTF measurement grid was uniform.");
                free(pData->weights);
                pData->weights = NULL;
            }
            for (int source = 0; source < pData->nSources; source++)
                diffuseFieldEqualiseHRTFs(pData->N_hrir_dirs, pData->itds_s[source], pData->freqVector, HYBRID_BANDS, pData->weights, 1, 0, (float_complex*)pData->hrtf_fb[source]);
        }
    }
    /* calculate magnitude responses */
    pData->hrtf_fb_mag = (float**)realloc2d((void**)pData->hrtf_fb_mag, pData->nSources, HYBRID_BANDS * NUM_EARS * (pData->N_hrir_dirs), sizeof(float));
    for (int source = 0; source < pData->nSources; source++)
        for (int i = 0; i < HYBRID_BANDS * NUM_EARS * (pData->N_hrir_dirs); i++)
            pData->hrtf_fb_mag[source][i] = cabsf(pData->hrtf_fb[source][i]);
    
    /* The HRTFs should be re-interpolated */
    pData->recalc_hrtf_interpFLAG = 1;
}

void roombinauraliser_initTFT
(
    void* const hBin
)
{
    roombinauraliser_data *pData = (roombinauraliser_data*)(hBin);
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), pData->nSources, NUM_EARS, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
}
