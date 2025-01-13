#ifndef LSDEPTHFIRST_H
#define LSDEPTHFIRST_H

#include <math.h>
#include <stdio.h>				
#include <stdlib.h>
#include <string.h>
#include <errno.h>				
#include <Float.h>
#include <Windows.h>
#include "functions.h"

#define MAX_PYRAMID_LAYER (10)
#define MIN_TIES_NUM (6)
#define MIN_TIES_NUM_POLYNOMIAL (12)
#define SATURATED_VALUE_MSS (153) // 4/20/2020: changed from 153 to 225
#define SATURATED_VALUE_MSS_DENSEMATCHING (200) // 7/22/2022, change from 235
#define MIN_SAM_THRESHOLD (0.9799)

// v1.1.1 (11/3/2017): all functions with input image size parameters have the order of image width followed by image height

typedef enum { enum_TRANSLATION = 1, enum_AFFINE = 2, enum_POLYNOMIAL = 3, enum_RBFs = 4, enum_AUTO = 10 } enum_tranformation_type;
typedef enum { enum_NN = 1, enum_BL = 2 } enum_resampling_method;

typedef struct {
	// image pair to be registered
	// notes
	//	1) images 1 and 2 have the same size and resolution, and "short int" data type; for registration between Sentinel and Lansat, they should be roughly geo-registered and resampled to the same resolution, e.g. resample Landsat to 10m
	//  2) image 1 should be the one with higher raw resolution, e.g. should be Sentinel for registration between Sentinel and Landsat
	//  3) in v1.1.1, images do not have to be stored in the work directory, see InputRegistrationProjectFile()
	char pacWorkDir[2000];				// directory to output results
	char pacImageFileName1[2000];		// image 1 file name
	char pacImageFileName2[2000];		// image 2 file name
	char pacImageFilePathName1[2000];	// image 1 file full path name
	char pacImageFilePathName2[2000];	// image 2 file full path name
	char pacRegID[2000];				// ID composed of pacImageFileName1 and pacImageFileName2
	int iRowMax;						// lines
	int iColMax;						// samples

	// registration parameters
	enum_tranformation_type enum_TransformationType;		// translation (1), affine (2), or polynomial (3); suggested value = 2
	bool bIsImage2BaseImage;								// always = false, indicating image 1 is base image, register image 2 to image 1
	char pacTargetImageFilePathName[2000];					// optional parameter, which is the (non-base) image to be registered to base image; images 2 is usually the NIR band, and here it be set to other bands rather than NIR; should have the same size (iColMax x iRowMax) as images 1 and 2; if not set, no registered images will be genereated (only the tranformation coefficients are output)
	int iTargetImageResolution;								// (not used) corresponds to iRawResolution; not used as resolution differences are handled in ResampleSentinel_JP2()
	int iMaxOffsetOnTopLayer;								// maximum misregistration offset on top layer pyramid image, used in InitMatch(); suggested value = 1

	// pyramid
	int piLayerScales[MAX_PYRAMID_LAYER];	// scales of upper pyramid layers; suggested value = [2, 2, 2] or [3, 2, 2] for iLayerNum = 3; for example, if piLayerScales = [3, 2, 2] and iRawResolution = 10 (m), the resolutions of the four layers are 10m, 30m, 60m and 120m; can use two layers if the images are less than 1000 x 1000 pixels, such as the two AVHRR images in the "data2" example 
	int iLayerNum;							// number of upper pyramid layers; should be larger than 3, which means the pyramid have at least 4 layers; should work together with piLayerScales to ensure proper (not too small) image size on top layer such that sufficient POIs can be detected; suggested value = 3
	int iRawResolution;						// (not used in Windows version) raw spaital resolution of images 1 and 2 in unit of meters

	// matching parameters
	int iDefaultHalfWindow;					// image matching half window on top layer; suggested value = 10
	int iMaxHalfWindow;						// maximum image matching half window on top layer; suggested value = 30
	float fSAMThreshold;					// SAM (spectral angle mapper) similarity threshold used for judging succssful SAM-based least squares matching; suggested value = 0.995 for Sentinel-2A to Sentinel-2A and Sentinel-2A to Landsat-8 (resampled to the same resolution); can use lower values like 0.990 if the two images have low SNR ratio, such as the two AVHRR images in the "data2" example
	int iDenseMatchingGridSize;				// optional parameter; sampling interval for dense grid point matching on bottom layer; if set to 0, no dense matching will be conducted; suggested value = 6
	
	//  added for Windows version; used in DenseMatching()
	short int shtBandValueThreshold;		// do not do matching if band value is smaller than this threshold; used to filter out water; a value of 500 is genrally suitable for NIR band; can be set to 0 for non-NIR images on general land covers
	short int shtMeanDiffThreshold;			// do not do matching if the mean difference of two matching windows is larger than the threshold; suggested value = 400									

	bool bShortData;
	unsigned char ucBandValueThreshold;		// do not do matching if band value is smaller than this threshold; used to filter out water; a value of 500 is genrally suitable for NIR band; can be set to 0 for non-NIR images on general land covers
	unsigned char ucMeanDiffThreshold;			// do not do matching if the mean difference of two matching windows is larger than the threshold; suggested value = 400									

	double pdCoefs[MAX_COEFs_NUM];			// transformation coefficients from base image to target image (used to register target image to base image) fitted using depth-first tie pioints
	int n_tie;								// nubmer of depth-first tie points obtained in DepthFirstMatch()
	double dRMSE;							// transformation fitting RMSE using depth-first tie points
	
	int iRBFs_K;							// set in InputRegistrationProjectFile()
	int iRBFs_K_adp;						// set in GetRegistrationTransformation_DenseMatching() by calling FitRBFsTransform_TPS_Poly_v2()

	enum_resampling_method enum_ResamplingMethod;

	double pdCoefs_DenseMatching[MAX_COEFs_NUM];		// transformation coefficients from base image to target image fitted using dense matching tie pioints
	int n_tie_DenseMatching;				// nubmer of depth-first tie points obtained in DepthFirstMatch()
	double dRMSE_DenseMatching;
} LSDepthFirstRegistration;

////////////////////////////////////////////////////////////////
// 7/26/2021, uint8
bool DepthFirstImageRegistration_UC_simple(LSDepthFirstRegistration* pLSReg, char* regfilename, double* pdMeanShift_x, double* pdMeanShift_y, double* pdRMSE, int* pi_ntie, int iHalfSearchWindow = 1, int iProcessingMode = 1);

void BuildImagePyramidPerLayer_UC(char* pacWorkDir, char* pacImageFileName, int iColMax, int iRowMax, int iLayerScale, int iProcType, int* piColMax_Layer, int* piRowMax_Layer);
bool GetMaskAndBaseImageForPOI_UC(unsigned char* pucImg1, unsigned char* pucImg2, int iColMax, int iRowMax, unsigned char* pucMask, bool* pbImage1);
int DenseMatching_UC_v2(LSDepthFirstRegistration* pLSReg, int iLayer, int iSamplingInterval, int iHalfWindow, int iHalfSearchWindow);
void GenerateRegisteredImage_v1_1_UC(LSDepthFirstRegistration* pLSReg);

////////////////////////////////////////////////////////////////
// Call this function to do registration

void SetRBFs_K(LSDepthFirstRegistration* pLSReg, int K);

int GetPyramidLayerResolution(LSDepthFirstRegistration *pLSReg, int iLayer);
void GetPyramidLayerInfo(LSDepthFirstRegistration *pLSReg, bool bImage1, int iLayer, char *pacPathName, int *piColMax_Layer, int *piRowMax_Layer);
int GetCumulativeScaleFromBottom(LSDepthFirstRegistration *pLSReg, int iLayer);

bool CreateRegitrationProjectFile(char* pacWorkDir, char* pacImagePathName1, char* pacImagePathName2, int iColMax, int iRowMax, int iDenseMatchingSamplingInterval, int iTransformationType, int iResamplingMethod, float fSAMThreshold, char* pacRegFileName);

bool InputRegistrationProjectFile(LSDepthFirstRegistration *pLSReg, char *regfilename);

void BuildImagePyramidPerLayer(char *pacWorkDir, char *pacImageFileName, int iColMax, int iRowMax, int iLayerScale, int iProcType, int *piColMax_Layer, int *piRowMax_Layer);
void BuildImagePyramid(LSDepthFirstRegistration *pLSReg, bool bImage1);

void InputPOI_v1(LSDepthFirstRegistration *pLSReg, int iLayer, int *piPOICols, int *piPOIRows, int iPOINum);

int DepthFirstMatch(LSDepthFirstRegistration *pLSReg, int iLayer);

enum_tranformation_type AutomaticTransformationType(int ntie);
int GetRegistrationTransformation_DepthFirst(LSDepthFirstRegistration *pLSReg, bool bImage2ToImage1, enum_tranformation_type enum_TransformationType, double *pdCoefs, double *pdRMSE);
int GetRegistrationTransformation_DenseMatching(LSDepthFirstRegistration *pLSReg, bool bIsImage2BaseImage, enum_tranformation_type enum_TransformationType, int iLayer, int iSamplingInterval, char *pacDenseMatchingFileName, double *pdCoefs, double *pdRMSE);
void OutputTransformations_v2(LSDepthFirstRegistration *pLSReg);
void OutputRegistrationNumaricals(LSDepthFirstRegistration *pLSReg, double *ptr_a0_translation, double *ptr_b0_translation, double *ptr_RMSE, int *ptr_ntie);
void OutputFittingResidualMaps(LSDepthFirstRegistration* pLSReg, bool FitTransform = false); 

void GenerateRegisteredImage_v1_1(LSDepthFirstRegistration *pLSReg);

void GetRegisteredImageFileName(LSDepthFirstRegistration* pLSReg, char* pacImageFilePathName);

#endif
