#include "LSDepthFirst.h"

/***
Created 7/28/2021
Simple version of registration with no depth first matching; directly do dense matching
iProcessingMode: 1 - just do the dense matching not based on depth-matching; 2 - generate registered image 
***/
bool DepthFirstImageRegistration_UC_simple(LSDepthFirstRegistration* pLSReg, char* regfilename, double* pdMeanShift_x, double* pdMeanShift_y, double* pdRMSE, int* pi_ntie, int iHalfSearchWindow, int iProcessingMode)
{
	int ntie_DenseMatching;

	pLSReg->bShortData = false;

	// 1. initialize project
	InputRegistrationProjectFile(pLSReg, regfilename);

	// 2. dense grid point matching on bottom layer; results stored at [pLSReg->pdCoefs_DenseMatching, pLSReg->n_tie_DenseMatching, pLSReg->dRMSE_DenseMatching]
	SetRBFs_K(pLSReg, 64);
	ntie_DenseMatching = DenseMatching_UC_v2(pLSReg, 0, pLSReg->iDenseMatchingGridSize, 31, iHalfSearchWindow);

	if (ntie_DenseMatching < MIN_TIES_NUM)
	{
		printf(" Insufficient dense matching tie points: %d\n\n", ntie_DenseMatching);
		return false;
	}
	*pdRMSE = pLSReg->dRMSE_DenseMatching;
	*pi_ntie = pLSReg->n_tie_DenseMatching;
	*pdMeanShift_x = 0;
	*pdMeanShift_y = 0;

	if (iProcessingMode == 1)
		return true;

	return true;
}

/***
Build image pyramids. New pyramid images are saved in "pyramid" sub-directory

Parameters
iProcType:	1 - input from base directory and output to "pyramids" sub-directory;
2 - input and output in "pyramids" direcotry sub-directory;
11/2/2017: fixed bug that iColMax and iRowMax were swapped

6/3/2019: return new layer sizes (previously returned as 0) when image exists
***/
void BuildImagePyramidPerLayer_UC(char* pacWorkDir, char* pacImageFileName, int iColMax, int iRowMax, int iLayerScale, int iProcType, int* piColMax_Layer, int* piRowMax_Layer)
{
	FILE* fout = NULL, * fin = NULL;
	char pacImageFileName_in[STRLEN] = "", pacImageFileName_out[STRLEN] = "";
	unsigned char* pucData = NULL, * pucDataNew = NULL, * pucData_Layer = NULL;
	int		FilterWidth;
	float* pfFilter = NULL;
	double dFilterSigma;
	int Row, Col, Row_Layer, Col_Layer;
	int iSamplingInterval;
	int iRowMax_Layer, iColMax_Layer;
	char pacDescription[STRLEN] = "";
	char  pac_null[2] = "";

	*piRowMax_Layer = 0;
	*piColMax_Layer = 0;

	// check if file exists
	if (iProcType == 1 || iProcType == 2)
	{
		sprintf_s(pacImageFileName_out, STRLEN, "%s/pyramids/%s_%d", pacWorkDir, pacImageFileName, iLayerScale);
		if (FileExist(pacImageFileName_out))
		{
			//printf(" File exists.");

			// added 6/3/2019 in MSS
			iRowMax_Layer = (int)(iRowMax * 1.0 / iLayerScale + 1e-10);
			iColMax_Layer = (int)(iColMax * 1.0 / iLayerScale + 1e-10);
			*piRowMax_Layer = iRowMax_Layer;
			*piColMax_Layer = iColMax_Layer;

			return;
		}
	}

	if (iProcType == 1)
		sprintf_s(pacImageFileName_in, STRLEN, "%s/%s", pacWorkDir, pacImageFileName); // in base directory
	else // iProcType == 2
		sprintf_s(pacImageFileName_in, STRLEN, "%s/pyramids/%s", pacWorkDir, pacImageFileName); // in "pyramid" sub directory

	pucData = (unsigned char*)ReadInputImage(pacImageFileName_in, iColMax, iRowMax, 0);

	iRowMax_Layer = (int)(iRowMax * 1.0 / iLayerScale + 1e-10);
	iColMax_Layer = (int)(iColMax * 1.0 / iLayerScale + 1e-10);

	if (!(pucDataNew = (unsigned char*)calloc(iRowMax * iColMax, sizeof(unsigned char)))
		|| !(pucData_Layer = (unsigned char*)calloc(iRowMax_Layer * iColMax_Layer, sizeof(unsigned char))))
	{
		printf("\nError in BuildImagePyramidPerLayer_UC(): insufficient memory.\n");
		scanf_s(" %d", &iColMax_Layer);
		exit(1);
	}

	dFilterSigma = sqrt(iLayerScale * 1.0 / 2);

	// get Gaussian filter
	pfFilter = GetGaussian(dFilterSigma, &FilterWidth);
	// apply filter
	Conv2same_UC(pucData, pucDataNew, iColMax, iRowMax, pfFilter, FilterWidth); // 11/2/2017: fixed bug that iColMax and iRowMax were swapped

	iSamplingInterval = iLayerScale;

	// resample
	Row_Layer = 0;
	for (Row = 0; Row <= iRowMax - iSamplingInterval; Row += iSamplingInterval)
	{
		Col_Layer = 0;
		for (Col = 0; Col <= iColMax - iSamplingInterval; Col += iSamplingInterval)
		{
			if (ISFILLVALUE(pucData[Row * iColMax + Col]))
				pucData_Layer[Row_Layer * iColMax_Layer + Col_Layer] = FILLVALUE;
			else
				pucData_Layer[Row_Layer * iColMax_Layer + Col_Layer] = pucDataNew[Row * iColMax + Col];

			Col_Layer += 1;
			if (Col_Layer >= iColMax_Layer)
				break;
		}
		Row_Layer += 1;
		if (Row_Layer >= iRowMax_Layer)
			break;
	}

	// output to "pyramids" sub directory
	sprintf_s(pacImageFileName_out, STRLEN, "%s/pyramids/%s_%d", pacWorkDir, pacImageFileName, iLayerScale);
	fout = WriteBinary(pacImageFileName_out);
	fwrite(pucData_Layer, sizeof(unsigned char), iRowMax_Layer * iColMax_Layer, fout);
	fclose(fout);

	// output ENVI hdr file
	sprintf_s(pacDescription, STRLEN, "\0");
	OutputEnviHDRFile(pacImageFileName_out, pacDescription, iColMax_Layer, iRowMax_Layer, 1, 1, pac_null);

	//	printf(" pyramid image output to\n  /%s\n", pacImageFileName);

	*piRowMax_Layer = iRowMax_Layer;
	*piColMax_Layer = iColMax_Layer;

	free(pucData);
	free(pucDataNew);
	free(pucData_Layer);
	free(pfFilter);

	return;
}

/***
Determine overlapping data mask and which image has smaller std
return false if no overlap
***/
bool GetMaskAndBaseImageForPOI_UC(unsigned char* pucImg1, unsigned char* pucImg2, int iColMax, int iRowMax, unsigned char* pucMask, bool* pbImage1)
{
	int i;
	unsigned char* pucImg1_ptr = NULL, * pucImg2_ptr = NULL;
	float fStd1, fStd2;
	int n;

	memset(pucMask, 0, iRowMax * iColMax * sizeof(unsigned char));
	*pbImage1 = true;

	pucImg1_ptr = pucImg1;
	pucImg2_ptr = pucImg2;

	n = 0;
	for (i = 0; i < iRowMax * iColMax; i++)
	{
		if (*pucImg1_ptr > 0 && *pucImg2_ptr > 0)
		{
			n += 1;
			pucMask[i] = 1;
		}

		pucImg1_ptr++;
		pucImg2_ptr++;
	}

	if (n == 0)
		return false;

	fStd1 = GetStd_UC(pucImg1, iRowMax * iColMax, pucMask);
	fStd2 = GetStd_UC(pucImg2, iRowMax * iColMax, pucMask);

	*pbImage1 = (fStd1 <= fStd2) ? true : false;

	return true;
}

//////////////////
void SetRBFs_K(LSDepthFirstRegistration* pLSReg, int K)
{
	pLSReg->iRBFs_K = K;
	return;
}

/***
Get resolution of a given pyramid layer
***/
int GetPyramidLayerResolution(LSDepthFirstRegistration *pLSReg, int iLayer)
{
	int i;
	int iResolution = pLSReg->iRawResolution;

	for (i = 0; i < iLayer; i++)
		iResolution *= pLSReg->piLayerScales[i];

	return iResolution;
}

/***
return pacPathName, piRowMax_Layer, piColMax_Layer
***/
void GetPyramidLayerInfo(LSDepthFirstRegistration *pLSReg, bool bImage1, int iLayer, char *pacPathName, int *piColMax_Layer, int *piRowMax_Layer)
{
	int iRowMax_Layer, iColMax_Layer;
	int i;

	// get project variables
	char pacWorkDir[STRLEN] = "", pacImageFileName[STRLEN] = "";
	int iRowMax, iColMax;
	int piLayerScales[MAX_PYRAMID_LAYER];
	int iLayerNum;
	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);
	if (bImage1)
		sprintf_s(pacImageFileName, STRLEN, "%s", pLSReg->pacImageFileName1);
	else
		sprintf_s(pacImageFileName, STRLEN, "%s", pLSReg->pacImageFileName2);
	iRowMax = pLSReg->iRowMax;
	iColMax = pLSReg->iColMax;
	iLayerNum = pLSReg->iLayerNum;
	memmove(piLayerScales, pLSReg->piLayerScales, iLayerNum*sizeof(int));

	*piRowMax_Layer = 0;
	*piColMax_Layer = 0;

	if (iLayer > iLayerNum)
	{
		printf("\nError in GetPyramidLayerInfo(): requested layer %d exceeds pyramid layer number %d.\n", iLayer + 1, iLayerNum + 1);
		scanf_s(" %d", &i);
		exit(1);
	}

	iRowMax_Layer = iRowMax;
	iColMax_Layer = iColMax;
	sprintf_s(pacPathName, STRLEN, "%s", pacImageFileName);
	for (i = 0; i < iLayer; i++)
	{
		iRowMax_Layer = (int)(iRowMax_Layer*1.0 / piLayerScales[i] + 1e-10);
		iColMax_Layer = (int)(iColMax_Layer*1.0 / piLayerScales[i] + 1e-10);
		sprintf_s(pacPathName, STRLEN, "%s_%d", pacPathName, piLayerScales[i]);
	}

	*piRowMax_Layer = iRowMax_Layer;
	*piColMax_Layer = iColMax_Layer;

	return;
}

/***
Get the cumulative scale from the bottom layer
***/
int GetCumulativeScaleFromBottom(LSDepthFirstRegistration *pLSReg, int iLayer)
{
	int iCurLayer;
	int iCumulativeScale = 1;
	for (iCurLayer = 0; iCurLayer < iLayer; iCurLayer++)
		iCumulativeScale *= pLSReg->piLayerScales[iCurLayer];

	return iCumulativeScale;
}

/***
Output registration ascii file of two tile images
iTransformationType: 1-translation; 2-affine; 3-polynomial
8/14/2019: adjust for MSS
***/
bool CreateRegitrationProjectFile(char* pacWorkDir, char* pacImagePathName1, char* pacImagePathName2, int iColMax, int iRowMax, int iDenseMatchingSamplingInterval, int iTransformationType, int iResamplingMethod, float fSAMThreshold, char* pacRegFileName)
{
	char pacRegFilePathName[STRLEN] = "";
	FILE* fout = NULL;
	//	char pacPathName[STRLEN];

	if (iTransformationType <= 0 || iTransformationType > 3)
	{
		iTransformationType = 2;
		printf("Invalid transforation type: %d. Reset to: 2 (affine).\n", iTransformationType);
	}
	sprintf_s(pacRegFilePathName, STRLEN, "%s/%s", pacWorkDir, pacRegFileName);

	// write to registration ascii file
	fout = Writetxt(pacRegFilePathName);

	// output dir
	fprintf(fout, "work directory                  %s\n", pacWorkDir);

	// image path names
	fprintf(fout, "image 1                         %s\n", pacImagePathName1);
	fprintf(fout, "image 2                         %s\n", pacImagePathName2);

	// others
	fprintf(fout, "samples                         %d\n", iColMax); //fixed wrong reversed output order for samples and lines
	fprintf(fout, "lines                           %d\n", iRowMax);
	fprintf(fout, "pyramid layers                  2\n");
	fprintf(fout, "pyramid layer scales            2 2\n");
	fprintf(fout, "raw resolution                  60\n");
	fprintf(fout, "maximum Offset On top layer     2\n");
	fprintf(fout, "tranformation type              %d\n", iTransformationType);
	fprintf(fout, "resampling method               %d\n", iResamplingMethod);
	fprintf(fout, "target image                    %s\n", pacImagePathName2);
	fprintf(fout, "target resolution               60\n");
	fprintf(fout, "default half window size        8\n");
	fprintf(fout, "maximum half window size        20\n");
	fprintf(fout, "SAM threshold                   %.3f\n", fSAMThreshold);
	fprintf(fout, "dense matching grid size        %d\n", iDenseMatchingSamplingInterval);
	fprintf(fout, "band value threshold            10\n");
	fprintf(fout, "mean diff threshold             150\n");

	fclose(fout);

	// create sub-directories
//	sprintf_s(pacPathName, STRLEN, "%s/matching", pacWorkDir); // POI and dig-matching results
//	_mkdir(pacPathName);

//	sprintf_s(pacPathName, STRLEN, "%s/pyramids", pacWorkDir);
//	_mkdir(pacPathName);

	return true;
}

/***
Inupt information from registration project file LSReg.txt located in the current working direcotry and create sub-directories
v1.1.1 (11/3/2017): the "pyramid" dir is made in the images data dirs rather than work dir, so images do not have to be stored in the work dir
1/7/2020: set pLSReg->iRBFs_K (= DEFAULT_RBFs_K * 2)
7/26/2021: set ucBandValueThreshold and ucMeanDiffThreshold if applicable
***/
bool InputRegistrationProjectFile(LSDepthFirstRegistration *pLSReg, char *regfilename)
{
	FILE* fin = NULL;
	char cLine[STRLEN] = "";
	char pacPathName[STRLEN] = "";
	char *pacPos = NULL;
	int i;
	char pacDataDir[STRLEN] = "", pacPyramidDir[STRLEN] = "";

	// check if file exists
	if (FileExist(regfilename) == false)
	{
		printf("%s does not exist.\n", regfilename);
		scanf_s(" %d", &i);
		exit(1);
		return false;
	}

	fopen_s(&fin, regfilename, "r");
	// read registration file
	while (fgets(cLine, STRLEN, fin) != NULL)
	{
		if (strstr(cLine, "work directory") != NULL)
			sscanf_s(cLine + strlen("work directory"), "%s", pLSReg->pacWorkDir, STRLEN);
		else if (strstr(cLine, "image 1") != NULL)
		{
			sscanf_s(cLine + strlen("image 1"), "%s", pLSReg->pacImageFilePathName1, STRLEN);
			pacPos = strrchr(pLSReg->pacImageFilePathName1, '/');
			pLSReg->pacImageFileName1[0] = '\0';
			strncpy_s(pLSReg->pacImageFileName1, STRLEN, pacPos + 1, strlen(pacPos) - 1);
			pLSReg->pacImageFileName1[pacPos + 1, strlen(pacPos) - 1] = '\0';
		}
		else if (strstr(cLine, "image 2") != NULL)
		{
			sscanf_s(cLine + strlen("image 2"), "%s", pLSReg->pacImageFilePathName2, STRLEN);
			pacPos = strrchr(pLSReg->pacImageFilePathName2, '/');
			pLSReg->pacImageFileName2[0] = '\0';
			strncpy_s(pLSReg->pacImageFileName2, STRLEN, pacPos + 1, strlen(pacPos) - 1);
			pLSReg->pacImageFileName2[pacPos + 1, strlen(pacPos) - 1] = '\0';
		}
		else if (strstr(cLine, "samples") != NULL)
			sscanf_s(cLine + strlen("samples"), "%d", &(pLSReg->iColMax));
		else if (strstr(cLine, "lines") != NULL)
			sscanf_s(cLine + strlen("lines"), "%d", &(pLSReg->iRowMax));
		else if (strstr(cLine, "pyramid layers") != NULL)
			sscanf_s(cLine + strlen("pyramid layers"), "%d", &(pLSReg->iLayerNum));
		else if (strstr(cLine, "pyramid layer scales") != NULL)
		{
			sscanf_s(cLine + strlen("pyramid layers scales"), "%d%d%d%d%d%d%d%d%d%d", pLSReg->piLayerScales, pLSReg->piLayerScales + 1, pLSReg->piLayerScales + 2, pLSReg->piLayerScales + 3, pLSReg->piLayerScales + 4, pLSReg->piLayerScales + 5, pLSReg->piLayerScales + 6, pLSReg->piLayerScales + 7, pLSReg->piLayerScales + 8, pLSReg->piLayerScales + 9);
			if (pLSReg->iLayerNum > 0)
				memset(pLSReg->piLayerScales + pLSReg->iLayerNum, 0, (MAX_PYRAMID_LAYER - pLSReg->iLayerNum)*sizeof(int));
		}
		else if (strstr(cLine, "raw resolution") != NULL)
			sscanf_s(cLine + strlen("raw resolution"), "%d", &(pLSReg->iRawResolution));
		else if (strstr(cLine, "maximum Offset On top layer") != NULL)
			sscanf_s(cLine + strlen("maximum Offset On top layer"), "%d", &(pLSReg->iMaxOffsetOnTopLayer));
		else if (strstr(cLine, "default half window size") != NULL)
			sscanf_s(cLine + strlen("default half window size"), "%d", &(pLSReg->iDefaultHalfWindow));
		else if (strstr(cLine, "maximum half window size") != NULL)
			sscanf_s(cLine + strlen("maximum half window size"), "%d", &(pLSReg->iMaxHalfWindow));
		else if (strstr(cLine, "SAM threshold") != NULL)
			sscanf_s(cLine + strlen("SAM threshold"), "%f", &(pLSReg->fSAMThreshold));
		else if (strstr(cLine, "tranformation type") != NULL)
		{
			sscanf_s(cLine + strlen("tranformation type"), "%d", &(pLSReg->enum_TransformationType));
			if (pLSReg->enum_TransformationType < 1 || pLSReg->enum_TransformationType > 4)
			{
				printf("Invalid transformation type %d.", pLSReg->enum_TransformationType);
				printf(" Reset to 2 (affine transformation).\n");
				pLSReg->enum_TransformationType = enum_AFFINE;
			}
			if (pLSReg->enum_TransformationType == 4)
				pLSReg->iRBFs_K = DEFAULT_RBFs_K;
			else
				pLSReg->iRBFs_K = 0;
		}
		else if (strstr(cLine, "resampling method") != NULL)
		{
			sscanf_s(cLine + strlen("resampling method"), "%d", &(pLSReg->enum_ResamplingMethod));
			if (pLSReg->enum_ResamplingMethod < 1 || pLSReg->enum_ResamplingMethod > 2)
			{
				printf("Invalid resampling method %d.", pLSReg->enum_ResamplingMethod);
				printf(" Reset to 2 (bilinear).\n");
				pLSReg->enum_ResamplingMethod = enum_BL;
			}
		}
		else if (strstr(cLine, "target image") != NULL)
		{
			if (sscanf_s(cLine + strlen("target image"), "%s", pLSReg->pacTargetImageFilePathName, STRLEN) < 0)
			{
				// no target image name found, set to nil
				sprintf_s(pLSReg->pacTargetImageFilePathName, STRLEN, "\0");
			}
		}
		else if (strstr(cLine, "target resolution") != NULL)
			sscanf_s(cLine + strlen("target resolution"), "%d", &(pLSReg->iTargetImageResolution));
		else if (strstr(cLine, "dense matching grid size") != NULL)
		{
			if (sscanf_s(cLine + strlen("dense matching grid size"), "%d", &(pLSReg->iDenseMatchingGridSize)) < 0)
			{
				// no dense matching grid size value found, set to 0
				pLSReg->iDenseMatchingGridSize = 0;
			}
		}
		else if (strstr(cLine, "band value threshold") != NULL)
		{
			sscanf_s(cLine + strlen("band value threshold"), "%hi", &(pLSReg->shtBandValueThreshold));
			if (pLSReg->bShortData == false)
				pLSReg->ucBandValueThreshold = (unsigned char)(pLSReg->shtBandValueThreshold);
		}
		else if (strstr(cLine, "mean diff threshold") != NULL)
		{
			sscanf_s(cLine + strlen("mean diff threshold"), "%hi", &(pLSReg->shtMeanDiffThreshold));
			if (pLSReg->bShortData == false)
				pLSReg->ucMeanDiffThreshold = (unsigned char)(pLSReg->shtMeanDiffThreshold);
		}
		else
			continue;
	}
	fclose(fin);

	pLSReg->bIsImage2BaseImage = false; // image 1 is always base image

	// set pacRegID
	sprintf_s(pLSReg->pacRegID, STRLEN, "%s_%s", pLSReg->pacImageFileName1, pLSReg->pacImageFileName2);

	// set output transformation fitting results
	memset(pLSReg->pdCoefs, 0, MAX_COEFs_NUM * sizeof(double));
	pLSReg->dRMSE = 0;
	pLSReg->n_tie = 0;

	memset(pLSReg->pdCoefs_DenseMatching, 0, MAX_COEFs_NUM * sizeof(double));
	pLSReg->dRMSE_DenseMatching = 0;
	pLSReg->n_tie_DenseMatching = 0;

	_mkdir(pLSReg->pacWorkDir); // added for MSS
								
	// create sub-directories (added for windows version)
	sprintf_s(pacPathName, STRLEN, "%s/matching", pLSReg->pacWorkDir); // POI and dig-matching results
	_mkdir(pacPathName);

	// create sub-directories (added for windows version)
//	sprintf_s(pacPathName, STRLEN, "%s/RBFs_stats", pLSReg->pacWorkDir); // POI and dig-matching results
//	_mkdir(pacPathName);

	// set directory for pyramids
	// 11/3/2017, the "pyramid" dir is made in the images data dirs rather than work dir
	GetFileDir(pLSReg->pacImageFilePathName1, pacDataDir);
	sprintf_s(pacPyramidDir, STRLEN, "%s/pyramids", pacDataDir);
	_mkdir(pacPyramidDir);

	GetFileDir(pLSReg->pacImageFilePathName2, pacDataDir);
	sprintf_s(pacPyramidDir, STRLEN, "%s/pyramids", pacDataDir);
	_mkdir(pacPyramidDir);

	sprintf_s(pacPathName, STRLEN, "%s/results", pLSReg->pacWorkDir);
	_mkdir(pacPathName);

	return true;
}

/***
Build image pyramid
***/
void BuildImagePyramid(LSDepthFirstRegistration *pLSReg, bool bImage1)
{
	int iLayer;
	int iRowMax_Layer_cur, iColMax_Layer_cur, iRowMax_Layer_new, iColMax_Layer_new;
	char pacImageFileName_cur[STRLEN] = "";
	int iLayerScale;

	// project variables
	char pacWorkDir[STRLEN] = "", pacImageFileName[STRLEN] = "", pacImageFilePathName[STRLEN] = "";
	int iRowMax, iColMax;
	int piLayerScales[MAX_PYRAMID_LAYER];
	int iLayerNum;

	char pacDataDir[STRLEN] = "", pacPyramidDir[STRLEN] = "";
	char *pacPos = NULL;

	// get project variables
	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);
	if (bImage1)
	{
		sprintf_s(pacImageFileName, STRLEN, "%s", pLSReg->pacImageFileName1);
		sprintf_s(pacImageFilePathName, STRLEN, "%s", pLSReg->pacImageFilePathName1);
	}
	else
	{
		sprintf_s(pacImageFileName, STRLEN, "%s", pLSReg->pacImageFileName2);
		sprintf_s(pacImageFilePathName, STRLEN, "%s", pLSReg->pacImageFilePathName2);
	}

	iRowMax = pLSReg->iRowMax;
	iColMax = pLSReg->iColMax;
	iLayerNum = pLSReg->iLayerNum;
	memmove(piLayerScales, pLSReg->piLayerScales, iLayerNum*sizeof(int));

	iRowMax_Layer_new = iRowMax;
	iColMax_Layer_new = iColMax;
	sprintf_s(pacImageFileName_cur, STRLEN, "%s", pacImageFileName);

	// set directory for pyramids
	GetFileDir(pacImageFilePathName, pacDataDir);
	sprintf_s(pacPyramidDir, STRLEN, "%s/pyramids", pacDataDir);

	// generate pyramids
	//	printf("Building pyramid for %s...", pacImageFileName);
	for (iLayer = 0; iLayer < iLayerNum; iLayer++)
	{
		iLayerScale = piLayerScales[iLayer];

		iRowMax_Layer_cur = iRowMax_Layer_new;
		iColMax_Layer_cur = iColMax_Layer_new;

		if (pLSReg->bShortData == true)
			BuildImagePyramidPerLayer(pacDataDir, pacImageFileName_cur, iColMax_Layer_cur, iRowMax_Layer_cur, iLayerScale, (iLayer == 0) ? 1 : 2, &iColMax_Layer_new, &iRowMax_Layer_new);
		else
			BuildImagePyramidPerLayer_UC(pacDataDir, pacImageFileName_cur, iColMax_Layer_cur, iRowMax_Layer_cur, iLayerScale, (iLayer == 0) ? 1 : 2, &iColMax_Layer_new, &iRowMax_Layer_new);

		sprintf_s(pacImageFileName_cur, STRLEN, "%s_%d", pacImageFileName_cur, iLayerScale);

		//printf(" %d m done.\n", GetPyramidLayerResolution(pLSReg, iLayer+1));
	}
	//	printf("done.\n");

	return;
}

/***
Build image pyramids. New pyramid images are saved in "pyramid" sub-directory

Parameters
iProcType:	1 - input from base directory and output to "pyramids" sub-directory;
2 - input and output in "pyramids" direcotry sub-directory;
11/2/2017: fixed bug that iColMax and iRowMax were swapped

6/3/2019: return new layer sizes (previously returned as 0) when image exists
***/
void BuildImagePyramidPerLayer(char *pacWorkDir, char *pacImageFileName, int iColMax, int iRowMax, int iLayerScale, int iProcType, int *piColMax_Layer, int *piRowMax_Layer)
{
	FILE *fout = NULL, *fin = NULL;
	char pacImageFileName_in[STRLEN] = "", pacImageFileName_out[STRLEN] = "";
	short int *pshtData = NULL, *pshtDataNew = NULL, *pshtData_Layer = NULL;
	int		FilterWidth;
	float	*pfFilter = NULL;
	double dFilterSigma;
	int Row, Col, Row_Layer, Col_Layer;
	int iSamplingInterval;
	int iRowMax_Layer, iColMax_Layer;
	char pacDescription[STRLEN] = "";
	char  pac_null[2] = "";

	*piRowMax_Layer = 0;
	*piColMax_Layer = 0;

	// check if file exists
	if (iProcType == 1 || iProcType == 2)
	{
		sprintf_s(pacImageFileName_out, STRLEN, "%s/pyramids/%s_%d", pacWorkDir, pacImageFileName, iLayerScale);
		if (FileExist(pacImageFileName_out))
		{
			//printf(" File exists.");

			// added 6/3/2019 in MSS
			iRowMax_Layer = (int)(iRowMax*1.0 / iLayerScale + 1e-10);
			iColMax_Layer = (int)(iColMax*1.0 / iLayerScale + 1e-10);
			*piRowMax_Layer = iRowMax_Layer;
			*piColMax_Layer = iColMax_Layer;

			return;
		}
	}

	if (iProcType == 1)
		sprintf_s(pacImageFileName_in, STRLEN, "%s/%s", pacWorkDir, pacImageFileName); // in base directory
	else // iProcType == 2
		sprintf_s(pacImageFileName_in, STRLEN, "%s/pyramids/%s", pacWorkDir, pacImageFileName); // in "pyramid" sub directory

	pshtData = (short int*)ReadInputImage(pacImageFileName_in, iColMax, iRowMax, 1);

	iRowMax_Layer = (int)(iRowMax*1.0 / iLayerScale + 1e-10);
	iColMax_Layer = (int)(iColMax*1.0 / iLayerScale + 1e-10);

	if (!(pshtDataNew = (short int*)calloc(iRowMax * iColMax, sizeof(short int)))
		|| !(pshtData_Layer = (short int*)calloc(iRowMax_Layer * iColMax_Layer, sizeof(short int))))
	{
		printf("\nError in BuildImagePyramidPerLayer(): insufficient memory.\n");
		scanf_s(" %d", &iColMax_Layer);
		exit(1);
	}

	dFilterSigma = sqrt(iLayerScale*1.0 / 2);

	// get Gaussian filter
	pfFilter = GetGaussian(dFilterSigma, &FilterWidth);
	// apply filter
	Conv2same(pshtData, pshtDataNew, iColMax, iRowMax, pfFilter, FilterWidth); // 11/2/2017: fixed bug that iColMax and iRowMax were swapped

	iSamplingInterval = iLayerScale;

	// resample
	Row_Layer = 0;
	for (Row = 0; Row <= iRowMax - iSamplingInterval; Row += iSamplingInterval)
	{
		Col_Layer = 0;
		for (Col = 0; Col <= iColMax - iSamplingInterval; Col += iSamplingInterval)
		{
			if (ISFILLVALUE(pshtData[Row*iColMax + Col]))
				pshtData_Layer[Row_Layer*iColMax_Layer + Col_Layer] = FILLVALUE;
			else
				pshtData_Layer[Row_Layer*iColMax_Layer + Col_Layer] = pshtDataNew[Row*iColMax + Col];

			Col_Layer += 1;
			if (Col_Layer >= iColMax_Layer)
				break;
		}
		Row_Layer += 1;
		if (Row_Layer >= iRowMax_Layer)
			break;
	}

	// output to "pyramids" sub directory
	sprintf_s(pacImageFileName_out, STRLEN, "%s/pyramids/%s_%d", pacWorkDir, pacImageFileName, iLayerScale);
	fout = WriteBinary(pacImageFileName_out);
	fwrite(pshtData_Layer, sizeof(short int), iRowMax_Layer*iColMax_Layer, fout);
	fclose(fout);

	// output ENVI hdr file
	sprintf_s(pacDescription, STRLEN, "\0");
	OutputEnviHDRFile(pacImageFileName_out, pacDescription, iColMax_Layer, iRowMax_Layer, 1, 2, pac_null);

	//	printf(" pyramid image output to\n  /%s\n", pacImageFileName);

	*piRowMax_Layer = iRowMax_Layer;
	*piColMax_Layer = iColMax_Layer;

	free(pshtData);
	free(pshtDataNew);
	free(pshtData_Layer);
	free(pfFilter);

	return;
}

/***
Input POIs from file in "matching" sub-directory
note: piPOIRows and piPOICols are pre-allocated
***/
void InputPOI_v1(LSDepthFirstRegistration *pLSReg, int iLayer, int *piPOICols, int *piPOIRows, int iPOINum)
{
	char pacWorkDir[STRLEN] = "", pacImageFileName[STRLEN] = "";
	int i;
	FILE *fin = NULL;

	// get POI file name
	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);
	sprintf_s(pacImageFileName, STRLEN, "%s/matching/poi_layer%d.txt", pacWorkDir, iLayer);

	if (FileExist(pacImageFileName) == false)
	{
		printf("Error in InputPOI_v1(): %s does not exist.", pacImageFileName);
		scanf_s(" %d", &i);
		exit(1);
	}

	memset(piPOIRows, 0, iPOINum*sizeof(int));
	memset(piPOICols, 0, iPOINum*sizeof(int));

	fin = Readtxt(pacImageFileName);
	for (i = 0; i < iPOINum; i++)
		fscanf_s(fin, "%d%d", piPOICols + i, piPOIRows + i);
	fclose(fin);

	return;
}

/***
Depth-first matching for initially matched points obtained by InitMatch(); output matching results on bottom layer to poi_layer*_matching_DF.txt in "matching" sub-directory

Parameters
iLayer:		the layer where depth-first matching starts; default is the number of pyramid layers indicating top layer; can be lower layers
***/
int DepthFirstMatch(LSDepthFirstRegistration *pLSReg, int iLayer)
{
	char pacWorkDir[STRLEN] = "", pacImageFileName1[STRLEN] = "", pacImageFileName2[STRLEN] = "";
	int iRowMax, iColMax;
	char pacImageFileName[STRLEN] = "";

	short int *pshtImg2 = NULL, *pshtImg1 = NULL;
	FILE *fin = NULL, *fout = NULL;
	double tmp_value;
	int ntie, n;
	int *piPreRows2 = NULL, *piPreCols2 = NULL;
	float *pfPreRows1 = NULL, *pfPreCols1 = NULL;
	float *pfPreCoefs = NULL;
	float *pfNewRows1 = NULL, *pfNewCols1 = NULL;		// matched coordinates on all lower layers
	float *pfNewCoefs = NULL;
	int iLayerScale, piLayerScales[MAX_PYRAMID_LAYER];
	int iCumulativeScale;
	unsigned char *pucTiesFlag = NULL;					// set to 1 if not passing depth-first match

	// variables for image matching
	float corr;
	float x1n, y1n, x2n, y2n;
	float xdif, ydif;
	int ws, w, ww;
	int h, max_h;
	int h_default, max_h_default;
	int x2, y2;
	float x1, y1;
	int ncol, nrow;
	double diff_threshold;
	float corr_threshold;

	int iCurLayer;
	int ntie_new;

	double pdCoefs[MAX_COEFs_NUM], dRMSE;

	char pacDataDir[STRLEN] = "";

	int value_threshold = 25; // changed to 25 for MSS.  500; // do not do matching if band value is smaller than this threshold; used to filter out water; a value of 500 is used for NIR band

	if (iLayer <= 0) // bug fix: changed from 1 to 0 in MSS.
		return 0;

	memmove(piLayerScales, pLSReg->piLayerScales, pLSReg->iLayerNum*sizeof(int));
	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);

	// check if output file exits
	sprintf_s(pacImageFileName, STRLEN, "%s/matching/poi_layer%d_matching_DF.txt", pacWorkDir, iLayer);
/*	if (FileExist(pacImageFileName)) // disabled for MSS
	{
		printf(" File exists.");
		ntie_new = Gettxtlines(pacImageFileName);
		printf(" %d tie points pass depth-first matching.\n", ntie_new);

		// save transformation fitting results
		ntie = GetRegistrationTransformation_DepthFirst(pLSReg, pLSReg->bIsImage2BaseImage, enum_AUTO, pdCoefs, &dRMSE);
		memmove(pLSReg->pdCoefs, pdCoefs, MAX_COEFs_NUM * sizeof(double));
		pLSReg->n_tie = ntie;
		pLSReg->dRMSE = dRMSE;

		//printf(" RMSE:\t%.3f pixels\n", dRMSE);
		return ntie_new;
	}*/

	// get initially matched points number
	sprintf_s(pacImageFileName, STRLEN, "%s/matching/poi_layer%d_matching_init.txt", pacWorkDir, iLayer);
	ntie = Gettxtlines(pacImageFileName);

	if (!(piPreRows2 = (int*)calloc(ntie, sizeof(int)))
		|| !(piPreCols2 = (int*)calloc(ntie, sizeof(int)))
		|| !(pfPreRows1 = (float*)calloc(ntie, sizeof(float)))
		|| !(pfPreCols1 = (float*)calloc(ntie, sizeof(float)))
		|| !(pfPreCoefs = (float*)calloc(ntie, sizeof(float)))
		|| !(pfNewRows1 = (float*)calloc(ntie*(iLayer + 1), sizeof(float)))
		|| !(pfNewCols1 = (float*)calloc(ntie*(iLayer + 1), sizeof(float)))
		|| !(pfNewCoefs = (float*)calloc(ntie*(iLayer + 1), sizeof(float)))
		|| !(pucTiesFlag = (unsigned char*)calloc(ntie, sizeof(unsigned char))))
	{
		printf("\nError in DepthFirstMatch(): insufficient memory.\n");
		scanf_s(" %d", &n);
		exit(1);
	}

	// input initially matched points on top level
	fin = Readtxt(pacImageFileName);
	for (n = 0; n < ntie; n++)
		fscanf_s(fin, "%d%d%f%f%lf%lf%f", piPreCols2 + n, piPreRows2 + n, pfPreCols1 + n, pfPreRows1 + n, &tmp_value, &tmp_value, pfPreCoefs + n);
	fclose(fin);

	// pdNewRows1, pdNewCols1 and pdNewCoefs store results on all layers including the previous layer "iLayer"
	// save results on top layer
	for (n = 0; n < ntie; n++)
	{
		pfNewRows1[n*(iLayer + 1)] = pfPreRows1[n];
		pfNewCols1[n*(iLayer + 1)] = pfPreCols1[n];
		pfNewCoefs[n*(iLayer + 1)] = pfPreCoefs[n];
	}

	// set matching parameters
	h_default = pLSReg->iDefaultHalfWindow;
	max_h_default = pLSReg->iMaxHalfWindow;
	corr_threshold = pLSReg->fSAMThreshold;

	// start depth-first matching
	ntie_new = ntie;
	memset(pucTiesFlag, 0, ntie*sizeof(unsigned char));
	iCumulativeScale = 1;
	for (iCurLayer = iLayer - 1; iCurLayer >= 0; iCurLayer--)
	{
		// get scale between current layer and higher (previous) layer
		iLayerScale = piLayerScales[iCurLayer]; 

		iCumulativeScale *= iLayerScale; // scale between current layer and top layer

		// set matching parameters
		h = (int)(h_default*sqrt((double)(iCumulativeScale)) + 0.5);	// enlarge matching window (half window)
		max_h = (int)(max_h_default*sqrt((double)(iCumulativeScale)) + 0.5);
		w = 2 * h + 1;
		ww = w*w;

		diff_threshold = (iCurLayer > 0) ? 0.35 : 0.4; // empiracal dislocation thresholds

		// get (pyramid) image names
		GetPyramidLayerInfo(pLSReg, true, iCurLayer, pacImageFileName1, &iColMax, &iRowMax);
		GetPyramidLayerInfo(pLSReg, false, iCurLayer, pacImageFileName2, &iColMax, &iRowMax);

		// input image 1
		GetFileDir(pLSReg->pacImageFilePathName1, pacDataDir);
		if (iCurLayer > 0)
			sprintf_s(pacImageFileName, STRLEN, "%s/pyramids/%s", pacDataDir, pacImageFileName1);
		else
			sprintf_s(pacImageFileName, STRLEN, "%s/%s", pacDataDir, pacImageFileName1);
		pshtImg1 = (short int*)ReadInputImage(pacImageFileName, iColMax, iRowMax, 1);

		// input image 2
		GetFileDir(pLSReg->pacImageFilePathName2, pacDataDir);
		if (iCurLayer > 0)
			sprintf_s(pacImageFileName, STRLEN, "%s/pyramids/%s", pacDataDir, pacImageFileName2);
		else
			sprintf_s(pacImageFileName, STRLEN, "%s/%s", pacDataDir, pacImageFileName2);
		pshtImg2 = (short int*)ReadInputImage(pacImageFileName, iColMax, iRowMax, 1);

		nrow = iRowMax;
		ncol = iColMax;

		// do matching on current layer
		for (n = 0; n < ntie; n++)
		{
			if (pucTiesFlag[n] == 1)
				continue;	// point n have been detected as mismatch

			x2 = piPreCols2[n] * iCumulativeScale;
			y2 = piPreRows2[n] * iCumulativeScale;

			x1 = pfPreCols1[n] * iLayerScale;	// note: pdPreCols1 and pdPreRows1 are the matching coordinates on the higher layer and are updated with new matchings
			y1 = pfPreRows1[n] * iLayerScale;

			// skip possible water pixels; added 1/27/2017
			if (pshtImg2[y2*iColMax + x2] <= value_threshold || pshtImg1[(int)(y1 + 1e-10)*iColMax + (int)(x1 + 1e-10)] <= value_threshold)
			{
				pucTiesFlag[n] = 1;
				ntie_new -= 1;
				continue;
			}

			// skip possible cloud pixels for MSS
			if (pshtImg2[y2*iColMax + x2] >= SATURATED_VALUE_MSS || pshtImg1[(int)(y1 + 1e-10)*iColMax + (int)(x1 + 1e-10)] >= SATURATED_VALUE_MSS)
			{
				pucTiesFlag[n] = 1;
				ntie_new -= 1;
				continue;
			}

			// (x1, y1) is initially matched to (x2n, y2n), i.e. (x2, y2)
			// do least-squares matching (LSM) to get new matched position (x1n, y1n)
			corr = -1;
			x2n = (float)(x2);
			y2n = (float)(y2);
			xdif = 0.f;
			ydif = 0.f;
			ws = w;
			do
			{
				// get LSM matched position (x1n, y1n)
				x1n = x1;
				y1n = y1;
				LSMatching_SAM(pshtImg2, ncol, nrow, pshtImg1, ncol, nrow, ws, ws, x2n, y2n, &x1n, &y1n, &corr, 32767);

				// increment matching window size
				ws = ws + 4 * (int)(sqrt((double)(iCumulativeScale)) + 0.5);

				// compare (x1n, y1n) with initial position (x1, y1)
				xdif = ABS(x1n - x1);
				ydif = ABS(y1n - y1);
			} while (xdif < 1e-10 && ydif < 1e-10 && ws / 2 < max_h && corr > 0);

			pfNewCols1[n*(iLayer + 1) + (iLayer - iCurLayer)] = x1n;
			pfNewRows1[n*(iLayer + 1) + (iLayer - iCurLayer)] = y1n;
			pfNewCoefs[n*(iLayer + 1) + (iLayer - iCurLayer)] = corr;

			if (((xdif<1e-10 && ydif<1e-10) || sqrt(xdif*xdif + ydif*ydif)>diff_threshold || corr < corr_threshold) && pucTiesFlag[n] == 0)
			{
				// mark point n as a mismatch
				pucTiesFlag[n] = 1;
				ntie_new -= 1;
			}

			// update pdPreCols1 and pdPreRows1, i.e. LSM matched positions
			pfPreCols1[n] = x1n;
			pfPreRows1[n] = y1n;
		}
		free(pshtImg2);
		pshtImg2 = NULL;
		free(pshtImg1);
		pshtImg1 = NULL;
	}

	// output matched points on bottom layer (10m resolution)
	sprintf_s(pacImageFileName, STRLEN, "%s/matching/poi_layer%d_matching_DF.txt", pacWorkDir, iLayer);
	fout = Writetxt(pacImageFileName);
	iCumulativeScale = 1;
	for (iCurLayer = iLayer; iCurLayer >= 0; iCurLayer--)
	{
		if (iCurLayer < iLayer)
			iCumulativeScale *= piLayerScales[iCurLayer];
	}
	iCurLayer = 0;
	for (n = 0; n < ntie; n++)
	{
		if (pucTiesFlag[n] == 1)
			continue;

		x2 = piPreCols2[n] * iCumulativeScale;
		y2 = piPreRows2[n] * iCumulativeScale;
		x1 = pfNewCols1[n*(iLayer + 1) + (iLayer - iCurLayer)];
		y1 = pfNewRows1[n*(iLayer + 1) + (iLayer - iCurLayer)];
		xdif = x1 - pfNewCols1[n*(iLayer + 1) + (iLayer - iCurLayer - 1)] * piLayerScales[iCurLayer];
		ydif = y1 - pfNewRows1[n*(iLayer + 1) + (iLayer - iCurLayer - 1)] * piLayerScales[iCurLayer];
		corr = (pfNewCoefs[n*(iLayer + 1) + (iLayer - iCurLayer)]);

		fprintf(fout, "%d\t%d\t%d\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%.3f\n", n + 1, x2, y2, x1, y1, x1 - x2, y1 - y2, sqrt(xdif*xdif + ydif*ydif), corr);
	}
	fclose(fout);
	fout = NULL;

	printf(" %d tie points pass depth-first matching.\n", ntie_new);

	// save transformation fitting results
	ntie = GetRegistrationTransformation_DepthFirst(pLSReg, pLSReg->bIsImage2BaseImage, enum_AUTO, pdCoefs, &dRMSE);
	memmove(pLSReg->pdCoefs, pdCoefs, MAX_COEFs_NUM * sizeof(double));
	pLSReg->n_tie = ntie;
	pLSReg->dRMSE = dRMSE;
	//printf(" RMSE:\t%.3f pixels\n", dRMSE);

	free(piPreRows2);
	free(piPreCols2);
	free(pfPreRows1);
	free(pfPreCols1);
	free(pfPreCoefs);
	free(pfNewRows1);
	free(pfNewCols1);
	free(pfNewCoefs);
	free(pucTiesFlag);

	return ntie_new;
}

/***
Fit tranformation f() between image 1 and image 2 using tie points for image registration (xa, ya) = f(xb, yb)
Tie points are from depth-first matching results ASCII file (*matching_DF.txt).

Return number of tie points after removing outliers

Parameters
bIsImage2BaseImage:			= false to register image 2 (e.g. Landsat-8) to image 1 (e.g. Sentinel-2); the transformation is (x2, y2) = f(x1, y1);
= true to register image 1 to image 2; the transformation is (x1, y1) = f(x2, y2)
enum_TransformationType:	= enum_TRANSLATION, enum_AFFINE, or enum_POLYNOMIAL
pdCoefs[12]:				transformation coefficients
***/
int GetRegistrationTransformation_DepthFirst(LSDepthFirstRegistration *pLSReg, bool bIsImage2BaseImage, enum_tranformation_type enum_TransformationType_input, double *pdCoefs, double *pdRMSE)
{
	char pacWorkDir[STRLEN] = "";
	char pacImageFileName[STRLEN] = "";

	FILE *fin = NULL;
	double tmp_value, corr, diff;
	int ntie, n, ntie_outliers;
	double *pdRows0 = NULL, *pdCols0 = NULL, *pdRows1 = NULL, *pdCols1 = NULL;
	int *piTieIndexs = NULL;
	double *pdFittingResiduals = NULL;
	int iLayer;
	double dMeanResidual, dFittingRMSE;
	double *pdRowsBase = NULL, *pdColsBase = NULL, *pdRowsTarget = NULL, *pdColsTarget = NULL;
	enum_tranformation_type enum_TransformationType;

	memset(pdCoefs, 0, MAX_COEFs_NUM * sizeof(double));

	iLayer = pLSReg->iLayerNum;
	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);

	// get tie points number
	sprintf_s(pacImageFileName, STRLEN, "%s/matching/poi_layer%d_matching_DF.txt", pacWorkDir, iLayer);
	ntie = Gettxtlines(pacImageFileName);
	if (ntie < 2)
	{
		// 		printf("\n Error in GetRegistrationTransformation_DepthFirst(): %d matched tie points are insufficient.\n", ntie); 
		return 0; // added 1/31/2017
	}

	// allocate memory
	if (!(pdRows0 = (double*)calloc(ntie, sizeof(double)))
		|| !(pdCols0 = (double*)calloc(ntie, sizeof(double)))
		|| !(pdRows1 = (double*)calloc(ntie, sizeof(double)))
		|| !(pdCols1 = (double*)calloc(ntie, sizeof(double)))
		|| !(piTieIndexs = (int*)calloc(ntie, sizeof(int)))
		|| !(pdFittingResiduals = (double*)calloc(ntie, sizeof(double))))
	{
		printf("\nError in GetRegistrationTransformation_DepthFirst(): insufficient memory.\n");
		scanf_s(" %d", &n);
		exit(1);
	}

	// input tie points coordinates
	sprintf_s(pacImageFileName, STRLEN, "%s/matching/poi_layer%d_matching_DF.txt", pacWorkDir, iLayer);
	fin = Readtxt(pacImageFileName);
	for (n = 0; n < ntie; n++)
		fscanf_s(fin, "%d%lf%lf%lf%lf%lf%lf%lf%lf", piTieIndexs + n, pdCols0 + n, pdRows0 + n, pdCols1 + n, pdRows1 + n, &tmp_value, &tmp_value, &diff, &corr);
	fclose(fin);

	// set base image and target image
	// note parallax are obtained by matching image 1 to image 2.
	if (bIsImage2BaseImage)
	{
		pdRowsBase = pdRows0;
		pdColsBase = pdCols0;
		pdRowsTarget = pdRows1;
		pdColsTarget = pdCols1;
	}
	else
	{
		pdRowsBase = pdRows1;
		pdColsBase = pdCols1;
		pdRowsTarget = pdRows0;
		pdColsTarget = pdCols0; 
	}

	enum_TransformationType = enum_TransformationType_input;
	if (enum_TransformationType_input == enum_AUTO)
		enum_TransformationType = AutomaticTransformationType(ntie);

	// fit transformation
	switch (enum_TransformationType)
	{
	case enum_TRANSLATION:
		FitTranslationTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_AFFINE:
		FitAffineTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_POLYNOMIAL:
		FitPolynomialTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	default:
		FitTranslationTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	}

	// remove outliers
	ntie_outliers = 0;
	for (n = 0; n < ntie; n++)
	{
		if (pdFittingResiduals[n] == -1)
			break;

		if (pdFittingResiduals[n] > 2 * dFittingRMSE)
		{
			memmove(pdCols1 + n, pdCols1 + n + 1, (ntie - n - 1 - ntie_outliers)*sizeof(double));
			memmove(pdRows1 + n, pdRows1 + n + 1, (ntie - n - 1 - ntie_outliers)*sizeof(double));
			memmove(pdCols0 + n, pdCols0 + n + 1, (ntie - n - 1 - ntie_outliers)*sizeof(double));
			memmove(pdRows0 + n, pdRows0 + n + 1, (ntie - n - 1 - ntie_outliers)*sizeof(double));
			memmove(pdFittingResiduals + n, pdFittingResiduals + n + 1, (ntie - n - 1 - ntie_outliers)*sizeof(double));
			pdCols1[ntie - ntie_outliers - 1] = -1;
			pdRows1[ntie - ntie_outliers - 1] = -1;
			pdCols0[ntie - ntie_outliers - 1] = -1;
			pdRows0[ntie - ntie_outliers - 1] = -1;
			pdFittingResiduals[ntie - ntie_outliers - 1] = -1;

			ntie_outliers += 1;
			n -= 1;
		}
	}

	// refit transformation
	ntie -= ntie_outliers;

	enum_TransformationType = enum_TransformationType_input;
	if (enum_TransformationType_input == enum_AUTO)
		enum_TransformationType = AutomaticTransformationType(ntie);

	switch (enum_TransformationType)
	{
	case enum_TRANSLATION:
		FitTranslationTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_AFFINE:
		FitAffineTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_POLYNOMIAL:
		FitPolynomialTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	default:
		FitTranslationTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	}

	free(pdRows0);
	free(pdCols0);
	free(pdRows1);
	free(pdCols1);
	free(piTieIndexs);
	free(pdFittingResiduals);

	*pdRMSE = dFittingRMSE;

	return ntie;
}

/***
Fit tranformation f() between image 1 and image 2 using tie points for image registration (xa, ya) = f(xb, yb)
v1: Tie points are from dense grid point matching results binary file

Return number of tie points after removing outliers

Parameters
bIsImage2BaseImage:         = false to register image 2 (e.g. Landsat-8) to image 1 (e.g. Sentinel-2); the transformation is (x2, y2) = f(x1, y1);
= true to register image 1 to image 2; the transformation is (x1, y1) = f(x2, y2)
enum_TransformationType:    = enum_TRANSLATION, enum_AFFINE, or enum_POLYNOMIAL
pacDenseMatchingFileName:   dense matching file (binary)
pdCoefs[12]:                fitted transformation coefficients
9/3/2017: fixed a memory leak
1/22/2020: utilize piTieIndexs (previously declared but not used) to handle outliers removal
3/10/2020: update pLSReg->iRBFs_K if FitRBFsTransform_TPS_Poly_v2() is used
***/
int GetRegistrationTransformation_DenseMatching(LSDepthFirstRegistration* pLSReg, bool bIsImage2BaseImage, enum_tranformation_type enum_TransformationType, int iLayer, int iSamplingInterval, char* pacDenseMatchingFileName, double* pdCoefs, double* pdRMSE)
{
	char pacWorkDir[STRLEN] = "";
	char pacImageFileName[STRLEN] = "";

	FILE* fin = NULL;
	int ntie, n, ntie_outliers, ntie_good;
	double* pdRows0 = NULL, * pdCols0 = NULL, * pdRows1 = NULL, * pdCols1 = NULL;
	int* piTieIndexs = NULL;	//
	double* pdFittingResiduals = NULL;
	double dMeanResidual, dFittingRMSE;
	double* pdRowsBase = NULL, * pdColsBase = NULL, * pdRowsTarget = NULL, * pdColsTarget = NULL;

	float* pfParalaxMap_x = NULL, * pfParalaxMap_y = NULL;
	float* pfCorrMap = NULL;
	int Row, Col;
	int iRowMax, iColMax;
	int iMatchingRowMax, iMatchingColMax;

	// added 1/6/2020
	int iRBFs_K = pLSReg->iRBFs_K;
	int iMin_ties_RBF;
	iMin_ties_RBF = 12 + 2 * iRBFs_K;

	int iRBFs_K_new = 0; // added 3/10/2020

	memset(pdCoefs, 0, MAX_COEFs_NUM * sizeof(double));

	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);
	GetPyramidLayerInfo(pLSReg, true, iLayer, pacImageFileName, &iColMax, &iRowMax);

	// input matched tie points coordinates
	iMatchingRowMax = iRowMax / iSamplingInterval;
	iMatchingColMax = iColMax / iSamplingInterval;
	pfParalaxMap_x = (float*)ReadInputImage(pacDenseMatchingFileName, iMatchingRowMax, iMatchingColMax * 4, 3);
	pfParalaxMap_y = pfParalaxMap_x + iMatchingRowMax * iMatchingColMax;
	pfCorrMap = pfParalaxMap_x + 3 * iMatchingRowMax * iMatchingColMax;

	// get number of matched points
	ntie = 0;
	for (n = 0; n < iMatchingRowMax * iMatchingColMax; n++)
	{
		if (pfCorrMap[n] > 0.1f)
			ntie += 1;
	}

	// 9/3/22017: moved up; originall below the following memory call
	if ((enum_TransformationType == enum_POLYNOMIAL && ntie < MIN_TIES_NUM_POLYNOMIAL)
		|| (enum_TransformationType == enum_POLYNOMIAL && ntie < iMin_ties_RBF)
		|| (enum_TransformationType != enum_POLYNOMIAL && ntie < MIN_TIES_NUM))
	{
		free(pfParalaxMap_x);	// added 9/3/2017
		return ntie;
	}

	if (!(pdRows0 = (double*)calloc(ntie, sizeof(double)))
		|| !(pdCols0 = (double*)calloc(ntie, sizeof(double)))
		|| !(pdRows1 = (double*)calloc(ntie, sizeof(double)))
		|| !(pdCols1 = (double*)calloc(ntie, sizeof(double)))
		|| !(piTieIndexs = (int*)calloc(ntie, sizeof(int)))
		|| !(pdFittingResiduals = (double*)calloc(ntie, sizeof(double))))
	{
		printf("\nError in GetRegistrationTransformation_DenseMatching(): insufficient memory.\n");
		scanf_s(" %d", &n);
		exit(1);
	}

	// get coordinates of matched points
	ntie = 0;
	for (Row = 0; Row < iMatchingRowMax; Row++)
	{
		for (Col = 0; Col < iMatchingColMax; Col++)
		{
			n = Row * iMatchingColMax + Col;
			if (pfCorrMap[n] > 0.1f)
			{
				pdRows0[ntie] = Row * iSamplingInterval;
				pdCols0[ntie] = Col * iSamplingInterval;
				pdRows1[ntie] = pdRows0[ntie] + pfParalaxMap_y[n];
				pdCols1[ntie] = pdCols0[ntie] + pfParalaxMap_x[n];
				ntie += 1;
			}
		}
	}

	// set base image and target image
	if (bIsImage2BaseImage)
	{
		pdRowsBase = pdRows0;
		pdColsBase = pdCols0;
		pdRowsTarget = pdRows1;
		pdColsTarget = pdCols1;
	}
	else
	{
		pdRowsBase = pdRows1;
		pdColsBase = pdCols1;
		pdRowsTarget = pdRows0;
		pdColsTarget = pdCols0;
	}

	// fit transformation
	switch (enum_TransformationType)
	{
	case enum_TRANSLATION:
		FitTranslationTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_AFFINE:
		FitAffineTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_POLYNOMIAL:
		FitPolynomialTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_RBFs:
		iRBFs_K_new = FitRBFsTransform_TPS_Poly_v2(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, iRBFs_K, iColMax, iRowMax, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		pLSReg->iRBFs_K_adp = iRBFs_K_new;
		break;
	default:
		FitAffineTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	}

	// remove outliers
	ntie_outliers = 0;
	ntie_good = 0; // added 1/22/2020
	for (n = 0; n < ntie; n++)
		piTieIndexs[n] = (int)(-1);
	for (n = 0; n < ntie; n++)
	{
		if (pdFittingResiduals[n] <= 2 * dFittingRMSE)
		{
			piTieIndexs[ntie_good] = n;
			ntie_good += 1;
		}
		else
			ntie_outliers += 1;
	}

	for (n = 0; n < ntie; n++)
	{
		if (piTieIndexs[n] >= 0)
		{
			pdCols1[n] = pdCols1[piTieIndexs[n]];
			pdRows1[n] = pdRows1[piTieIndexs[n]];
			pdCols0[n] = pdCols0[piTieIndexs[n]];
			pdRows0[n] = pdRows0[piTieIndexs[n]];
			pdFittingResiduals[n] = pdFittingResiduals[piTieIndexs[n]];
		}
		else
		{
			pdCols1[n] = -1;
			pdRows1[n] = -1;
			pdCols0[n] = -1;
			pdRows0[n] = -1;
			pdFittingResiduals[n] = (int)(-1);
		}
	}

	// refit transformation
	ntie -= ntie_outliers;
	if ((enum_TransformationType == enum_POLYNOMIAL && ntie < MIN_TIES_NUM_POLYNOMIAL)
		|| (enum_TransformationType == enum_POLYNOMIAL && ntie < iMin_ties_RBF)
		|| (enum_TransformationType != enum_POLYNOMIAL && ntie < MIN_TIES_NUM))
		return ntie;

	switch (enum_TransformationType)
	{
	case enum_TRANSLATION:
		FitTranslationTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_AFFINE:
		FitAffineTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_POLYNOMIAL:
		FitPolynomialTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	case enum_RBFs:
		iRBFs_K_new = FitRBFsTransform_TPS_Poly_v2(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, iRBFs_K, iColMax, iRowMax, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		pLSReg->iRBFs_K_adp = iRBFs_K_new;
		break;
	default:
		FitAffineTransform(pdColsBase, pdRowsBase, pdColsTarget, pdRowsTarget, ntie, pdCoefs, pdFittingResiduals, &dMeanResidual, &dFittingRMSE);
		break;
	}

	free(pdRows0);
	free(pdCols0);
	free(pdRows1);
	free(pdCols1);
	free(piTieIndexs);
	free(pdFittingResiduals);
	free(pfParalaxMap_x);

	*pdRMSE = dFittingRMSE;

	return ntie;
}

/***
Created 1/7/2020. Revised from GetRegistrationTransformation_DenseMatching()
FitTransform = true, fit transformation from v11 dense matching file and save in pLSReg->pdCoefs_DenseMatching
FitTransform = false, use transformation coefficients fitted in DenseMatching() (on v11 file)

then input dense matching file v11_ext to calculate the residual map
***/
void OutputFittingResidualMaps(LSDepthFirstRegistration* pLSReg, bool FitTransform)
{
	char pacImageFileName[STRLEN] = "";
	FILE* fout = NULL;
	double pdCoefs[MAX_COEFs_NUM];
	double dRMSE;
	int ntie;

	char pacImageFileName_[STRLEN] = "";

	enum_tranformation_type enum_TransformationType = pLSReg->enum_TransformationType;
	int iLayer = 0;
	int iSamplingInterval = pLSReg->iDenseMatchingGridSize;
	char pacDenseMatchingFileName[STRLEN];

	char pacWorkDir[STRLEN] = "";

	FILE* fin = NULL;
	int n;

	float* pfParalaxMap_x = NULL, * pfParalaxMap_y = NULL;
	float* pfCorrMap = NULL;
	int Row, Col;
	int iRowMax, iColMax;
	int iMatchingRowMax, iMatchingColMax;

	int iRBFs_K = pLSReg->iRBFs_K_adp;
	int iMin_ties_RBF = 12 + 2 * iRBFs_K;

	float* pfResidualsMap = NULL;

	int iTransformationType;
	double dTargetCol, dTargetRow;
	double dMatchedCol, dMatchedRow;
	float fRow0, fCol0;
	double dRes_x, dRes_y, dRes;
	char pacDescription[STRLEN] = "";
	char pacBandNames[STRLEN] = "";
	int ntie_fitting;

	iTransformationType = pLSReg->enum_TransformationType;

	// fit designated transformation
	if (FitTransform == true)
	{
		sprintf_s(pacImageFileName_, STRLEN, "%s/matching/layer%d_dense_matching_rsp%d.v11", pLSReg->pacWorkDir, 0, pLSReg->iDenseMatchingGridSize);
		sprintf_s(pacDenseMatchingFileName, STRLEN, "%s", pacImageFileName_);
		ntie_fitting = GetRegistrationTransformation_DenseMatching(pLSReg, pLSReg->bIsImage2BaseImage, pLSReg->enum_TransformationType, 0, pLSReg->iDenseMatchingGridSize, pacImageFileName_, pdCoefs, &dRMSE);
		memmove(pLSReg->pdCoefs_DenseMatching, pdCoefs, MAX_COEFs_NUM * sizeof(double));
		pLSReg->dRMSE_DenseMatching = dRMSE;
		pLSReg->n_tie_DenseMatching = ntie_fitting;
	}
	else
	{
		// get fitting coefficients obtained in DenseMatching()
		memmove(pdCoefs, pLSReg->pdCoefs_DenseMatching, MAX_COEFs_NUM * sizeof(double));
		dRMSE = pLSReg->dRMSE_DenseMatching;
		ntie_fitting = pLSReg->n_tie_DenseMatching;
	}

	sprintf_s(pacImageFileName_, STRLEN, "%s/matching/layer%d_dense_matching_rsp%d.v11_ext", pLSReg->pacWorkDir, 0, pLSReg->iDenseMatchingGridSize);
	sprintf_s(pacDenseMatchingFileName, STRLEN, "%s", pacImageFileName_);

	sprintf_s(pacImageFileName_, STRLEN, "%s/RBFs_stats/layer%d_dense_matching_rsp%d_residuals.v11_K-%d_ext", pLSReg->pacWorkDir, 0, pLSReg->iDenseMatchingGridSize, iRBFs_K);
	fout = WriteBinary(pacImageFileName_);

	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);
	GetPyramidLayerInfo(pLSReg, true, iLayer, pacImageFileName, &iColMax, &iRowMax);

	// input matched tie points coordinates
	iMatchingRowMax = iRowMax / iSamplingInterval;
	iMatchingColMax = iColMax / iSamplingInterval;
	pfParalaxMap_x = (float*)ReadInputImage(pacDenseMatchingFileName, iMatchingRowMax, iMatchingColMax * 4, 3);
	pfParalaxMap_y = pfParalaxMap_x + iMatchingRowMax * iMatchingColMax;
	pfCorrMap = pfParalaxMap_x + 3 * iMatchingRowMax * iMatchingColMax;

	// get number of matched points
	ntie = 0;
	for (n = 0; n < iMatchingRowMax * iMatchingColMax; n++)
	{
		if (pfCorrMap[n] > 0.1f)
			ntie += 1;
	}

	// 9/3/22017: moved up; originall below the following memory call
	if ((enum_TransformationType == enum_POLYNOMIAL && ntie < MIN_TIES_NUM_POLYNOMIAL)
		|| (enum_TransformationType == enum_POLYNOMIAL && ntie < iMin_ties_RBF)
		|| (enum_TransformationType != enum_POLYNOMIAL && ntie < MIN_TIES_NUM))
	{
		free(pfParalaxMap_x);	// added 9/3/2017
		return;
	}

	if (!(pfResidualsMap = (float*)calloc(iMatchingRowMax*iMatchingColMax, sizeof(float))))
	{
		printf("\nError in OutputFittingResidualMaps(): insufficient memory.\n");
		scanf_s(" %d", &n);
		exit(1);
	}

	// get coordinates of matched points
	ntie = 0;
	for (Row = 0; Row < iMatchingRowMax; Row++)
	{
		for (Col = 0; Col < iMatchingColMax; Col++)
		{
			n = Row * iMatchingColMax + Col;
			if (pfCorrMap[n] > 0.1f)
			{
				fRow0 = (float)(Row) * iSamplingInterval;
				fCol0 = (float)(Col) * iSamplingInterval;

				dMatchedCol = (double)(fCol0 + pfParalaxMap_x[n]);
				dMatchedRow = (double)(fRow0 + pfParalaxMap_y[n]);

				// note parallax are obtained by matching image 1 to image 2.
				GetTransformedCoords((double)dMatchedCol, (double)dMatchedRow, iTransformationType, pdCoefs, &dTargetCol, &dTargetRow, pLSReg->iRBFs_K_adp);
				
				dRes_x = dTargetCol - fCol0;
				dRes_y = dTargetRow - fRow0;
				dRes = (float)(sqrt(dRes_x * dRes_x + dRes_y * dRes_y));

				if (dRes > 1.0 && pfCorrMap[n] < 0.995f) // 1/24/2020, change dRes from 0.5 to 1.0
					continue;

				pfResidualsMap[n] = (float)dRes;
			}
		}
	}

	fwrite(pfResidualsMap, sizeof(float), iMatchingColMax* iMatchingRowMax, fout);
	fclose(fout);

	sprintf_s(pacDescription, STRLEN, "%d tie points. RMSE: %.3f pixels.", ntie_fitting, dRMSE);

	// output ENVI hdr file
	OutputEnviHDRFile(pacImageFileName_, pacDescription, iMatchingColMax, iMatchingRowMax, 1, 4, pacBandNames);
	free(pfParalaxMap_x);
	free(pfResidualsMap);

	return;
}

/***
4/21/2021: change from 3 times to 1.5 times
***/
enum_tranformation_type AutomaticTransformationType(int ntie)
{
	if (ntie >= 18)
		return enum_POLYNOMIAL;
	else if (ntie >= 9)
		return enum_AFFINE;
	else // (ntie > 6)
		return enum_TRANSLATION;
}

/***
Fit the three types of transformations and output fitted transformation coefficients in file poi_layer3_matching_DF_coefs_v2.txt in intermediate-results folder

note:
The output coefficients represent the transformations from the base image to the target image (pLSReg->pacTargetImageFilePathName), i.e. (x_target, y_target) = f(x_base, y_base).
They are to be used to register the target image to the base image.
***/
void OutputTransformations_v2(LSDepthFirstRegistration *pLSReg)
{
	char pacWorkDir[STRLEN] = "", pacImageFileName[STRLEN] = "";
	FILE *fout = NULL;
	double pdCoefs[MAX_COEFs_NUM];
	double dRMSE;
	int ntie;

	char pacImageFileName_[STRLEN] = "";

	// output transformation coefficients
	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);
	sprintf_s(pacImageFileName, STRLEN, "%s/results/poi_layer%d_matching_DF_coefs_v2.txt", pacWorkDir, pLSReg->iLayerNum);
	fout = Writetxt(pacImageFileName);

	// image pair information
	fprintf(fout, "Image 1: %s\n", pLSReg->pacImageFilePathName1);
	fprintf(fout, "Image 2: %s\n", pLSReg->pacImageFilePathName2);

	// v2
	sprintf_s(pacImageFileName_, STRLEN, "%s/matching/layer%d_dense_matching_rsp%d", pacWorkDir, 0, pLSReg->iDenseMatchingGridSize);

	// fit translation transformation
	ntie = GetRegistrationTransformation_DenseMatching(pLSReg, pLSReg->bIsImage2BaseImage, enum_TRANSLATION, 0, pLSReg->iDenseMatchingGridSize, pacImageFileName_, pdCoefs, &dRMSE);
	fprintf(fout, "\nTranslation transformation (dense matching)\n");
	fprintf(fout, "a0:\t%.9f\nb0:\t%.9f\n", pdCoefs[0], pdCoefs[1]);
	fprintf(fout, "RMSE:\t%.3f\n", dRMSE);
	fprintf(fout, "n:\t%d\n", ntie);

	//	printf(" Translation transformation (dense matching).");
	//	printf("\tRMSE: %.3f pixels.", dRMSE);
	//	printf("\tn: %d\n", ntie);

	// fit affine transformation
	ntie = GetRegistrationTransformation_DenseMatching(pLSReg, pLSReg->bIsImage2BaseImage, enum_AFFINE, 0, pLSReg->iDenseMatchingGridSize, pacImageFileName_, pdCoefs, &dRMSE);
	fprintf(fout, "\nAffine transformation (dense matching)\n");
	fprintf(fout, "a0:\t%.9f\na1:\t%.9f\na2:\t%.9f\nb0:\t%.9f\nb1:\t%.9f\nb2:\t%.9f\n", pdCoefs[0], pdCoefs[1], pdCoefs[2], pdCoefs[3], pdCoefs[4], pdCoefs[5]);
	fprintf(fout, "RMSE:\t%.3f\n", dRMSE);
	fprintf(fout, "n:\t%d\n", ntie);

	//	printf(" Affine transformation (dense matching).");
	//	printf("\tRMSE: %.3f pixels.", dRMSE);
	//	printf("\tn: %d\n", ntie);
	//	printf(" a0 (x shift): %.9f   b0 (y shift): %.9f\n", pdCoefs[0], pdCoefs[3]);

	// fit polynomial transformation
	ntie = GetRegistrationTransformation_DenseMatching(pLSReg, pLSReg->bIsImage2BaseImage, enum_POLYNOMIAL, 0, pLSReg->iDenseMatchingGridSize, pacImageFileName_, pdCoefs, &dRMSE);
	fprintf(fout, "\nPolynomial transformation (dense matching)\n");
	fprintf(fout, "a0:\t%.9f\na1:\t%.9f\na2:\t%.9f\na3:\t%.9f\na4:\t%.9f\na5:\t%.9f\nb0:\t%.9f\nb1:\t%.9f\nb2:\t%.9f\nb3:\t%.9f\nb4:\t%.9f\nb5:\t%.9f\n", pdCoefs[0], pdCoefs[1], pdCoefs[2], pdCoefs[3], pdCoefs[4], pdCoefs[5], pdCoefs[6], pdCoefs[7], pdCoefs[8], pdCoefs[9], pdCoefs[10], pdCoefs[11]);
	fprintf(fout, "RMSE:\t%.3f\n", dRMSE);
	fprintf(fout, "n:\t%d\n", ntie);
	//	printf(" Polynomial transformation (dense matching).");
	//	printf("\tRMSE: %.3f pixels.", dRMSE);
	//	printf("\tn: %d\n", ntie);
	//	printf(" a0 (x shift): %.9f   b0 (y shift): %.9f\n", pdCoefs[0], pdCoefs[6]);

	fclose(fout);

	return;
}

/***
formal registration output file: [image 1 name]-[image 2 name]_coefs.txt in pacWorkDir

Output file contents:
x shift, y shift
transformation type
transformation fitting RMSE
transformation fitting tie points number
transformation coefficients

Example:
-0.136119748 0.477286985
2
0.094
770
-0.102495755	0.999999644	-0.000006559	0.482854684	0.000000906	0.999997698

note:
x shift, y shift are obtained by fitting translation transformation coefficients a0, b0
***/
void OutputRegistrationNumaricals(LSDepthFirstRegistration *pLSReg, double *ptr_a0_translation, double *ptr_b0_translation, double *ptr_RMSE, int *ptr_ntie)
{
	char pacImageFileName[STRLEN] = "";
	FILE *fout = NULL;
	double pdCoefs[MAX_COEFs_NUM];
	double dRMSE;
	int ntie;

	char pacImageFileName_[STRLEN] = "";
	char pacUpperDir[STRLEN] = "";

	*ptr_a0_translation = 99999;
	*ptr_b0_translation = 99999;
	*ptr_RMSE = 99999;
	*ptr_ntie = 0;

	// output transformation coefficients
//	GetUpperDir(pLSReg->pacWorkDir, pacUpperDir);
	sprintf_s(pacImageFileName, STRLEN, "%s/results/%s-%s_coefs.txt", pLSReg->pacWorkDir, pLSReg->pacImageFileName1, pLSReg->pacImageFileName2);
	fout = Writetxt(pacImageFileName)	;

	sprintf_s(pacImageFileName_, STRLEN, "%s/matching/layer%d_dense_matching_rsp%d", pLSReg->pacWorkDir, 0, pLSReg->iDenseMatchingGridSize);

	// fit translation transformation, output fitted a0 and b0
	ntie = GetRegistrationTransformation_DenseMatching(pLSReg, pLSReg->bIsImage2BaseImage, enum_TRANSLATION, 0, pLSReg->iDenseMatchingGridSize, pacImageFileName_, pdCoefs, &dRMSE);
	fprintf(fout, "%.9f %.9f\n", pdCoefs[0], pdCoefs[1]);
	*ptr_a0_translation = pdCoefs[0];
	*ptr_b0_translation = pdCoefs[1];

	// fit designated transformation
	ntie = GetRegistrationTransformation_DenseMatching(pLSReg, pLSReg->bIsImage2BaseImage, pLSReg->enum_TransformationType, 0, pLSReg->iDenseMatchingGridSize, pacImageFileName_, pdCoefs, &dRMSE);
	fprintf(fout, "%d\n", pLSReg->enum_TransformationType);
	fprintf(fout, "%.3f\n", dRMSE);
	fprintf(fout, "%d\n", ntie);
	switch (pLSReg->enum_TransformationType)
	{
	case 1:
		fprintf(fout, "%.9f\t%.9f\n", pdCoefs[0], pdCoefs[1]);
		break;
	case 2:
		fprintf(fout, "%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n", pdCoefs[0], pdCoefs[1], pdCoefs[2], pdCoefs[3], pdCoefs[4], pdCoefs[5]);
		break;
	case 3:
		fprintf(fout, "%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n", pdCoefs[0], pdCoefs[1], pdCoefs[2], pdCoefs[3], pdCoefs[4], pdCoefs[5], pdCoefs[6], pdCoefs[7], pdCoefs[8], pdCoefs[9], pdCoefs[10], pdCoefs[11]);
		break;
	default:
		fprintf(fout, "%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n", pdCoefs[0], pdCoefs[1], pdCoefs[2], pdCoefs[3], pdCoefs[4], pdCoefs[5], pdCoefs[6], pdCoefs[7], pdCoefs[8], pdCoefs[9], pdCoefs[10], pdCoefs[11]);
		break;
	}
	fclose(fout);

	*ptr_RMSE = dRMSE;
	*ptr_ntie = ntie;

	return;
}


/***
Register target image (pLSReg->pacTargetImageFilePathName) to base image; output registered binary image to intermediate-results folder
***/
void GenerateRegisteredImage_v1_1(LSDepthFirstRegistration *pLSReg)
{
	int iRowMax, iColMax;
	char pacImageFileName[STRLEN] = "";
	char pacImageFilePathName[STRLEN] = "";

	FILE *fin = NULL, *fout = NULL;
	double pdCoefs[MAX_COEFs_NUM], dRMSE;
	int ntie;
	int Row, Col;
	short int *pshtImg0 = NULL, *pshtImg0_new = NULL;

	double dTargetRow, dTargetCol;
	int iTargetRow, iTargetCol;

	float fTargetValue;
	float fWeight, fWeightSum;
	float fDX, fDY;
	int iTransformationType;

	char pacDescription[STRLEN] = "";
	char pacTargetImageFileName[STRLEN] = "";
	bool bUseDFCoefs = false;  // alwayes use densing matching tie points rather than depth-first matching tie points
	char pacResamplingMethod[STRLEN] = "NNBL"; // NN and BL
	char pac_null[2] = "";

	int iRBFs_K = pLSReg->iRBFs_K_adp;

	iTransformationType = pLSReg->enum_TransformationType;

	if (strlen(pLSReg->pacTargetImageFilePathName) == 0)
	{
		printf("No output registered image.");
		return;
	}

	GetFileName(pLSReg->pacTargetImageFilePathName, pacTargetImageFileName);

	pacResamplingMethod[pLSReg->enum_ResamplingMethod * 2] = '\0';

	// get output file name, "BL" stands for bilinear
	switch (pLSReg->enum_TransformationType)
	{
	case enum_TRANSLATION:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_translation_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_AFFINE:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_affine_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_POLYNOMIAL:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_polynomial_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_RBFs:
		sprintf_s(pacImageFileName, STRLEN, "%s.v11_reg_RBFs_poly_K-%d_%s", pacTargetImageFileName, iRBFs_K, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;	
	default:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_affine_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	}

	sprintf_s(pacImageFilePathName, STRLEN, "%s/results/%s", pLSReg->pacWorkDir, pacImageFileName);
	// check whether output file exists
	/*	if (FileExist(pacImageFilePathName))
	{
	printf("\nFile exists. Registered binary image output as\n  %s\n", pacImageFileName);
	return;
	}*/

	iRowMax = pLSReg->iRowMax;
	iColMax = pLSReg->iColMax;

	if (bUseDFCoefs)
	{
		memmove(pdCoefs, pLSReg->pdCoefs, MAX_COEFs_NUM * sizeof(double));
		ntie = pLSReg->n_tie;
		dRMSE = pLSReg->dRMSE;
	}
	else
	{
		memmove(pdCoefs, pLSReg->pdCoefs_DenseMatching, MAX_COEFs_NUM * sizeof(double));
		ntie = pLSReg->n_tie_DenseMatching;
		dRMSE = pLSReg->dRMSE_DenseMatching;
	}
	// input target image to be registered to base image
	pshtImg0 = (short int*)ReadInputImage(pLSReg->pacTargetImageFilePathName, iColMax, iRowMax, 1);

	// allocate memory for new registered (warped) image
	if (!(pshtImg0_new = (short int*)calloc(iRowMax*iColMax, sizeof(short int))))
	{
		printf("\nError in GenerateRegisteredImage(): insufficient memory.\n");
		scanf_s(" %d", &Row);
		exit(1);
	}

	//	printf("\nGenerating registered image\n");
	memset(pshtImg0_new, 0, iRowMax*iColMax*sizeof(short int));
	for (Row = 0; Row < iRowMax; Row++)
	{
		for (Col = 0; Col < iColMax; Col++)
		{
			GetTransformedCoords((double)Col, (double)Row, iTransformationType, pdCoefs, &dTargetCol, &dTargetRow, pLSReg->iRBFs_K_adp);

			if (pLSReg->enum_ResamplingMethod == 1)
			{
				iTargetCol = (int)(dTargetCol + 0.5 + 1e-10);
				iTargetRow = (int)(dTargetRow + 0.5 + 1e-10);

				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
					pshtImg0_new[Row*iColMax + Col] = pshtImg0[iTargetRow*iColMax + iTargetCol];
			}
			else // pLSReg->enum_ResamplingMethod == 2
			{
				// for bilinear resampling
				fWeightSum = 0;
				fTargetValue = 0;

				// top left point
				iTargetCol = (int)dTargetCol;
				iTargetRow = (int)dTargetRow;

				fDX = (float)(dTargetCol - iTargetCol);
				fDY = (float)(dTargetRow - iTargetRow);
				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
				{
					fWeight = (1 - fDX)*(1 - fDY);
					fTargetValue += pshtImg0[iTargetRow*iColMax + iTargetCol] * fWeight;
					fWeightSum += fWeight;
				}

				// top right point
				iTargetCol += 1;
				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
				{
					fWeight = fDX*(1 - fDY);
					fTargetValue += pshtImg0[iTargetRow*iColMax + iTargetCol] * fWeight;
					fWeightSum += fWeight;
				}

				// bottom right point
				iTargetRow += 1;
				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
				{
					fWeight = fDX*fDY;
					fTargetValue += pshtImg0[iTargetRow*iColMax + iTargetCol] * fWeight;
					fWeightSum += fWeight;
				}

				// bottom left point
				iTargetCol -= 1;
				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
				{
					fWeight = (1 - fDX)*fDY;
					fTargetValue += pshtImg0[iTargetRow*iColMax + iTargetCol] * fWeight;
					fWeightSum += fWeight;
				}

				if (fWeightSum > 0)
				{
					fTargetValue /= fWeightSum;
					pshtImg0_new[Row*iColMax + Col] = (short)(fTargetValue + 0.5f);
				}
			}
		}
	}

	// output
	fout = WriteBinary(pacImageFilePathName);
	fwrite(pshtImg0_new, sizeof(short int), iRowMax*iColMax, fout);
	fclose(fout);

	printf(" Registered binary image output in intermediate-results folder: %s\n", pacImageFileName);

	// get description text
	if (pLSReg->enum_TransformationType == enum_TRANSLATION)
		sprintf_s(pacDescription, STRLEN, "Registered image by translation transformation. %d tie points. RMSE: %.3f pixels.", ntie, dRMSE);
	else if (pLSReg->enum_TransformationType == enum_AFFINE)
		sprintf_s(pacDescription, STRLEN, "Registered image by affine transformation. %d tie points. RMSE: %.3f pixels.", ntie, dRMSE);
	else
		sprintf_s(pacDescription, STRLEN, "Registered image by polynomial transformation. %d tie points. RMSE: %.3f pixels.", ntie, dRMSE);

	// output ENVI hdr file
	OutputEnviHDRFile(pacImageFilePathName, pacDescription, iColMax, iRowMax, 1, 2, pac_null);

	// release memory
	free(pshtImg0);
	free(pshtImg0_new);

	return;
}

/***
Register target image (pLSReg->pacTargetImageFilePathName) to base image; output registered binary image to intermediate-results folder
***/
void GenerateRegisteredImage_v1_1_UC(LSDepthFirstRegistration* pLSReg)
{
	int iRowMax, iColMax;
	char pacImageFileName[STRLEN] = "";
	char pacImageFilePathName[STRLEN] = "";

	FILE* fin = NULL, * fout = NULL;
	double pdCoefs[MAX_COEFs_NUM], dRMSE;
	int ntie;
	int Row, Col;
	unsigned char* pucImg0 = NULL, * pucImg0_new = NULL;

	double dTargetRow, dTargetCol;
	int iTargetRow, iTargetCol;

	float fTargetValue;
	float fWeight, fWeightSum;
	float fDX, fDY;
	int iTransformationType;

	char pacDescription[STRLEN] = "";
	char pacTargetImageFileName[STRLEN] = "";
	bool bUseDFCoefs = false;  // alwayes use densing matching tie points rather than depth-first matching tie points
	char pacResamplingMethod[STRLEN] = "NNBL"; // NN and BL
	char pac_null[2] = "";

	int iRBFs_K = pLSReg->iRBFs_K;

	iTransformationType = pLSReg->enum_TransformationType;

	if (strlen(pLSReg->pacTargetImageFilePathName) == 0)
	{
		printf("No output registered image.");
		return;
	}

	GetFileName(pLSReg->pacTargetImageFilePathName, pacTargetImageFileName);

	pacResamplingMethod[pLSReg->enum_ResamplingMethod * 2] = '\0';

	// get output file name, "BL" stands for bilinear
	switch (pLSReg->enum_TransformationType)
	{
	case enum_TRANSLATION:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_translation_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_AFFINE:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_affine_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_POLYNOMIAL:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_polynomial_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_RBFs:
		sprintf_s(pacImageFileName, STRLEN, "%s.v11_reg_RBFs_poly_K-%d_%s", pacTargetImageFileName, iRBFs_K, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	default:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_affine_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	}

	sprintf_s(pacImageFilePathName, STRLEN, "%s/results/%s", pLSReg->pacWorkDir, pacImageFileName);
	// check whether output file exists
	/*	if (FileExist(pacImageFilePathName))
	{
	printf("\nFile exists. Registered binary image output as\n  %s\n", pacImageFileName);
	return;
	}*/

	iRowMax = pLSReg->iRowMax;
	iColMax = pLSReg->iColMax;

	if (bUseDFCoefs)
	{
		memmove(pdCoefs, pLSReg->pdCoefs, MAX_COEFs_NUM * sizeof(double));
		ntie = pLSReg->n_tie;
		dRMSE = pLSReg->dRMSE;
	}
	else
	{
		memmove(pdCoefs, pLSReg->pdCoefs_DenseMatching, MAX_COEFs_NUM * sizeof(double));
		ntie = pLSReg->n_tie_DenseMatching;
		dRMSE = pLSReg->dRMSE_DenseMatching;
	}
	// input target image to be registered to base image
	pucImg0 = (unsigned char*)ReadInputImage(pLSReg->pacTargetImageFilePathName, iColMax, iRowMax, 0);

	// allocate memory for new registered (warped) image
	if (!(pucImg0_new = (unsigned char*)calloc(iRowMax * iColMax, sizeof(unsigned char))))
	{
		printf("\nError in GenerateRegisteredImage(): insufficient memory.\n");
		scanf_s(" %d", &Row);
		exit(1);
	}

	//	printf("\nGenerating registered image\n");
	memset(pucImg0_new, 0, iRowMax * iColMax * sizeof(unsigned char));
	for (Row = 0; Row < iRowMax; Row++)
	{
		for (Col = 0; Col < iColMax; Col++)
		{
			GetTransformedCoords((double)Col, (double)Row, iTransformationType, pdCoefs, &dTargetCol, &dTargetRow, pLSReg->iRBFs_K_adp);

			if (pLSReg->enum_ResamplingMethod == 1)
			{
				iTargetCol = (int)(dTargetCol + 0.5 + 1e-10);
				iTargetRow = (int)(dTargetRow + 0.5 + 1e-10);

				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
					pucImg0_new[Row * iColMax + Col] = pucImg0[iTargetRow * iColMax + iTargetCol];
			}
			else // pLSReg->enum_ResamplingMethod == 2
			{
				// for bilinear resampling
				fWeightSum = 0;
				fTargetValue = 0;

				// top left point
				iTargetCol = (int)dTargetCol;
				iTargetRow = (int)dTargetRow;

				fDX = (float)(dTargetCol - iTargetCol);
				fDY = (float)(dTargetRow - iTargetRow);
				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
				{
					fWeight = (1 - fDX) * (1 - fDY);
					fTargetValue += pucImg0[iTargetRow * iColMax + iTargetCol] * fWeight;
					fWeightSum += fWeight;
				}

				// top right point
				iTargetCol += 1;
				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
				{
					fWeight = fDX * (1 - fDY);
					fTargetValue += pucImg0[iTargetRow * iColMax + iTargetCol] * fWeight;
					fWeightSum += fWeight;
				}

				// bottom right point
				iTargetRow += 1;
				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
				{
					fWeight = fDX * fDY;
					fTargetValue += pucImg0[iTargetRow * iColMax + iTargetCol] * fWeight;
					fWeightSum += fWeight;
				}

				// bottom left point
				iTargetCol -= 1;
				if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
				{
					fWeight = (1 - fDX) * fDY;
					fTargetValue += pucImg0[iTargetRow * iColMax + iTargetCol] * fWeight;
					fWeightSum += fWeight;
				}

				if (fWeightSum > 0)
				{
					fTargetValue /= fWeightSum;
					pucImg0_new[Row * iColMax + Col] = (unsigned char)(fTargetValue + 0.5f);
				}
			}
		}
	}

	// output
	fout = WriteBinary(pacImageFilePathName);
	fwrite(pucImg0_new, sizeof(unsigned char), iRowMax * iColMax, fout);
	fclose(fout);

	printf(" Registered binary image output in intermediate-results folder: %s\n", pacImageFileName);

	// get description text
	if (pLSReg->enum_TransformationType == enum_TRANSLATION)
		sprintf_s(pacDescription, STRLEN, "Registered image by translation transformation. %d tie points. RMSE: %.3f pixels.", ntie, dRMSE);
	else if (pLSReg->enum_TransformationType == enum_AFFINE)
		sprintf_s(pacDescription, STRLEN, "Registered image by affine transformation. %d tie points. RMSE: %.3f pixels.", ntie, dRMSE);
	else
		sprintf_s(pacDescription, STRLEN, "Registered image by polynomial transformation. %d tie points. RMSE: %.3f pixels.", ntie, dRMSE);

	// output ENVI hdr file
	OutputEnviHDRFile(pacImageFilePathName, pacDescription, iColMax, iRowMax, 1, 1, pac_null);

	// release memory
	free(pucImg0);
	free(pucImg0_new);

	return;
}

void GetRegisteredImageFileName(LSDepthFirstRegistration* pLSReg, char* pacImageFilePathName)
{
	char pacImageFileName[STRLEN] = "";
	char pacTargetImageFileName[STRLEN] = "";
	int iTransformationType;

	char pacResamplingMethod[STRLEN] = "NNBL"; // NN and BL
	char pac_null[2] = "";

	int iRBFs_K = pLSReg->iRBFs_K;

	iTransformationType = pLSReg->enum_TransformationType;

	if (strlen(pLSReg->pacTargetImageFilePathName) == 0)
	{
		printf("No output registered image.");
		return;
	}

	GetFileName(pLSReg->pacTargetImageFilePathName, pacTargetImageFileName);

	pacResamplingMethod[pLSReg->enum_ResamplingMethod * 2] = '\0';

	// get output file name, "BL" stands for bilinear
	switch (pLSReg->enum_TransformationType)
	{
	case enum_TRANSLATION:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_translation_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_AFFINE:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_affine_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_POLYNOMIAL:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_polynomial_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	case enum_RBFs:
		sprintf_s(pacImageFileName, STRLEN, "%s.v11_reg_RBFs_poly_K-%d_%s", pacTargetImageFileName, iRBFs_K, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	default:
		sprintf_s(pacImageFileName, STRLEN, "%s.reg_affine_%s", pacTargetImageFileName, pacResamplingMethod + (pLSReg->enum_ResamplingMethod - 1) * 2);
		break;
	}

	sprintf_s(pacImageFilePathName, STRLEN, "%s/results/%s", pLSReg->pacWorkDir, pacImageFileName);

	return;
}

/***
Dense grid point matching. For a sampled grid point on image 2, get least squares matched position on image 1 to get the x shift (x2-x1) and y shift (y2-y1)

Parameters
iLayer:				layer on which dense matching is undertaken; usually the bottom layer with highest resolution
iSamplingInterval:	grid size by pixels
iHalfWindow:		half matching window size; if not set (default = 0), the window size used in depth-first matching will be used

Note
1. For all LSM matching, image 1 is held still and image 2 is resampled (warped) to fit image 1
2. Two thresholds (shtMeanDiffThreshold and shtBandValueThreshold) are used. shtMeanDiffThreshold is used to exclude matching on cloudy patches.
shtBandValueThreshold is used to excluded water. See notes for LSMatching_SAM().
3. shtBandValueThreshold is used assuming the input image1 and image 2 are NIR images.
4. The dense matching results are theoretically independent on the transformation type, which is automatically selected.
1/9/2020: output pfParalaxMap as the 5th layer
1/27/2020: add parameter bExtensive; default = false and output .v11_ext; = true to apply LSMatching_SAM_v1() and output .v11
1/28/2020: add corr_threshold_ext and diff_threshold_ext
1/28/2020: make corr_threshold_ext adaptive, [0.985, 0.995] <-> [2*dRMSE, 3*dRMSE]
_UC: created 7/26/2021
_UC_v1 (7/28/2021)
-no bExtensive parameter
-add iMatchingMode prameter
	if = 0, normal depth-first-based dense matching
	if = 1, do not use depth-first matching results, just do pixel-to-pixel matching with no shifts; add "_1" in output file name
	if = 2, same as 1 except add "_2" in output file name
-use LSMatching_SAM_UC_precise() and scale-distort criteria
-use shrinking matching windows instead of enlarging matching windows
7/22/2022: added (dif<= 4 && corr>=0.995) condition
v2 (7/27/2022): used for MSS L1TP dense matching
- do pixel-by-pixel dense matching with initial window matching
- remove parameter 'iMatchingMode'
- add parameter iHalfSearchWindow
- no iterations with decreased match window
***/
int DenseMatching_UC_v2(LSDepthFirstRegistration* pLSReg, int iLayer, int iSamplingInterval, int iHalfWindow, int iHalfSearchWindow)
{
	char pacWorkDir[STRLEN] = "", pacImageFileName1[STRLEN] = "", pacImageFileName2[STRLEN] = "";
	int iRowMax, iColMax;
	char pacImageFileName[STRLEN] = "";

	unsigned char* pucImg2 = NULL, * pucImg1 = NULL;
	FILE* fin = NULL, * fout = NULL;
	int iCumulativeScale;

	// variables for image matching
	float corr, corr_real, corr_real_org;
	float x1n, y1n, x2n, y2n, x1d, y1d;
	double xdif, ydif, dif;
	int ws, w, ww;
	int h;
	int h_default;
	float x1, y1;
	int ncol, nrow;
	double diff_threshold;
	double corr_threshold;
	unsigned char value_threshold;
	unsigned char ucMeanDiffThreshold;		// do not do matching if the mean difference of two matching windows is larger than the threshold
	unsigned char ucBandValueThreshold;	// do not do matching if band value is smaller than this threshold; used to filter out water; a value of 500 is used for NIR band

	ucMeanDiffThreshold = pLSReg->ucMeanDiffThreshold;
	ucBandValueThreshold = pLSReg->ucBandValueThreshold;

	int Row, Col;
	int Row_new, Col_new;
	float* pfParalaxMap_x = NULL, * pfParalaxMap_y = NULL, * pfParalaxMap = NULL;
	float* pfCorrMap = NULL, * pfPredictErrorMap = NULL;
	int iRowMax_new, iColMax_new;
	int iMatchedNum, iPointNum;

	int iTransformationType;

	char pacDescription[STRLEN] = "", pacImageFileName_[STRLEN] = "", pacBandNames[STRLEN] = "";

	int ntie_DenseMatching;
	double pdCoefs_DenseMatching[MAX_COEFs_NUM], dRMSE_DenseMatching;
	int ntie;

	char pacDataDir[STRLEN] = "";

	float pfLSM_paras[8];
	float* pfLSM_parasMap = NULL;
	int i;

	bool bSuc; // v1
	float fScale_distort, fShear_distort;
	int iIter, iMaxIter;

	// v2
	int x, y, x2, y2;
	int dx, dy;
	double cor, maxcor;
	unsigned char* sub1 = NULL, * sub2 = NULL;

	dx = iHalfSearchWindow;
	dy = iHalfSearchWindow;

	if (iSamplingInterval == 0)
		return 0;

	printf(" Dense matching. ");

	// get output file name
	iCumulativeScale = GetCumulativeScaleFromBottom(pLSReg, pLSReg->iLayerNum) / GetCumulativeScaleFromBottom(pLSReg, iLayer); // scale between current layer and top layer
	if (iHalfWindow == 0)
	{
		// get half window size used in depth-first matching
		h_default = pLSReg->iDefaultHalfWindow;
		h = (int)(h_default * sqrt((double)(iCumulativeScale)) + 0.5);
	}
	else
		h = iHalfWindow; // use input value
	w = 2 * h + 1;	// matching window size
	sprintf_s(pacWorkDir, STRLEN, "%s", pLSReg->pacWorkDir);
	sprintf_s(pacImageFileName_, STRLEN, "%s/matching/layer%d_dense_matching_rsp%d", pacWorkDir, iLayer, iSamplingInterval);

	// check whether output file exists
	if (FileExist(pacImageFileName_)) // disabled for MSS
	{
		//		printf("File exists. Dense matching results output to\n  %s\n", pacImageFileName_);
		printf("File exists.\n");

		// save transformation fitting results
		ntie_DenseMatching = GetRegistrationTransformation_DenseMatching(pLSReg, pLSReg->bIsImage2BaseImage, pLSReg->enum_TransformationType, iLayer, iSamplingInterval, pacImageFileName_, pdCoefs_DenseMatching, &dRMSE_DenseMatching);
		memmove(pLSReg->pdCoefs_DenseMatching, pdCoefs_DenseMatching, MAX_COEFs_NUM * sizeof(double));
		pLSReg->n_tie_DenseMatching = ntie_DenseMatching;
		pLSReg->dRMSE_DenseMatching = dRMSE_DenseMatching;

		return ntie_DenseMatching;
	}

	ntie = 0;
	iTransformationType = 0;
	
	// pixel-to-pixel matching with no shifts
	corr_threshold = 0.96;	
	diff_threshold = 2;
	
	// get (pyramid) image names
	GetPyramidLayerInfo(pLSReg, true, iLayer, pacImageFileName1, &iColMax, &iRowMax);
	GetPyramidLayerInfo(pLSReg, false, iLayer, pacImageFileName2, &iColMax, &iRowMax);

	// input image 1
	GetFileDir(pLSReg->pacImageFilePathName1, pacDataDir);
	if (iLayer > 0)
		sprintf_s(pacImageFileName, STRLEN, "%s/pyramids/%s", pacDataDir, pacImageFileName1);
	else
		sprintf_s(pacImageFileName, STRLEN, "%s", pLSReg->pacImageFilePathName1);
	pucImg1 = (unsigned char*)ReadInputImage(pacImageFileName, iColMax, iRowMax, 0);

	// input image 2
	GetFileDir(pLSReg->pacImageFilePathName2, pacDataDir);
	if (iLayer > 0)
		sprintf_s(pacImageFileName, STRLEN, "%s/pyramids/%s", pacWorkDir, pacImageFileName2);
	else
		sprintf_s(pacImageFileName, STRLEN, "%s", pLSReg->pacImageFilePathName2);
	pucImg2 = (unsigned char*)ReadInputImage(pacImageFileName, iColMax, iRowMax, 0);

	iRowMax_new = iRowMax / iSamplingInterval;
	iColMax_new = iColMax / iSamplingInterval;

	// set matching parameters
	iCumulativeScale = GetCumulativeScaleFromBottom(pLSReg, pLSReg->iLayerNum) / GetCumulativeScaleFromBottom(pLSReg, iLayer); // scale between current layer and top layer
	if (iHalfWindow == 0)
	{
		// get half window size used in depth-first matching
		h_default = pLSReg->iDefaultHalfWindow;
		h = (int)(h_default * sqrt((double)(iCumulativeScale)) + 0.5);	// same as the window size on this layer in depth-first matching
	}
	else
		h = iHalfWindow; // use input value
	w = 2 * h + 1;
	ww = w * w;

	value_threshold = ucBandValueThreshold;

	nrow = iRowMax;
	ncol = iColMax;

	if (!(pfParalaxMap_x = (float*)calloc(iRowMax_new * iColMax_new, sizeof(float)))
		|| !(pfParalaxMap_y = (float*)calloc(iRowMax_new * iColMax_new, sizeof(float)))
		|| !(pfParalaxMap = (float*)calloc(iRowMax_new * iColMax_new, sizeof(float)))
		|| !(pfCorrMap = (float*)calloc(iRowMax_new * iColMax_new, sizeof(float)))
		|| !(pfLSM_parasMap = (float*)calloc(iRowMax_new * iColMax_new * 4, sizeof(float)))
		|| !(pfPredictErrorMap = (float*)calloc(iRowMax_new * iColMax_new, sizeof(float)))
		|| !(sub1 = (unsigned char*)calloc(ww, sizeof(unsigned char)))
		|| !(sub2 = (unsigned char*)calloc(ww, sizeof(unsigned char))))
	{
		printf("\nError in DenseMatching(): insufficient memory.\n");
		scanf_s(" %d", &Row);
		exit(1);
	}

	iMatchedNum = 0;
	iPointNum = 0;
	//	printf(" (progress,matching ratio)\n");
	for (Row_new = 0; Row_new < iRowMax_new; Row_new++)
	{
		if ((Row_new + 1) % MAX((10 * iRowMax_new / 100), 1) == 0 && iPointNum > 0)
			printf(" (%d%%,%.1f%%)", (int)((Row_new + 1) * 100.0 / iRowMax_new) + 1, iMatchedNum * 1.0f / iPointNum * 100); // the second number is sucessful matching ratio; if this number is low, it means the matching criteria parameters are set too high, or one or two images are highly cloudy

		for (Col_new = 0; Col_new < iColMax_new; Col_new++)
		{
			// get a sampled grid point (Row, Col) on image 2
			Row = Row_new * iSamplingInterval;
			Col = Col_new * iSamplingInterval;

			// skip possible water pixels
			if (pucImg2[Row * iColMax + Col] <= value_threshold || pucImg1[Row * iColMax + Col] <= value_threshold)
				continue;

			// skip possible cloud pixels for MSS
			if (pucImg2[Row * iColMax + Col] >= SATURATED_VALUE_MSS_DENSEMATCHING || pucImg1[Row * iColMax + Col] >= SATURATED_VALUE_MSS_DENSEMATCHING)
				continue;

			// skip if matching window contain fill values
			if (FindTargetValueInWindow_UC(pucImg2, ncol, nrow, Col, Row, h, FILLVALUE))
				continue;

			// v2
			x2 = Col;
			y2 = Row;
			if (!imsub_UC(pucImg2, ncol, nrow, x2, y2, h, sub2))
				continue;

			maxcor = -1; // 1/9/2020: changed from 0 to -1 because corr2() return value can be negative
			for (y = y2 - dy; y <= y2 + dy; y++)
			{
				for (x = x2 - dx; x <= x2 + dx; x++)
				{
					if (imsub_UC(pucImg1, ncol, nrow, x, y, h, sub1))
					{
						cor = corr2_UC(sub2, sub1, ww);
						if (cor > maxcor)
						{
							x1 = (float)x;
							y1 = (float)y;
							maxcor = cor;
						}
					}
				}
			}
			if (maxcor < 0)
				continue; // added 7/27/2021

			iPointNum += 1; // base for matching ratio calculation (water pixels and fill value pixels excluded)

			// (x1, y1) is initially matched to (x2n, y2n), i.e. (x2, y2)
			// least-square matching
			corr = -1;
			x2n = (float)(Col);
			y2n = (float)(Row);
			xdif = 0.f;
			ydif = 0.f;
			ws = w;
			iIter = 0;
			iMaxIter = 3;
			do
			{
				iIter += iMaxIter;

				x1n = x1;
				y1n = y1;

				LSMatching_SAM_UC_precise_fast(pucImg2, ncol, nrow, pucImg1, ncol, nrow, ws, ws, x2n, y2n, &x1n, &y1n, &corr, &corr_real, &corr_real_org, ucMeanDiffThreshold, pfLSM_paras, 1, 20, 0.01); // iMinIter = 1, iMaxIter = 20, h1_threshold = 0.01

				x1d = x1n;
				y1d = y1n;

				xdif = ABS(x1d - x1);
				ydif = ABS(y1d - y1);
				dif = sqrt(xdif * xdif + ydif * ydif);

				bSuc = Check_LSM_components(pfLSM_paras, corr, corr_real, corr_real_org, &fScale_distort, &fShear_distort, true);
				// pixel-to-pixel matching with no shifts
				if (bSuc == true && corr >= corr_threshold && ( (dif > 1e-10 && dif <= 1.51)|| (dif <= 2.51 && corr >= 0.995) ) )
					break;
				else
					bSuc = false;

				ws = ws - 8;

			} while (xdif < 1e-10 && ydif < 1e-10 && corr > 0 && iIter < iMaxIter); // corr condition added 3/3/2016
			if (bSuc == true)
			{
				pfParalaxMap_x[Row_new * iColMax_new + Col_new] = x1n - x2n;
				pfParalaxMap_y[Row_new * iColMax_new + Col_new] = y1n - y2n;
				pfParalaxMap[Row_new * iColMax_new + Col_new] = sqrt((x1n - x2n) * (x1n - x2n) + (y1n - y2n) * (y1n - y2n));
				pfCorrMap[Row_new * iColMax_new + Col_new] = corr;
				pfPredictErrorMap[Row_new * iColMax_new + Col_new] = (float)(dif);
				for (i = 0; i < 2; i++)
					pfLSM_parasMap[iColMax_new * iRowMax_new * i + Row_new * iColMax_new + Col_new] = pfLSM_paras[i];
				pfLSM_parasMap[iColMax_new * iRowMax_new * 2 + Row_new * iColMax_new + Col_new] = fScale_distort;
				pfLSM_parasMap[iColMax_new * iRowMax_new * 3 + Row_new * iColMax_new + Col_new] = fShear_distort;

				iMatchedNum += 1;
			}
		}
	}

	// output
	fout = WriteBinary(pacImageFileName_);
	fwrite(pfParalaxMap_x, sizeof(float), iRowMax_new * iColMax_new, fout);
	fwrite(pfParalaxMap_y, sizeof(float), iRowMax_new * iColMax_new, fout);
	fwrite(pfPredictErrorMap, sizeof(float), iRowMax_new * iColMax_new, fout);
	fwrite(pfCorrMap, sizeof(float), iRowMax_new * iColMax_new, fout);
	fwrite(pfParalaxMap, sizeof(float), iRowMax_new * iColMax_new, fout);
	fwrite(pfLSM_parasMap, sizeof(float), iRowMax_new * iColMax_new * 4, fout);
	fclose(fout);

	// get description text
	if (pLSReg->enum_TransformationType == enum_TRANSLATION)
		sprintf_s(pacDescription, STRLEN, "Dense matching results (translation transformation used). Grid size: %d pixels. %d grid points considered. %.2f%% matched.", iSamplingInterval, iPointNum, iMatchedNum * 1.0f / iPointNum * 100);
	else if (pLSReg->enum_TransformationType == enum_AFFINE)
		sprintf_s(pacDescription, STRLEN, "Dense matching results (affine transformation used). Grid size: %d pixels. %d grid points considered. %.2f%% matched.", iSamplingInterval, iPointNum, iMatchedNum * 1.0f / iPointNum * 100);
	else
		sprintf_s(pacDescription, STRLEN, "Dense matching results (polynomial transformation used). Grid size: %d pixels. %d grid points considered. %.2f%% matched.", iSamplingInterval, iPointNum, iMatchedNum * 1.0f / iPointNum * 100);

	// write band names
	sprintf_s(pacBandNames, STRLEN, "x shift, y shift, prediction error, SAM, shift, h0, h1, scale distort, shear distort");
	// output ENVI hdr file
	OutputEnviHDRFile(pacImageFileName_, pacDescription, iColMax_new, iRowMax_new, 9, 4, pacBandNames);

	// save transformation fitting results
	ntie_DenseMatching = GetRegistrationTransformation_DenseMatching(pLSReg, pLSReg->bIsImage2BaseImage, pLSReg->enum_TransformationType, iLayer, iSamplingInterval, pacImageFileName_, pdCoefs_DenseMatching, &dRMSE_DenseMatching);
	memmove(pLSReg->pdCoefs_DenseMatching, pdCoefs_DenseMatching, MAX_COEFs_NUM * sizeof(double));
	pLSReg->n_tie_DenseMatching = ntie_DenseMatching;
	pLSReg->dRMSE_DenseMatching = dRMSE_DenseMatching;

	// release memory
	free(sub1);
	free(sub2);
	free(pucImg2);
	free(pucImg1);
	free(pfParalaxMap_x);
	free(pfParalaxMap_y);
	free(pfParalaxMap);
	free(pfCorrMap);
	free(pfPredictErrorMap);
	free(pfLSM_parasMap);

	//	printf("\n Dense matching results output to\n  %s\n", pacImageFileName);
	printf(" %d (%.2f%%) grid points matched\n", ntie_DenseMatching, iMatchedNum * 1.0f / iPointNum * 100);

	return ntie_DenseMatching;
}
