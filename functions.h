#ifndef FUNCTIONS_H
#define FUNCTIONS_H
			
#include <math.h>
#include <stdio.h>				
#include <stdlib.h>
#include <string.h>
#include <errno.h>				
#include <Float.h>
#include <direct.h>
#include <iostream>
#include <fstream>
#include <stdint.h>

#define DELTA_LIMIT			(0.001f)
#define FILLVALUE			(0)
#define FILLVALUE_1			(1)
#define FILLVALUE_255			(255)
#define MIN_MEAN_PATCH_VALUE	(5)
#define ISFILLVALUE(x)		((ABS(x)==0)? 1 : 0)

#define ABS(x)				(((x)>=0)? (x):-(x))
#define MAX(a,b) (((a)>(b))? (a) : (b))
#define MIN(a,b) (((a)<(b))? (a) : (b))
#define IS_OVERLAPPLING(x1_min, x1_max, x2_min, x2_max) (MAX((x1_min), (x2_min)) <= MIN((x1_max), (x2_max))? 1 : 0)
#define EQUAL_FLT(x1, x2) ((ABS((x1)-(x2))<1e-6)? (true):(false))

#define STRLEN (200)
#define STRLEN_LONG (2000)

#define S2A_L1C_ROWS			10980 // band 8 (NIR)
#define S2A_L1C_COLS			10980
#define S2A_L1C_RESOLUTION		10

#define MAX_SENTINEL2A_SHIFT_THRESHOLD (6)		// experience value, used to filter out (rare) bad matchings; currently the max observed shift is over 5 pixels found 
												//  between S2A_OPER_PRD_MSIL1C_PDMC_20160829T222904_R050_V20160829T090552_20160829T093120.SAFE\GRANULE\S2A_OPER_MSI_L1C_TL_SGS__20160829T143818_A006196_T33LTD_N02.04 
												//      and S2A_OPER_PRD_MSIL1C_PDMC_20160816T213752_R007_V20160816T090022_20160816T091559.SAFE\GRANULE\S2A_OPER_MSI_L1C_TL_SGS__20160816T143027_A006010_T33LTD_N02.04

#define MAX_RBFs_K				(100)
#define DEFAULT_RBFs_K			(64)
#define MAX_COEFs_NUM			(12+MAX_RBFs_K*4)
#define RBF_TIE_NUMBER_THRESHOLD (50)

/// new functions for GCP chips matching, added 11/3/2020
int LSMatching_SAM_UC_precise_fast(unsigned char* imagel, int coll, int rowl, unsigned char* imager, int colr, int rowr, int img_coll, int img_rowl, float x1, float y1, float* ptr_x2, float* ptr_y2, float* ptr_cr, float* ptr_cr_real, float* ptr_cr_real_org, unsigned char ucMeanDiffThreshold, float* ptr_paras = NULL, int iMinIter = 0, int iMaxIter = 20, double h1_threshold = 0.5);
inline int resampling_UC_precise(unsigned char* image, int col, int row, float* img, int img_col, int img_row, int x1, int y1, double* par_x, double* par_y);
float GetMeanDiff_FLT(float* img1, int col1, int row1, float* img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height);
void SAM_FLT(float* img1, int col1, int row1, float* img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height, float* ptr_cr);
inline void correlation_coefficient_FLT_v2(float* img1, int col1, int row1, float* img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height, float* ptr_cr, float h0 = 0.0f, float h1 = 1.0f);
inline int imsub_UC(unsigned char* im, int ncol, int nrow, int x, int y, int h, unsigned char* sub);

inline bool normalized_FLT_simple(float* data, int n, float fFillValue = 0.0f, float fShift = 127, bool bMedian = false);
double corr2_UC(unsigned char* x, unsigned char* y, int n, unsigned char fillvalue = 0);

void GetDecomposedAffineComponents(float* pfAffineCoefs, float* pfScale_x, float* pfScale_y, float* pfShear);
bool Check_LSM_components(float* pfLSM_paras, float fSAM, float corr_real, float corr_real_org, float* pfScale_distort = NULL, float* pfShear_distort = NULL, bool bL1TP_matching = false);

int ConvertToDOY(int iYear, int iMonth, int iDay);

int GetConnectionGraphFromMatchingFile(char* pacMatchingFileName, int n_Images, bool* ptr_bConnectionMap, bool bIndexFrom_0, float fRMSE_threshold, float fShift_threshold, int ntie_threshold);
int FindUnconnectedImages(char* pacMatchingFileName, int n_Images, int iConnectNumThreshold, bool* ptr_bUnconnectedImagesFlag, bool bIndexFrom_0, float fRMSE_threshold, float fShift_threshold, int ntie_threshold);
int GetShifts_v2(char* pacDenseMatchingFileName, int iColMax, int iRowMax, int* pdCols, int* pdRows, float fSAMThreshold, float* ptr_fShifts_x, float* ptr_fShifts_y);
int Get_DenseMatching_tie_count_map_PerImage(char* pacMatchingSummaryFileName, int n_Images, int matching_width, int matching_height, int * piImages_Control_Flag,  bool bIndexFrom_0, float fRMSE_threshold, float fShift_threshold, int ntie_threshold, float fSAM_threshold, int* pi_tie_CountMapsPerImage);

// least squares matching
void LSMatching_SAM(short int* imagel, int coll, int rowl, short int* imager, int colr, int rowr, int img_coll, int img_rowl, float x1, float y1, float* ptr_x2, float* ptr_y2, float* ptr_cr, short int shtMeanDiffThreshold, float* ptr_paras = NULL, int iMinIter = 0, double h1_threshold = 0.5);
void SAM(short int *img1, int col1, int row1, short int *img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height, float *ptr_cr);
int resampling(short int *image, int col, int row, short int *img, int img_col, int img_row, int x1, int y1, double *par_x, double *par_y);
double corr2(short int *x, short int *y, int n);
int imsub(short int* im, int ncol, int nrow, int x, int y, int h, short int* sub);
short int GetMeanDiff(short int *img1, int col1, int row1, short int *img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height);

// statistics
double Mean1(double *x, int n);
void RMSE1(double *x, int n, int t, double *ptr_rmse, double *ptr_avg);
void Std1(double* x, int n, double* ptr_stdv, double* ptr_ave);
float Mean1_FLT(float* x, int n);
unsigned char Mean1_UC(unsigned char* x, int n);
void Std1_FLT(float* x, int n, float* ptr_stdv, float* ptr_ave);
inline int GetStatistics_FLT(float* pInputData, int iNum, float FillValue, float* pStats);
int LargerOrEqualTo_Num_INT(int *piData, int n, int iValue);
int LargerOrEqualTo_Num_FLT(float* pData, int n, float value);
float GetStd_UC(unsigned char* pucData, int n, unsigned char* pucMask);
void GetMinMax(double* pData, int n, bool bExludeZero, double* pMin, double* pMax);

// least squares fitting for tranformations
void FitTranslationTransform(double *x1, double *y1, double *x2, double *y2, int n, double Coefs[], double *ptr_Errors, double *ptr_errors_mean, double *ptr_fitting_rmse);
void FitAffineTransform(double *x1, double *y1, double *x2, double *y2, int n, double U[], double *ptr_Errors, double *ptr_errors_mean, double *ptr_fitting_rmse);
void FitPolynomialTransform(double *x1, double *y1, double *x2, double *y2, int n, double Coefs[], double *ptr_Errors, double *ptr_errors_mean, double *ptr_fitting_rmse);
void GetTransformedCoords(double x1, double y1, int iTransformationType, double *Coefs, double *ptr_x2, double *ptr_y2, int iRBFs_K = 0);

int FitRBFsTransform_TPS_Poly_v2(double* x1, double* y1, double* x2, double* y2, int n, int K, int iWidth, int iHeight, double Coefs[], double* ptr_Errors, double* ptr_errors_mean, double* ptr_fitting_rmse);
int GetGridded_RBF_Centers(int iWidth, int iHeight, int iRBFs_K, double* pd_x_k, double* pd_y_k, int* pi_grid_width, int* pi_grid_height);

void TransformImage_sht2char(short* pshtData, int iColMax, int iRowMax, int iTransformationType, double* pdCoefs, unsigned char* pucData_trsf);
void TransformImage_UC(unsigned char *pucData, int iColMax, int iRowMax, int iTransformationType, double* pdCoefs, unsigned char* pucData_trsf, int iRBFs_K = 0);
void CopySht2UChar(short* pusData, int iSize, unsigned char* pucData_cpy);

// solve normal equation in least squares adjustment
int INVSQR1(double* A, double* B, int n);
int INVSQR1_FLT(float* A, float* B, int n);
int INVSQR2(double* A, double* B, int n, bool* pbSolvedUnknownsFlags);

// POI detection
void CalcGradient2D(float* pImageTemp, float* pImageX, float* pImageY, int iWidth, int iHeight);
void CalcMultiply(float* pImageTemp1, float* pImageTemp2, float* pImageNew, int iSize);
void CalcDivide(float* pImageTemp1, float* pImageTemp2, float* pImageNew, int iSize);
void CalcAdd(float* pImageTemp1, float* pImageTemp2, float* pImageNew, int iSize);
void CalcSubtract(float* pImageTemp1, float* pImageTemp2, float* pImageNew, int iSize);
void CalcMinMaxMean(float* pImageTemp, float *pdMin, float *pdMax, float *pdMean, int iSize);
void CalcMinMaxMeanWithMask(float* pImg, unsigned char *pucMask, float *pdMin, float *pdMax, float *pdMean, int iSize);
void ZeroBuffer(float *pBuffer, int iSize);
void AddConst(float *pBuffer, float c, int iSize);

float* GetGaussian(double dSigma, int *iFilterWidth);
void Conv2same(short int* pImg, short int* pImgNew, int iWidth, int iHeight, float* dFilter, int w);
void Conv2same_UC(unsigned char* pImg, unsigned char* pImgNew, int iWidth, int iHeight, float* dFilter, int w);

// I/O
FILE* Readtxt(char* file);
FILE* Writetxt(char* file);
FILE* WriteBinary(char* file);
int Gettxtlines(char* file);
bool FileExist(char* file);
bool CheckFileExist(char* file);
void* ReadInputImage(char* pacImageFileName, int iColMax, int iRowMax, int iDataType);
void* ReadInputImage_v1(char* pacImageFileName, int iColMax, int iRowMax, int iDataType);
void* ReadInputImage_v2(char* pacImageFileName, int iColMax, int iRowMax, int iDataType, unsigned long ulOffsetElementsNum = 0);
void OutputEnviHDRFile(char *pacPathName, char *pacDescription, int iSamples, int iLines, int iBands, int iDataType, char *pacBandNames);
bool InputEnviHDRFile(char* pacPathName, int* ptr_iSamples, int* ptr_iLines, int* ptr_iBands, int* ptr_iDataType);

bool FindTargetValueInWindow_UC(unsigned char* pImg, int iWidth, int iHeight, int iTargetCol, int iTargetRow, int w, unsigned char targetValue);

bool GetFileDir(char *pacFilePathName, char *pacFileDir);
bool GetUpperDir(char *pacFilePathName, char *pacFileDir);
bool GetFileName(char *pacFilePathName, char *pacFileName);
bool GetFileName_NoExtension(char *pacFilePathName, char *pacFileName);

// pStats[4]: min, max, mean, median
inline int GetStatistics_FLT(float* pInputData, int iNum, float FillValue, float* pStats)
{
	int i, k, n;
	float Max, Min;
	float Value;
	float* pData = NULL;
	int iValidNum;
	double sum;

	if (iNum <= 0)
		return 0;

	if (!(pData = (float*)calloc(iNum, sizeof(float))))
	{
		printf("\nError in GetStatistics_FLT(): insufficient memory.\n");
		scanf_s(" %d", &i);
		return 0;
	}

	memset(pData, 0, iNum * sizeof(float));

	pStats[0] = 0;
	pStats[1] = 0;
	pStats[2] = 0;
	pStats[3] = 0;

	Max = 0;
	Min = FLT_MAX;
	sum = 0;
	iValidNum = 0;
	for (i = 0; i < iNum; i++)
	{
		Value = pInputData[i];

		if (ABS(Value - FillValue) < 1e-6)
			continue;

		if (Value > Max)
			Max = Value;

		if (Value < Min)
			Min = Value;

		sum += Value;
		pData[iValidNum] = Value;
		iValidNum += 1;
		for (k = 0; k < iValidNum - 1; k++)
		{
			if (Value < pData[k])
			{
				for (n = iValidNum - 1; n > k; n--)
					pData[n] = pData[n - 1];

				pData[k] = Value;
				break;
			}
		}
	}

	if (iValidNum > 0)
	{
		pStats[0] = Min;
		pStats[1] = Max;
		pStats[2] = (float)(sum / iValidNum);
		pStats[3] = pData[iValidNum / 2];
	}

	free(pData);

	return iValidNum;
}

/***
Bilinear resampling using geometric affine transformation parameteres
11/6/2020
-change
if (x20 < 0 || x20 >= col - 1 || y20 < 0 || y20 >= row - 1)
TO
if (x20 < 0 || x20 > col - 1 || y20 < 0 || y20 > row - 1)
_precise(12/1/2020): input and return float type rather than unsigned char
***/
inline int resampling_UC_precise(unsigned char* image, int col, int row, float* img, int img_col, int img_row, int x1, int y1, double* par_x, double* par_y)
{
	int i, j, i1, i2, j1, j2;
	float g, * img0 = NULL;
	unsigned char* image0 = NULL;
	double x2, y2, p, q;
	int		x20, y20;
	double   bit;

	img0 = img;
	i1 = -img_row / 2;
	i2 = img_row / 2;
	j1 = -img_col / 2;
	j2 = img_col / 2;
	for (i = i1; i <= i2; i++)
	{
		for (j = j1; j <= j2; j++)
		{
			x2 = par_x[0] + j * par_x[1] + i * par_x[2];		//	x2=x+a0+x*a1+y*a2
			y2 = par_y[0] + j * par_y[1] + i * par_y[2];		//	y2=y+b0+y*b1+x*b2
			x20 = (int)(x2 + 1e-10);
			y20 = (int)(y2 + 1e-10);

			if (x20 < 0 || x20 > col - 1 || y20 < 0 || y20 > row - 1)
				return 0;

			image0 = image + y20 * col + x20;
			p = x2 - x20;
			q = y2 - y20;
			g = 0;
			bit = 0;
			if ((1 - p) * (1 - q) > 1e-10) // top left (x20, y20)
				bit += (*image0) * (1 - p) * (1 - q);
			if (p * (1 - q) > 1e-10 && x20 + 1 <= col - 1) // top right (x20+1, y20)
				bit += (*(image0 + 1)) * p * (1 - q);
			if ((1 - p) * q > 1e-10 && y20 + 1 <= row - 1) // bottom left (x20, y20+1)
				bit += (*(image0 + col)) * (1 - p) * q;
			if (p * q > 1e-10 && x20 + 1 <= col - 1 && y20 + 1 <= row - 1) // bottom right (x20+1, y20+1)
				bit += (*(image0 + col + 1)) * p * q;
			g = (float)(bit);
			*img0++ = g;
		}
	}

	return 1;
}

/***
v2 skip 0 values
***/
inline void correlation_coefficient_FLT_v2(float* img1, int col1, int row1, float* img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height, float* ptr_cr, float h0, float h1)
{
	int i, j, n;
	int	height2, width2;
	float tpi;
	double Sx, Sxx, Sy, Syy, Sxy, v1, v2, cv;
	float* image1, * image10, * image2, * image20;
	float fValue2;
	bool* bMask = (bool*)calloc(height * width, sizeof(bool));

	height2 = height / 2;
	width2 = width / 2;

	n = width * height;

	// v2: get valid data mask
	for (i = 0; i < n; i++)
		bMask[i] = true;

	image10 = image1 = img1 + (y1 - height2) * col1 + (x1 - width2);
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (ABS((*image1)) < 1e-6)
				bMask[i * width + j] = false;
			image1++;
		}
		image10 += col1;
		image1 = image10;
	}

	image20 = image2 = img2 + (y2 - height2) * col2 + (x2 - width2);
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (ABS((*image2)) < 1e-6)
				bMask[i * width + j] = false;
			image2++;
		}
		image20 += col2;
		image2 = image20;
	}

	Sx = Sxx = 0.0;
	image10 = image1 = img1 + (y1 - height2) * col1 + (x1 - width2);
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (bMask[i * width + j] == true)
			{
				Sx += *image1;
				Sxx += (*image1) * (*image1);
			}
			image1++;
		}
		image10 += col1;
		image1 = image10;
	}
	v1 = (Sxx - (Sx * Sx) / n) / n;

	Sy = Syy = 0.0;
	image20 = image2 = img2 + (y2 - height2) * col2 + (x2 - width2);
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (bMask[i * width + j] == true)
			{
				fValue2 = (*image2) * h1 + h0;
				Sy += fValue2;
				Syy += fValue2 * fValue2;
			}
			image2++;
		}
		image20 += col2;
		image2 = image20;
	}
	v2 = (Syy - (Sy * Sy) / n) / n;

	Sxy = 0;
	image10 = image1 = img1 + (y1 - height2) * col1 + (x1 - width2);
	image20 = image2 = img2 + (y2 - height2) * col2 + (x2 - width2);

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (bMask[i * width + j] == true)
			{
				fValue2 = (*image2) * h1 + h0;
				Sxy += (*image1) * fValue2;
			}
			image1++;
			image2++;
		}
		image10 += col1;  image1 = image10;
		image20 += col2;  image2 = image20;
	}
	cv = (Sxy - (Sx * Sy) / n) / n;
	if (v1 < 1 || v2 < 1 || cv <= 0.0)  tpi = 0;
	else
		tpi = (float)(cv / sqrt(v1 * v2));

	*ptr_cr = tpi;

	free(bMask);

	return;
}

/***
do not used std; add fShift
***/
inline bool normalized_FLT_simple(float* data, int n, float fFillValue, float fShift, bool bMedian)
{
	float* data_valid = NULL;
	int i;
	int n_valid, valid_index;
	float mean;
	float pfStats[4];

	data_valid = (float*)calloc(n, sizeof(float));

	n_valid = 0;
	for (i = 0; i < n; i++)
	{
		if (ABS(data[i] - fFillValue) > 1e-6)
		{
			data_valid[n_valid] = data[i];
			n_valid += 1;
		}
	}

	if (n_valid <= 2)
	{
		free(data_valid);
		return false;
	}

	// get mean and std
	if (bMedian == false)
		mean = Mean1_FLT(data_valid, n_valid);
	else
	{
		GetStatistics_FLT(data_valid, n_valid, fFillValue, pfStats);
		mean = pfStats[3]; // median
	}

	// normalize
	for (i = 0; i < n_valid; i++)
		data_valid[i] = (data_valid[i] - mean) + fShift;

	// save data back
	valid_index = 0;
	for (i = 0; i < n; i++)
	{
		if (ABS(data[i] - fFillValue) > 1e-6)
		{
			data[i] = data_valid[valid_index];
			valid_index += 1;
		}
	}

	free(data_valid);
	return true;
}

/***
Get an image subset
***/
inline int imsub_UC(unsigned char* im, int ncol, int nrow, int x, int y, int h, unsigned char* sub)
{
	int i, j;
	int w;

	w = h * 2 + 1;
	memset(sub, 0, w * w * sizeof(unsigned char));

	if (x<h || x>ncol - h - 1 || y<h || y>nrow - h - 1)
		return 0;

	for (i = 0; i < w; i++)
	{
		for (j = 0; j < w; j++)
			sub[i * w + j] = im[(y - h + i) * ncol + (x - h + j)];
	}

	return 1;
}
#endif