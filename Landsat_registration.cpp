#include "Landsat_registration.h"
#include "LSDepthFirst.h"

void SetNewLandsatImage(Landsat_image* pLandsatImage)
{
	pLandsatImage->path = -1;
	pLandsatImage->row = -1;
	pLandsatImage->DOY = -1;
	pLandsatImage->resolution = -1;
	pLandsatImage->width = -1;
	pLandsatImage->height = -1;
	pLandsatImage->year = -1;
	pLandsatImage->proc_level = enum_unknown;
	pLandsatImage->ID = -1;

	return;
}

void InitializeLandsatImage(Landsat_image* pLandsatImage, char* cLandsatImagePathName)
{
	char cData[STRLEN];
	char cFileName[STRLEN];

	strncpy_s(pLandsatImage->cImageFilePathName, STRLEN, cLandsatImagePathName, strlen(cLandsatImagePathName));
	GetFileDir(pLandsatImage->cImageFilePathName, pLandsatImage->cImageDir);
	GetFileName(pLandsatImage->cImageFilePathName, pLandsatImage->cImageFileName);

	strncpy_s(cFileName, STRLEN, pLandsatImage->cImageFileName, strlen(pLandsatImage->cImageFileName));

	// path
	strncpy_s(cData, STRLEN, cFileName + 10, 3);
	pLandsatImage->path = atoi(cData);
	// row
	strncpy_s(cData, STRLEN, cFileName + 13, 3);
	pLandsatImage->row = atoi(cData);
	// year
	strncpy_s(cData, STRLEN, cFileName + 17, 4);
	pLandsatImage->year = atoi(cData);
	// month
	strncpy_s(cData, STRLEN, cFileName + 21, 2);
	pLandsatImage->month = atoi(cData);
	// day
	strncpy_s(cData, STRLEN, cFileName + 23, 2);
	pLandsatImage->day = atoi(cData);
	// DOY
	pLandsatImage->DOY = ConvertToDOY(pLandsatImage->year, pLandsatImage->month, pLandsatImage->day);
	// sensor
	pLandsatImage->sensor = GetLandsatSensorFromFileName(cFileName);
	// processing level
	pLandsatImage->proc_level = GetLandsatProcLevelFromFileName(cFileName);
	// resolution
	strncpy_s(cData, STRLEN, cFileName + strlen(cFileName)-1, 1);
	pLandsatImage->resolution = GetLandsatResolution(pLandsatImage->sensor, atoi(cData));

	pLandsatImage->ID = pLandsatImage->year * 1000 + pLandsatImage->DOY;
	pLandsatImage->ID_date = pLandsatImage->year * 10000 + pLandsatImage->month*100 + pLandsatImage->day;

	return;
}

bool IsLandsatImageValid(Landsat_image* pLandsatImage)
{
	if (pLandsatImage->path > 0
		&& pLandsatImage->row > 0
		&& pLandsatImage->DOY > 0
		&& pLandsatImage->resolution > 0
		&& pLandsatImage->width > 0
		&& pLandsatImage->height > 0
		&& pLandsatImage->year > 0)
		return true;
	else
		return false;
}

enum_Landsat_sensor GetLandsatSensorFromFileName(char* cFileName)
{
	enum_Landsat_sensor sensor = enum_unknown_sensor;
	if (strstr(cFileName, "LM01") != NULL)
		sensor = enum_L1MSS;
	else if (strstr(cFileName, "LM02") != NULL)
		sensor = enum_L2MSS;
	else if (strstr(cFileName, "LM03") != NULL)
		sensor = enum_L3MSS;
	else if (strstr(cFileName, "LM04") != NULL)
		sensor = enum_L4MSS;
	else if (strstr(cFileName, "LM05") != NULL)
		sensor = enum_L5MSS;
	else if (strstr(cFileName, "LC08") != NULL)
		sensor = enum_L8;
	else if (strstr(cFileName, "LC09") != NULL)
		sensor = enum_L9;

	return sensor;
}

enum_Landsat_proc_level GetLandsatProcLevelFromFileName(char* cFileName)
{
	enum_Landsat_proc_level level = enum_unknown;
	if (strstr(cFileName, "L1GS") != NULL)
		level = enum_L1GS;
	else if (strstr(cFileName, "L1TP") != NULL)
		level = enum_L1TP;

	return level;
}

int GetLandsatResolution(enum_Landsat_sensor sensor, int band)
{
	int resolution;

	resolution = -1;
	if (sensor == enum_L4MSS)
	{
		if (band == 4)
			resolution = 60;
	}

	return resolution;
}

bool IsGoodMatching_MSS(bool bSamePathRow, double RMSE, double dMeanShift_x, double dMeanShift_y)
{
	double meanshift;
	meanshift = sqrt(dMeanShift_x * dMeanShift_x + dMeanShift_y * dMeanShift_y);

	if (bSamePathRow && RMSE < 0.7 && meanshift < 20)
		return true;
	else
		return false;
}

int GetShifts(char* pacDenseMatchingFileName, int iColMax, int iRowMax,	int* pdCols, int* pdRows, int iSamplingInterval, float* ptr_fShifts_x, float* ptr_fShifts_y)
{
	float* pfParalaxMap_x = NULL, * pfParalaxMap_y = NULL, * pfCorrMap = NULL;
	int n;
	int ntie;
	int Row, Col;
	//	int iSamplingInterval = DENSE_MATCHING_SAMPLING_INTERVAL_MSS;

	// input densing matching tie points coordinates and coefficients
	pfParalaxMap_x = (float*)ReadInputImage(pacDenseMatchingFileName, iRowMax, iColMax * 4, 3);
	pfParalaxMap_y = pfParalaxMap_x + iRowMax * iColMax;
	pfCorrMap = pfParalaxMap_x + 3 * iRowMax * iColMax;

	// get ntie
	ntie = 0;
	for (n = 0; n < iRowMax * iColMax; n++)
	{
		if (pfCorrMap[n] > 0.1f)
			ntie += 1;
	}

	// set matched coordinates
	ntie = 0;
	for (Row = 0; Row < iRowMax; Row++)
	{
		for (Col = 0; Col < iColMax; Col++)
		{
			n = Row * iColMax + Col;
			if (pfCorrMap[n] > 0.1f)
			{
				pdRows[ntie] = Row * iSamplingInterval;
				pdCols[ntie] = Col * iSamplingInterval;
				ptr_fShifts_x[ntie] = pfParalaxMap_x[n];
				ptr_fShifts_y[ntie] = pfParalaxMap_y[n];
				ntie += 1;
			}
		}
	}

	free(pfParalaxMap_x);

	return ntie;
}

/***
v1 (9/4/2019):
-add parameter fSAMThreshold;
-get threhsold that is the smaller of fSAMThreshold and SAM stat (max of median and mean), and discard matchings with SAMs smaller than the threshold
***/
int GetShifts_v1(char* pacDenseMatchingFileName, int iColMax, int iRowMax, int* pdCols, int* pdRows, int iSamplingInterval, float fSAMThreshold, float* ptr_fShifts_x, float* ptr_fShifts_y)
{
	float* pfParalaxMap_x = NULL, * pfParalaxMap_y = NULL, * pfCorrMap = NULL;
	int n;
	int ntie;
	int Row, Col;
	float pfSAM_stats[4];
	float median_SAM;
	float SAM_threshold;

	// input densing matching tie points coordinates and coefficients
	pfParalaxMap_x = (float*)ReadInputImage(pacDenseMatchingFileName, iRowMax, iColMax * 4, 3);
	pfParalaxMap_y = pfParalaxMap_x + iRowMax * iColMax;
	pfCorrMap = pfParalaxMap_x + 3 * iRowMax * iColMax;

	if (GetStatistics_FLT(pfCorrMap, iRowMax*iColMax, 0, pfSAM_stats) == false)
	{
		printf("Error in GetShifts_v1(): unable to get SAM statistics.\n");
		scanf_s(" %d", &n);
		exit(1);
	}

	// get SAM threshold
	median_SAM = MAX(pfSAM_stats[2], pfSAM_stats[3]);
	SAM_threshold = MIN(median_SAM, fSAMThreshold);

	// set matched coordinates
	ntie = 0;
	for (Row = 0; Row < iRowMax; Row++)
	{
		for (Col = 0; Col < iColMax; Col++)
		{
			n = Row * iColMax + Col;
			if (pfCorrMap[n] > SAM_threshold)
			{
				pdRows[ntie] = Row * iSamplingInterval;
				pdCols[ntie] = Col * iSamplingInterval;
				ptr_fShifts_x[ntie] = pfParalaxMap_x[n];
				ptr_fShifts_y[ntie] = pfParalaxMap_y[n];
				ntie += 1;
			}
		}
	}

	free(pfParalaxMap_x);

	return ntie;
}

///////////////////////////////////////////////////////////////////////////
Landsat_registration::Landsat_registration()
{
	int i;
	
	m_ImagesNum_L1TP = 0;
	m_width = -1;
	m_height = -1;

	memset(m_pImages_L1TP, 0, sizeof(Landsat_image) * MAX_PROC_LANDSAT_IMAGES);
	for (i = 0; i < MAX_PROC_LANDSAT_IMAGES; i++)
	{
		SetNewLandsatImage(m_pImages_L1TP + i);
	}

	memset(m_misregistrations_L1TP, 0, sizeof(double) * MAX_PROC_LANDSAT_IMAGES);

	return;
}

Landsat_registration::Landsat_registration(char* cFileList_L1TP, char * cFileList_L1GS, char *cOutputDir, int width, int height, int iDenseMatchingInterval)
{
	int i;

	m_width = width;
	m_height = height;
	m_iDenseMatchingInterval = iDenseMatchingInterval;

	sprintf_s(m_cOutputDir, STRLEN, "%s", cOutputDir);

	for (i = 0; i < MAX_PROC_LANDSAT_IMAGES; i++)
	{
		SetNewLandsatImage(m_pImages_L1TP + i);
		SetNewLandsatImage(m_pImages_L1GS + i);
	}

	m_ImagesNum_L1TP = 0;
	m_ImagesNum_L1GS = 0;
	if (strlen(cFileList_L1TP) > 1)
		InputImagesList(cFileList_L1TP, enum_L1TP);
	if (strlen(cFileList_L1GS) > 1)
		InputImagesList(cFileList_L1GS, enum_L1GS);

	// put m_pImages_L1TP and m_ImagesNum_L1GS in m_pImages (pointers only)
	for (i = 0; i < MAX_PROC_LANDSAT_IMAGES; i++)
		m_pImages[i] = NULL;
	for (i = 0; i < m_ImagesNum_L1TP; i++)
		m_pImages[i] = m_pImages_L1TP + i;
	for (i = 0; i < m_ImagesNum_L1GS; i++)
		m_pImages[i+ m_ImagesNum_L1TP] = m_pImages_L1GS + i;

	m_ImagesNum = m_ImagesNum_L1TP + m_ImagesNum_L1GS;

	memset(m_misregistrations_L1TP, 0, sizeof(double) * MAX_PROC_LANDSAT_IMAGES);

	return;
}

Landsat_registration::~Landsat_registration()
{
	return;
}

// find index in m_pImages
int Landsat_registration::GetImageIndex(int imageID)
{
	int i, index;
	index = -1;
	for (i=0; i< m_ImagesNum; i++)
	{
		if (m_pImages[i] == NULL)
			continue;

		if (m_pImages[i]->ID == imageID)
		{
			index = i;
			return index;
		}
	}

	return index;
}

// find index in m_pImages_L1TP
int Landsat_registration::GetL1TPImageIndex(int imageID)
{
	int i, index;
	index = -1;
	for (i = 0; i < m_ImagesNum_L1TP; i++)
	{
		if (m_pImages_L1TP[i].ID < 0)
			continue;

		if (m_pImages_L1TP[i].ID == imageID)
		{
			index = i;
			return index;
		}
	}

	return index;
}

Landsat_image* Landsat_registration::GetImage(int index)
{
	if (index < 0 || index >= m_ImagesNum)
		return NULL;

	return m_pImages[index];
}

int Landsat_registration::InputImagesList(char* cFileList, enum_Landsat_proc_level proc_level)
{
	FILE* fin = NULL;
	char cLine[STRLEN];
	int iData;
	char* pacPos = NULL;
	Landsat_image* pImages = NULL;
	int* pImage_num = NULL;

	if (!FileExist(cFileList))
	{
		printf("%s not exist.\n", cFileList);
		scanf_s(" %d", &iData);
		exit(0);
		return -1;
	}
	
	if (Gettxtlines(cFileList) == 0)
		return 0;

	if (proc_level == enum_L1TP)
	{
		pImage_num = &m_ImagesNum_L1TP;
		pImages = m_pImages_L1TP;
	}
	else
	{
		pImage_num = &m_ImagesNum_L1GS;
		pImages = m_pImages_L1GS;
	}

	// read registration file
	fin = Readtxt(cFileList);
	*pImage_num = 0;
	while (fgets(cLine, STRLEN, fin) != NULL)
	{
		cLine[strlen(cLine) - 1] = '\0';

		InitializeLandsatImage(pImages + *pImage_num, cLine);
		pImages[*pImage_num].width = this->m_width;
		pImages[*pImage_num].height = this->m_height;

		*pImage_num += 1;
	}
	fclose(fin);

	return *pImage_num;
}

/***
* v1:
-skip if i1 is later than i2
-no longer have cases of reversed image 1 and image 2 order
-7/28/2021: call DepthFirstImageRegistration_UC() instead of DepthFirstImageRegistration_v1()
8/17/2022
- change bSameSensor to iSensorCombine (1: free combine; 3: absolute shift matched to L8/9
***/
int Landsat_registration::L1TPScanning_v1(bool bOutputMatchingSummaryFile, int iSensorCombine, int iHalfSearchWindow)
{
	int matched_L1TP_num;
	int i1, i2;
	Landsat_image* pImage1 = NULL, * pImage2 = NULL;
	char cRegFileName[STRLEN], cRegFilePathName[STRLEN], cRegDir[STRLEN], cDir_output[STRLEN];

	LSDepthFirstRegistration* pLSReg = new LSDepthFirstRegistration;
	double dMeanShift_x, dMeanShift_y, RMSE, RMSE_affine = 0.0;
	int ntie;
	FILE* fout = NULL;
	char pacMatchingSummaryFileName[STRLEN];

	double dMeanShifts_x[MAX_PROC_LANDSAT_IMAGES], dMeanShifts_y[MAX_PROC_LANDSAT_IMAGES];
	int nCount[MAX_PROC_LANDSAT_IMAGES];
	int iHalfSearchWindow_;

	// registration/machings
	sprintf_s(cRegDir, STRLEN, "%s/matchings", this->m_cOutputDir);
	_mkdir(cRegDir);

	if (iSensorCombine == 1)
		sprintf_s(pacMatchingSummaryFileName, STRLEN, "%s/L1TP_matching_rsp%d.txt", this->m_cOutputDir, this->m_iDenseMatchingInterval);
	else if (iSensorCombine == 3)
		sprintf_s(pacMatchingSummaryFileName, STRLEN, "%s/L1TP_matching_rsp%d_abs.txt", this->m_cOutputDir, this->m_iDenseMatchingInterval); // rsp100 is L8/9 matching at 30m resolution

	if (bOutputMatchingSummaryFile == true)
		fout = Writetxt(pacMatchingSummaryFileName);

	memset(dMeanShifts_x, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
	memset(dMeanShifts_y, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
	memset(nCount, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(int));

	for (i1 = 0; i1 < this->m_ImagesNum_L1TP; i1++)
	{
		pImage1 = this->m_pImages_L1TP + i1;
		if (pImage1->proc_level != enum_L1TP)
			continue;

		for (i2 = 0; i2 < this->m_ImagesNum_L1TP; i2++)
		{
			if (i1 == i2)
				continue;

			pImage2 = this->m_pImages_L1TP + i2;
			if (pImage2->proc_level != enum_L1TP)
				continue;

			if (iSensorCombine == 3 && (pImage1->sensor != 8 && pImage1->sensor != 9 && pImage2->sensor != 8 && pImage2->sensor != 9))
				continue;

			if (pImage1->year * 1000 + pImage1->DOY >= pImage2->year * 1000 + pImage2->DOY)
				continue;

			if (pImage1->sensor >= enum_L8 && pImage2->sensor >= enum_L8)
				iHalfSearchWindow_ = 1;
			else
				iHalfSearchWindow_ = iHalfSearchWindow;
			
			sprintf_s(cDir_output, STRLEN, "%s/%d%02d%02d-L%d-%d_%d%02d%02d-L%d-%d", cRegDir, pImage1->year, pImage1->month, pImage1->day, pImage1->sensor, pImage1->proc_level, pImage2->year, pImage2->month, pImage2->day, pImage2->sensor, pImage2->proc_level);
			_mkdir(cDir_output);
			sprintf_s(cRegFileName, STRLEN, "LSReg.txt");
			CreateRegitrationProjectFile(cDir_output, pImage1->cImageFilePathName, pImage2->cImageFilePathName, m_width, m_height, this->m_iDenseMatchingInterval, 3, 2, 0.985f, cRegFileName);
			sprintf_s(cRegFilePathName, STRLEN, "%s/LSReg.txt", cDir_output);
			if (DepthFirstImageRegistration_UC_simple(pLSReg, cRegFilePathName, &dMeanShift_x, &dMeanShift_y, &RMSE, &ntie, iHalfSearchWindow_))
			{
				if (bOutputMatchingSummaryFile == true)
					fprintf(fout, "%d %d\t%.3f\t%.3f\t%.3f %d %s\n", i1 + 1, i2 + 1, dMeanShift_x, dMeanShift_y, RMSE, ntie, cDir_output);

				dMeanShifts_x[i1] += dMeanShift_x;
				dMeanShifts_y[i1] += dMeanShift_y;
				nCount[i1] += 1;

				dMeanShifts_x[i2] += (-dMeanShift_x);
				dMeanShifts_y[i2] += (-dMeanShift_y);
				nCount[i2] += 1;
			}
		}
	}

	//	fprintf(fout, "\n");
	matched_L1TP_num = 0;
	for (i1 = 0; i1 < this->m_ImagesNum_L1TP; i1++)
	{
		if (nCount[i1] > 0)
		{
			dMeanShifts_x[i1] /= nCount[i1];
			dMeanShifts_y[i1] /= nCount[i1];
		}

		if (nCount[i1] > 0)
			matched_L1TP_num += 1;
	}
	if (bOutputMatchingSummaryFile == true)
		fclose(fout);

	delete pLSReg;

	return matched_L1TP_num;
}

/***
Check L1TP scanning results, get m_misregistrations_L1TP and output to L1TP_misregistration.txt 
1/3/2020: use strict dMeanShifts_abs as in paper to replace dMeanShifts
***/
int Landsat_registration::GetL1TPMisregistrationsFromScanningResults()
{
	FILE* file_matching = NULL;
	char pacMatchingFileDir[STRLEN];
	char pacMatchingSummaryFileName[STRLEN];
	int iMatchedTilesNum, iMatch;
	int n_misalign;
	double dMeanShift_x, dMeanShift_y, RMSE;
	int ntie;
	int i1, i2;
	double dMeanShifts_x[MAX_PROC_LANDSAT_IMAGES], dMeanShifts_y[MAX_PROC_LANDSAT_IMAGES], dMeanShifts[MAX_PROC_LANDSAT_IMAGES], dMeanShifts_abs[MAX_PROC_LANDSAT_IMAGES];
	int nCount[MAX_PROC_LANDSAT_IMAGES];
	double dMaxShift;
	int maxShiftIndex;

	bool flags_misalign[MAX_PROC_LANDSAT_IMAGES];
	memset(flags_misalign, 0, sizeof(bool) * MAX_PROC_LANDSAT_IMAGES);

	sprintf_s(pacMatchingSummaryFileName, STRLEN, "%s/L1TP_matching_rsp%d.txt", this->m_cOutputDir, this->m_iDenseMatchingInterval);
	iMatchedTilesNum = Gettxtlines(pacMatchingSummaryFileName);	// number of matchings

	memset(dMeanShifts, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
	do
	{
		memset(dMeanShifts_x, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
		memset(dMeanShifts_y, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
		memset(dMeanShifts_abs, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
		memset(nCount, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(int));

		file_matching = Readtxt(pacMatchingSummaryFileName);

		// input x and y shifts
		for (iMatch = 0; iMatch < iMatchedTilesNum; iMatch++)
		{
			// note: dMeanShift_x and dMeanShift_y are fitted translation tranformation coefficients a0 and b0
			fscanf_s(file_matching, "%d%d%lf%lf%lf%d", &i1, &i2, &dMeanShift_x, &dMeanShift_y, &RMSE, &ntie);
			fscanf_s(file_matching, "%s", pacMatchingFileDir, STRLEN);

			if (RMSE > 0.7)
				continue;

			i1 = i1 - 1;
			i2 = i2 - 1;

			if (flags_misalign[i1] == false && flags_misalign[i2] == false)
			{
				// both good
				dMeanShifts_x[i1] += dMeanShift_x;
				dMeanShifts_y[i1] += dMeanShift_y;
				dMeanShifts_abs[i1] += sqrt(dMeanShift_x * dMeanShift_x + dMeanShift_y * dMeanShift_y);
				nCount[i1] += 1;

				dMeanShifts_x[i2] += (-dMeanShift_x);
				dMeanShifts_y[i2] += (-dMeanShift_y);
				dMeanShifts_abs[i2] += sqrt(dMeanShift_x * dMeanShift_x + dMeanShift_y * dMeanShift_y);
				nCount[i2] += 1;
			}
			else if (flags_misalign[i1] == false && flags_misalign[i2] == true)
			{
				// i1 is good and i2 is misaligned, count i2 w.r.t. i1
				dMeanShifts_x[i2] += (-dMeanShift_x);
				dMeanShifts_y[i2] += (-dMeanShift_y);
				dMeanShifts_abs[i2] += sqrt(dMeanShift_x * dMeanShift_x + dMeanShift_y * dMeanShift_y);
				nCount[i2] += 1;
			}
			else if (flags_misalign[i1] == true && flags_misalign[i2] == false)
			{
				// i2 is good and i1 is misaligned, count i1 w.r.t. i2
				dMeanShifts_x[i1] += dMeanShift_x;
				dMeanShifts_y[i1] += dMeanShift_y;
				dMeanShifts_abs[i1] += sqrt(dMeanShift_x * dMeanShift_x + dMeanShift_y * dMeanShift_y);
				nCount[i1] += 1;
			}
		}
		fclose(file_matching);

		// calculate mean shifts per image
		dMaxShift = 0;
		maxShiftIndex = -1;
		for (i1 = 0; i1 < this->m_ImagesNum_L1TP; i1++)
		{
			if (nCount[i1] > 0)
			{
				dMeanShifts_x[i1] /= nCount[i1];
				dMeanShifts_y[i1] /= nCount[i1];
				dMeanShifts_abs[i1] /= nCount[i1];
				dMeanShifts[i1] = sqrt(dMeanShifts_x[i1] * dMeanShifts_x[i1] + dMeanShifts_y[i1] * dMeanShifts_y[i1]);
				if (flags_misalign[i1] == false && dMeanShifts_abs[i1] > dMaxShift) // 1/3/2020: replace dMeanShifts with dMeanShifts_abs
				{
					dMaxShift = dMeanShifts_abs[i1];
					maxShiftIndex = i1;
				}
			}
		}	
		
		// check whether the max-shift image is misaligned; if so, flag it
		n_misalign = 0;
		if (dMaxShift > LANDSAT_MSS_MISALIGN_THRESHOLD)
		{
			n_misalign += 1;
			flags_misalign[maxShiftIndex] = true; // flag misaligned L1TP
		}
	} while (n_misalign > 0);

	// set misregistrations
	n_misalign = 0;
	for (i1 = 0; i1 < this->m_ImagesNum_L1TP; i1++)
	{
	//	m_misregistrations_L1TP[i1] = sqrt(dMeanShifts_x[i1] * dMeanShifts_x[i1] + dMeanShifts_y[i1] * dMeanShifts_y[i1]);
		m_misregistrations_L1TP[i1] = dMeanShifts_abs[i1];
		n_misalign += (int)flags_misalign[i1];
	}

	return n_misalign;
}

/***
Revised from Sentinel_UTM_Registration::GetConnectionGraphFromMatchingFile()
Note: pbConnectionMap is sized iSAFE_Image_Num x iSAFE_Image_Num
***/
int Landsat_registration::GetConnectionGraphFromMatchingFile(char* pacMatchingFileName, bool* ptr_bConnectionMap, bool bIndexFrom_0)
{
	int iMatchedTilesNum;
	FILE* file_matching = NULL;
	int iMatch;
	int index1, index2;
	int ntie;
	double RMSE, dMeanShift_x, dMeanShift_y;
	char pacMatchingFileDir[STRLEN];
	Landsat_image* pImage1 = NULL, * pImage2 = NULL;
	int iMatch_valid;
	int n_Images;
	int i;

	n_Images = this->m_ImagesNum;

	for (i = 0; i < n_Images * n_Images; i++)
		ptr_bConnectionMap[i] = false;

	iMatchedTilesNum = Gettxtlines(pacMatchingFileName);	// number of matched tiles
	file_matching = Readtxt(pacMatchingFileName);
	iMatch_valid = 0;
	for (iMatch = 0; iMatch < iMatchedTilesNum; iMatch++)
	{
		// note: dMeanShift_x and dMeanShift_y are fitted translation tranformation coefficients a0 and b0	
		fscanf_s(file_matching, "%d%d%lf%lf%lf%d", &index1, &index2, &dMeanShift_x, &dMeanShift_y, &RMSE, &ntie);
		fscanf_s(file_matching, "%s", pacMatchingFileDir, STRLEN);
		if (bIndexFrom_0 == false)
		{
			index1 -= 1;
			index2 -= 1;
		}

		pImage1 = this->GetImage(index1);
		pImage2 = this->GetImage(index2);
		
		if (pImage1 == NULL || pImage2 == NULL)
		{
			printf("Error in GetConnectionGraphFromMatchingFile(): image %d or %d not found\n", index1, index2);
			scanf_s(" %d", &i);
			exit(1);
		}

		if (IsGoodMatching_MSS(true, RMSE, dMeanShift_x, dMeanShift_y) == false)
			continue;

		ptr_bConnectionMap[index1 * n_Images + index2] = true;
		ptr_bConnectionMap[index2 * n_Images + index1] = true;
		iMatch_valid += 1;
	}

	fclose(file_matching);

	return iMatch_valid;
}

/***
Created 7/29/2021
Revised from GetConnectionGraphFromMatchingFile()
For each image, get tie points number in each RBF grid, considering all matched images
piRBFs_ntie_Map is sized: this->m_ImagesNum x iRBFs_K
iRBFs_K = n x n
v1 (8/2/2021): return xs_k_valid_nImages, ys_k_valid_nImages, K_valid_nImages
***/
int Landsat_registration::FindValid_RBF_Centers_PerImage_v1(char* pacMatchingSummaryFileName, int RBFs_K, unsigned int* piRBFs_ntie_Map, unsigned int iRBF_ties_num_threshold, double* xs_k_valid_nImages, double* ys_k_valid_nImages, int* K_valid_nImages)
{
	int iMatchedTilesNum;
	FILE* file_matching = NULL;
	int iMatch;
	int index1, index2;
	int ntie;
	double RMSE, dMeanShift_x, dMeanShift_y;
	char pacMatchingFileDir[STRLEN];
	Landsat_image* pImage1 = NULL, * pImage2 = NULL;
	int iMatch_valid;
	int n_Images;
	int i, m, k;
	float* pfShifts_x = NULL, * pfShifts_y = NULL;
	int* piMatchingCols = NULL, * piMatchingRows = NULL;
	int iMatchingColMax, iMatchingRowMax;
	int iColMax, iRowMax;
	double x_k[MAX_RBFs_K], y_k[MAX_RBFs_K];
	int grid_width, grid_height;
	int iDenseMatchingInterval;
	char pacMatchingFileName[STRLEN], pacOutputFileName[STRLEN];
	unsigned int* piRBFs_ntie_Map_input = NULL;
	FILE* fout = NULL;

	int iPtNum;
	double xs_k[MAX_RBFs_K], ys_k[MAX_RBFs_K]; // 4*K coefficients	int grid_width, grid_height;
	int n_K_total, n_K_cur;

	if (RBFs_K <= 0)
		return 0;

	iColMax = this->m_width;
	iRowMax = this->m_height;

	// get RBF centers and grid sizes
	iPtNum = GetGridded_RBF_Centers(iColMax, iRowMax, RBFs_K, xs_k, ys_k, &grid_width, &grid_height);

	iDenseMatchingInterval = this->m_iDenseMatchingInterval;
	n_Images = this->m_ImagesNum;

	iMatchingColMax = iColMax / iDenseMatchingInterval;
	iMatchingRowMax = iRowMax / iDenseMatchingInterval;

	sprintf_s(pacOutputFileName, STRLEN, "%s/L1TP_K%d_rsp%d_grid_ties_stack_%dimages", this->m_cOutputDir, RBFs_K, iDenseMatchingInterval, n_Images);
	if (FileExist(pacOutputFileName))
	{
		printf("File exists.\n");
		piRBFs_ntie_Map_input = (unsigned int*)ReadInputImage(pacOutputFileName, RBFs_K, n_Images, 2);
		memmove(piRBFs_ntie_Map, piRBFs_ntie_Map_input, RBFs_K * n_Images * sizeof(unsigned int));

		// get the numers of the valid K centers per image, and corresponding xs and ys
		n_K_total = 0;
		for (i = 0; i < n_Images; i++)
		{
			n_K_cur = 0;
			for (k = 0; k < RBFs_K; k++)
			{
				if (piRBFs_ntie_Map[i * RBFs_K + k] >= iRBF_ties_num_threshold)
				{
					xs_k_valid_nImages[i * RBFs_K + n_K_cur] = xs_k[k];
					ys_k_valid_nImages[i * RBFs_K + n_K_cur] = ys_k[k];
					n_K_cur += 1;
				}
			}
			K_valid_nImages[i] = n_K_cur;
			n_K_total += n_K_cur;
		}
		free(piRBFs_ntie_Map_input);
		return n_K_total;
	}

	// allocate memory
	if (!(pfShifts_x = (float*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(float)))
		|| !(pfShifts_y = (float*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(float)))
		|| !(piMatchingCols = (int*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(int)))
		|| !(piMatchingRows = (int*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(int))))
	{
		printf("\nError in FindValid_RBF_Centers_PerImage(): insufficient memory.\n");
		scanf_s(" %d", &i);
		exit(1);
	}

	// get RBF centers and grid sizes
	GetGridded_RBF_Centers(iColMax, iRowMax, RBFs_K, x_k, y_k, &grid_width, &grid_height);

	iMatchedTilesNum = Gettxtlines(pacMatchingSummaryFileName);	// number of matched tiles
	file_matching = Readtxt(pacMatchingSummaryFileName);
	iMatch_valid = 0;
	memset(piRBFs_ntie_Map, 0, n_Images * RBFs_K * sizeof(unsigned int)); // e.g. piRBFs_ntie_Map is sized: 10 (images) x 64
	for (iMatch = 0; iMatch < iMatchedTilesNum; iMatch++)
	{
		// note: dMeanShift_x and dMeanShift_y are fitted translation tranformation coefficients a0 and b0	
		fscanf_s(file_matching, "%d%d%lf%lf%lf%d", &index1, &index2, &dMeanShift_x, &dMeanShift_y, &RMSE, &ntie);
		fscanf_s(file_matching, "%s", pacMatchingFileDir, STRLEN);
		index1 = index1 - 1; // L1TP_matching.txt is indexed starting from 1; L1TP_L1GS_matching is indexed from 0
		index2 = index2 - 1;
		if (index1 < 0 || index2 < 0 || index1 > n_Images - 1 || index2 > n_Images - 1)
		{
			printf("\nError in FindValid_RBF_Centers_PerImage(): wrong image index.\n");
			scanf_s(" %d", &i);
			exit(1);
		}
		sprintf_s(pacMatchingFileName, STRLEN, "%s/matching/layer0_dense_matching_rsp%d", pacMatchingFileDir, iDenseMatchingInterval);
		ntie = GetShifts_v1(pacMatchingFileName, iMatchingColMax, iMatchingRowMax, piMatchingCols, piMatchingRows, iDenseMatchingInterval, 0.1f, pfShifts_x, pfShifts_y);

		// assign dense matching points to grids of images index1 and index2 (does not matter if either is a control)
		for (i = 0; i < ntie; i++)
		{
			for (m = 0; m < RBFs_K; m++)
			{
				// 9/18/2021, changed < to <=
				if (ABS(piMatchingCols[i] - x_k[m]) <= grid_width / 2 && ABS(piMatchingRows[i] - y_k[m]) <= grid_height / 2)
					piRBFs_ntie_Map[RBFs_K * index1 + m] += 1;

				if (ABS(piMatchingCols[i] + pfShifts_x[i] - x_k[m]) <= grid_width / 2 && ABS(piMatchingRows[i] + pfShifts_y[i] - y_k[m]) <= grid_height / 2)
					piRBFs_ntie_Map[RBFs_K * index2 + m] += 1;
			}
		}
	}
	fclose(file_matching);

	fout = WriteBinary(pacOutputFileName);
	fwrite(piRBFs_ntie_Map, sizeof(unsigned int), RBFs_K * n_Images, fout);
	fclose(fout);

	// get the numers of the valid K centers per image, and corresponding xs and ys
	n_K_total = 0;
	for (i = 0; i < n_Images; i++)
	{
		n_K_cur = 0;
		for (k = 0; k < RBFs_K; k++)
		{
			if (piRBFs_ntie_Map[i * RBFs_K + k] >= iRBF_ties_num_threshold)
			{
				xs_k_valid_nImages[i * RBFs_K + n_K_cur] = xs_k[k];
				ys_k_valid_nImages[i * RBFs_K + n_K_cur] = ys_k[k];
				n_K_cur += 1;
			}
		}
		K_valid_nImages[i] = n_K_cur;
		n_K_total += n_K_cur;
	}

	free(pfShifts_x);
	free(pfShifts_y);
	free(piMatchingCols);
	free(piMatchingRows);

	return n_K_total;
}

/***
Revised from Sentinel_UTM_Registration::CalculateAdjustmentResiduals_Affine_v2_2_5()
Fit to polynomial
Use bForwardModel = true such that the output coefficients do not need to be reversed
9/4/2019: use GetShifts_v1() to remove low-SAM matching points
v1 (1/3/2020): add return parameters ptr_dPre_MisregistrationPerImage and ptr_dPost_MisregistrationPerImage
v2 (8/9/2021)
- add parameters int iTransformationType, piK_valid_nImages, RBFs_K, bIndexedFrom_0
- allow RBFs transformation
- if iTransformationType == 4, piK_valid_nImages is not NULL and RBFs_K > 0, and each image has 12 + 4*RBFs_K coefficients in pCoefs_all
8/10/2021 bug fix: in y residual calculation, changed sqrt(ptr_dResiduals[2 * iMatch] to sqrt(ptr_dResiduals[2 * iMatch + 1]
8/11/2021: if neither two matched images are controls, still consider in the calculation of the per-image mean pre- and post-registration misregistration, i.e. ptr_dPre_MisregistrationPerImage and ptr_dPost_MisregistrationPerImage
***/
bool Landsat_registration::CalculateAdjustmentResiduals_v2(double* pCoefs_all, int* piImages_Control_Flag, int iTransformationType, int * piK_valid_nImages, int RBFs_K, double* ptr_dResiduals, double& ptr_dRMSE, double* ptr_dRMSEsPerImage, double* ptr_dMisregistrationRes_stats, double* ptr_dPre_MisregistrationPerImage, double* ptr_dPost_MisregistrationPerImage, bool bIndexedFrom_0)
{
	FILE* file_matching = NULL;
	char pacMatchingFileName[STRLEN], pacMatchingFileDir[STRLEN];
	char pacMatchingSummaryFileName[STRLEN];
	int iMatchedPairsNum, iMatch;

	int index_image1, index_image2;
	int iRowMax, iColMax;
	int iMatchingRowMax, iMatchingColMax;
	FILE* fout = NULL;

	double* coords = NULL, * coords_unadjusted = NULL;		// unknowns
	int ntie;
	double dRes, dRes_SQR_x, dRes_SQR_y, dSum;

	float* pfShifts_x = NULL, * pfShifts_y = NULL;
	int* piMatchingCols = NULL, * piMatchingRows = NULL;

	Landsat_image* pImage1 = NULL, * pImage2 = NULL, * pImage = NULL;

	double* dCoefs_i1 = NULL, * dCoefs_i2 = NULL;
	int i, n_Images;
	int ntie_total;
	int n_unknowns, n_unknowns_total;
	double RMSE, dMeanShift_x, dMeanShift_y;
	int ntie_;
	int nMatchedPairs;

	double pd_Sum_per_Image[MAX_PROC_LANDSAT_IMAGES];
	int pi_ntie_per_Image[MAX_PROC_LANDSAT_IMAGES];
	double dMisregistrationRes_mean, dMisregistrationRes_means_mean, dMisregistrationRes_means_std;
	double* pdMisregistrationRes_means = NULL;
	double x1, y1, x2, y2;
	double x1_0, y1_0, x2_0, y2_0;

	int nCount[MAX_PROC_LANDSAT_IMAGES];
	int n_coefs_per_image;
	int iDenseMatchingInterval = this->m_iDenseMatchingInterval;

	if (iTransformationType == 1)
		n_coefs_per_image = 2;
	else if (iTransformationType == 2)
		n_coefs_per_image = 6;
	else if (iTransformationType == 3)
		n_coefs_per_image = 12;
	else if (iTransformationType == 4)
		n_coefs_per_image = 12 + 4* RBFs_K;

	n_Images = this->m_ImagesNum;

	iColMax = this->m_width;
	iRowMax = this->m_height;
	iMatchingColMax = iColMax / iDenseMatchingInterval;
	iMatchingRowMax = iRowMax / iDenseMatchingInterval;

	// allocate memory
	if (!(pfShifts_x = (float*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(float)))
		|| !(pfShifts_y = (float*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(float)))
		|| !(piMatchingCols = (int*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(int)))
		|| !(piMatchingRows = (int*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(int))))
	{
		printf("\nError in Sentinel_UTM_Registration::CalculateAdjustmentResiduals_Affine(): insufficient memory.\n");
		scanf_s(" %d", &i);
		exit(1);
	}

	// get to matched images summary file, e.g. UTM_11_matching.txt
	sprintf_s(pacMatchingSummaryFileName, STRLEN, "%s/L1TP_matching_rsp%d.txt", this->m_cOutputDir, this->m_iDenseMatchingInterval);
	iMatchedPairsNum = Gettxtlines(pacMatchingSummaryFileName);	// number of matched tiles

	if (!(pdMisregistrationRes_means = (double*)calloc(iMatchedPairsNum, sizeof(double))))
	{
		printf("\nError in Sentinel_UTM_Registration::CalculateAdjustmentResiduals_Affine(): insufficient memory.\n");
		scanf_s(" %d", &i);
		exit(1);
	}

	file_matching = Readtxt(pacMatchingSummaryFileName);
	dSum = 0;
	ntie_total = 0;
	memset(ptr_dResiduals, 0, 2 * iMatchedPairsNum * sizeof(double));
	memset(pd_Sum_per_Image, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
	memset(pi_ntie_per_Image, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(int));
	nMatchedPairs = 0;

	memset(ptr_dPre_MisregistrationPerImage, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
	memset(ptr_dPost_MisregistrationPerImage, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(double));
	memset(nCount, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(int));
	for (iMatch = 0; iMatch < iMatchedPairsNum; iMatch++)
	{
		fscanf_s(file_matching, "%d%d%lf%lf%lf%d", &index_image1, &index_image2, &dMeanShift_x, &dMeanShift_y, &RMSE, &ntie);
		fscanf_s(file_matching, "%s", pacMatchingFileDir, STRLEN);
		if (bIndexedFrom_0 == false)
		{
			index_image1 -= 1;
			index_image2 -= 1;
		}

		if (IsGoodMatching_MSS(true, RMSE, dMeanShift_x, dMeanShift_y) == false)
			continue; // used same criteria as in Adjustment_v1()

		if (piImages_Control_Flag[index_image1] > 0 && piImages_Control_Flag[index_image2] > 0)
			continue; // added 9/14/2017; both are controls		

		if (piImages_Control_Flag[index_image1] < 0 || piImages_Control_Flag[index_image2] < 0) // 7/30/2020, fixed first "index_image2" to "index_image1"
			continue; // v2_2_5, one SAFE is unconnected

		//if (piImages_Control_Flag[index_image1] == 0 && piImages_Control_Flag[index_image2] == 0)
		//	continue; // added 8/3/2020, does not affect results

		// v1: get pre-registration shifts per image (do not count L1GS <-> L1GS)
		// 8/11/2021: still count if neither images are controls
//		if (piImages_Control_Flag[index_image1] > 0 || piImages_Control_Flag[index_image2] > 0) // at least one image is L1TP reference (control)
		{
			if (piImages_Control_Flag[index_image1] == 0)
			{
				ptr_dPre_MisregistrationPerImage[index_image1] += sqrt(dMeanShift_x * dMeanShift_x + dMeanShift_y * dMeanShift_y);
				nCount[index_image1] += 1;
			}
			if (piImages_Control_Flag[index_image2] == 0)
			{
				ptr_dPre_MisregistrationPerImage[index_image2] += sqrt(dMeanShift_x * dMeanShift_x + dMeanShift_y * dMeanShift_y);
				nCount[index_image2] += 1;
			}
		}

		pImage1 = this->GetImage(index_image1);
		pImage2 = this->GetImage(index_image2);

		// get shift x and y for the tile
		sprintf_s(pacMatchingFileName, STRLEN, "%s/matching/layer0_dense_matching_rsp%d", pacMatchingFileDir, iDenseMatchingInterval);
		ntie = GetShifts_v1(pacMatchingFileName, iMatchingColMax, iMatchingRowMax, piMatchingCols, piMatchingRows, iDenseMatchingInterval, 0.95f, pfShifts_x, pfShifts_y);

		// get coefficients
		dCoefs_i1 = pCoefs_all + index_image1 * n_coefs_per_image;
		dCoefs_i2 = pCoefs_all + index_image2 * n_coefs_per_image;

		// form observation equation for each matched point
		ntie_ = 0;
		dMisregistrationRes_mean = 0;
		for (i = 0; i < ntie; i++)
		{
			//			if (sqrt((dMeanShift_x + pfShifts_x[i]) * (dMeanShift_x + pfShifts_x[i]) + (dMeanShift_y + pfShifts_y[i]) * (dMeanShift_y + pfShifts_y[i]) > 2 * RMSE))
			//				continue;	// added 2/21/2017; disabled 8/30/2019

			ntie_ += 1;

			// image 1
			// get affine-transformed coordinates (in pixels) within the SAFE
			// (tile_h_offset_pixels1 + piMatchingCols[i],  tile_v_offset_pixels1 + piMatchingRows[i]) is the untransformed coordinates with SAFE
			x1_0 = piMatchingCols[i];
			y1_0 = piMatchingRows[i];
			GetTransformedCoords(x1_0, y1_0, iTransformationType, dCoefs_i1, &x1, &y1, piK_valid_nImages[index_image1]);

			// image 2
			// get affine-transformed coordinates (in pixels) within the SAFE
			// (tile_h_offset_pixels2 + piMatchingCols[i],  tile_v_offset_pixels2 + piMatchingRows[i]) is the untransformed coordinates within SAFE
			x2_0 = piMatchingCols[i];
			y2_0 = piMatchingRows[i];
			GetTransformedCoords(x2_0, y2_0, iTransformationType, dCoefs_i2, &x2, &y2, piK_valid_nImages[index_image2]);

			// calculate residuals for matched point i
			dRes = x1 - x2 - pfShifts_x[i];
			dRes_SQR_x = dRes * dRes;
			ptr_dResiduals[2 * iMatch] += dRes_SQR_x;
			dSum += dRes_SQR_x;
			pd_Sum_per_Image[index_image1] += dRes_SQR_x; // v2
			pd_Sum_per_Image[index_image2] += dRes_SQR_x; // v2

			dRes = y1 - y2 - pfShifts_y[i];
			dRes_SQR_y = dRes * dRes;
			ptr_dResiduals[2 * iMatch + 1] += dRes_SQR_y;
			dSum += dRes_SQR_y;
			pd_Sum_per_Image[index_image1] += dRes_SQR_y; // v2
			pd_Sum_per_Image[index_image2] += dRes_SQR_y; // v2

			dMisregistrationRes_mean += sqrt(dRes_SQR_x + dRes_SQR_y);
		}

		// calcualte RMSE for x and y for current matched tile
		if (iTransformationType != 4)
			n_unknowns = n_coefs_per_image;
		else
			n_unknowns = 12 + 2*MAX(piK_valid_nImages[index_image1], piK_valid_nImages[index_image2]); // in theory, the K should be calculated for each matched pair

		// 9/15/2017, use ntie_ instead of ntie
		ptr_dResiduals[2 * iMatch] = (ntie_ > n_unknowns) ? sqrt(ptr_dResiduals[2 * iMatch] / (ntie_ - n_unknowns)) : 0;
		ptr_dResiduals[2 * iMatch + 1] = (ntie_ > n_unknowns) ? sqrt(ptr_dResiduals[2 * iMatch + 1] / (ntie_ - n_unknowns)) : 0;

		ntie_total += ntie_; // v2: use ntie_ instead of ntie
		pi_ntie_per_Image[index_image1] += ntie_;
		pi_ntie_per_Image[index_image2] += ntie_;

		dMisregistrationRes_mean /= ntie_;

		pdMisregistrationRes_means[nMatchedPairs] = dMisregistrationRes_mean;
		nMatchedPairs += 1;

		// v1: get post-registration shifts per image
		// 8/11/2021: still count if neither images are controls
//		if (piImages_Control_Flag[index_image1] > 0 || piImages_Control_Flag[index_image2] > 0) // at least one image is L1TP reference (control)
		{
			if (piImages_Control_Flag[index_image1] == 0)
				ptr_dPost_MisregistrationPerImage[index_image1] += dMisregistrationRes_mean;
			if (piImages_Control_Flag[index_image2] == 0)
				ptr_dPost_MisregistrationPerImage[index_image2] += dMisregistrationRes_mean;
		}
	}
	fclose(file_matching); //added 2/27/2017

	Std1(pdMisregistrationRes_means, nMatchedPairs, &dMisregistrationRes_means_std, &dMisregistrationRes_means_mean);
	if (ptr_dMisregistrationRes_stats != NULL)
	{
		ptr_dMisregistrationRes_stats[0] = dMisregistrationRes_means_mean;
		ptr_dMisregistrationRes_stats[1] = dMisregistrationRes_means_std;
	}

	// calculate overall RMSE
	n_unknowns_total = 0;
	for (i = 0; i < n_Images; i++)
	{
		if (piImages_Control_Flag[i] != 0)
			continue; // = 1 (control), or -1 (unconnected or unsolved)

		if (iTransformationType != 4)
			n_unknowns_total += n_coefs_per_image;
		else
			n_unknowns_total += (12 + piK_valid_nImages[i] * 2);
	}
	ptr_dRMSE = (2 * ntie_total > n_unknowns_total) ? sqrt(dSum / (2 * ntie_total - n_unknowns_total)) : 0;

	// v2
	memset(ptr_dRMSEsPerImage, 0, n_Images * sizeof(double));
	for (i = 0; i < n_Images; i++)
	{
		if (iTransformationType != 4)
			ptr_dRMSEsPerImage[i] = (2 * pi_ntie_per_Image[i] > n_coefs_per_image) ? sqrt(pd_Sum_per_Image[i] / (2 * pi_ntie_per_Image[i] - n_coefs_per_image)) : 0;
		else
		{
			n_unknowns = (12 + piK_valid_nImages[i] * 2);
			ptr_dRMSEsPerImage[i] = (2 * pi_ntie_per_Image[i] > n_unknowns) ? sqrt(pd_Sum_per_Image[i] / (2 * pi_ntie_per_Image[i] - n_unknowns)) : 0;
		}
	}

	// v1: get post-registration shifts per image
	for (i = 0; i < n_Images; i++)
	{
		if (nCount[i] > 0)
		{
			ptr_dPre_MisregistrationPerImage[i] /= nCount[i];
			ptr_dPost_MisregistrationPerImage[i] /= nCount[i];
		}
	}

	free(pfShifts_x);
	free(pfShifts_y);
	free(piMatchingCols);
	free(piMatchingRows);
	free(pdMisregistrationRes_means);

	return true;
}

/***
Created 1/30/2020
Revised from Adjustment_v1()
Use bForwardModel = true such that the output coefficients do not need to be reversed
8/2/2021
-temp: use one image as control
-do not well-registered L1TP images as controls
v1 (8/10/2021): 
-when dealing with (x1, y1) and (x2, y2), apply shifts to make (x2, y2), which is similar to those in one-to-one image matching
note (11/17/2022): not sure what the above was for; but shifts are applied to both (x1, y1) and (x2, y2), and so should not make much difference
7/25/2022
-use MAX_MSS_SHIFT_THRESHOLD_L1TP (= 4) rather than MAX_MSS_SHIFT_THRESHOLD (= 15)
-add A_flag
-if iControl_idx < 0, use L8 and L9 images as controls
***/
double Landsat_registration::Adjustment_RBFs_v1(int RBFs_K, int iControl_idx)
{
	FILE* file_matching = NULL;
	char pacMatchingFileName[STRLEN], pacMatchingFileDir[STRLEN];
	char pacMatchingSummaryFileName[STRLEN];
	int iMatchedPairsNum, iMatch;
	int iRowMax, iColMax;
	int iMatchingRowMax, iMatchingColMax;

	Landsat_image* pImage1 = NULL, * pImage2 = NULL, * pImage = NULL;
	int index_image, index_image1, index_image2;

	FILE* fout = NULL;
	char pacAdjustmentFileName[STRLEN];

	int n;		// number of unknowns
	double* U = NULL, * A = NULL, * N = NULL;
	bool* A_flag = NULL;

	int iIterNum;
	int i, k, r, m;
	double L;
	int ntie;
	double* pdResiduals = NULL;		// defined per matched tile pair, not per matched pixel
	double dRMSE;

	double* dCoefs = NULL;	// unknowns
	double* dCoefs_all = NULL;
	int n_Images, n_unknowns;
	float* pfShifts_x = NULL, * pfShifts_y = NULL;
	int* piMatchingCols = NULL, * piMatchingRows = NULL;

	bool bSuc;

	double dMeanShift_x, dMeanShift_y, RMSE;

	double pdRMSEsPerImage[MAX_PROC_LANDSAT_IMAGES];
	double pdMisregistrationRes_stats[2];
	double x1, y1, x2, y2;
	double x1_0, y1_0, x2_0, y2_0;

	double pdPre_MisregistrationPerImage[MAX_PROC_LANDSAT_IMAGES];
	double pdPost_MisregistrationPerImage[MAX_PROC_LANDSAT_IMAGES];

	// note: piCoefsIndexes_all corresponds to dCoefs_all; some SAFEs are controls or partial controls, their coefficients are not indexed or not all indexed
	bool bImage1_control, bImage2_control;
	int piImages_Control_Flag[MAX_PROC_LANDSAT_IMAGES]; // v2_2_5, change from bool to int; 0: not control; 1: control; -1: not connected

	int* piCoefsIdx = NULL;
	int CoefsIdx_idx = 0;

	// Adjustment_RBFs
	int n_coefs_per_image, n_coefs_per_image_RBFs;
	int* piCoefsIndexes_all = NULL; // indexes of unknown coefficients
	bool* pbSolvedCoefsFlag = NULL; // used in INVSQR2(); corresponds to dCoefs
	double D_k;
	double RBFs_k_1[MAX_RBFs_K], RBFs_k_2[MAX_RBFs_K];

	int K1, K2;
	int iWidth = this->m_width;
	int iHeight = this->m_height;

	unsigned int* piRBFs_ntie_Map = NULL;
	unsigned int iRBF_ties_num_threshold = 50;// RBF_TIE_NUMBER_THRESHOLD;
	double* xs_k_valid_nImages = NULL, * ys_k_valid_nImages = NULL; // valid ones
	double* xs_k_i = NULL, * ys_k_i = NULL;
	int* K_valid_nImages = NULL; // number of valid K
	double* dCoefs_i1 = NULL, * dCoefs_i2 = NULL;
	double* dCoefs_all_RBFs = NULL; // per image: 12 polynomial, 2 x K center weights, 2 x K center coordinates

	n_Images = this->m_ImagesNum;

	this->GetL1TPMisregistrationsFromScanningResults();

	n_coefs_per_image = 12 + 2 * RBFs_K;
	n_coefs_per_image_RBFs = 12 + 4 * RBFs_K;

	if (!(piCoefsIndexes_all = (int*)calloc(n_Images * n_coefs_per_image, sizeof(int)))
		|| !(pbSolvedCoefsFlag = (bool*)calloc(n_Images * n_coefs_per_image, sizeof(bool)))
		|| !(piRBFs_ntie_Map = (unsigned int*)calloc(n_Images * n_coefs_per_image, sizeof(unsigned int)))
		|| !(xs_k_valid_nImages = (double*)calloc(n_Images * RBFs_K, sizeof(double)))
		|| !(ys_k_valid_nImages = (double*)calloc(n_Images * RBFs_K, sizeof(double)))
		|| !(K_valid_nImages = (int*)calloc(n_Images, sizeof(int)))
		|| !(dCoefs_i1 = (double*)calloc(12 + RBFs_K * 4, sizeof(double)))
		|| !(dCoefs_i2 = (double*)calloc(12 + RBFs_K * 4, sizeof(double))))

	{
		printf("\nError in Adjustment_RBFs(): insufficient memory.\n");
		scanf_s(" %d", &n);
		exit(1);
	}

	sprintf_s(pacMatchingSummaryFileName, STRLEN, "%s/L1TP_matching_rsp%d.txt", this->m_cOutputDir, this->m_iDenseMatchingInterval);

	// get the numers of the valid K centers per image, and corresponding xs and ys
	FindValid_RBF_Centers_PerImage_v1(pacMatchingSummaryFileName, RBFs_K, piRBFs_ntie_Map, iRBF_ties_num_threshold, xs_k_valid_nImages, ys_k_valid_nImages, K_valid_nImages);

	// Determine controls and indexes of the unknown coefficients that are to be solved
	// initializations
	n_unknowns = 0; // number of unknown coefficients
	memset(piImages_Control_Flag, 0, MAX_PROC_LANDSAT_IMAGES * sizeof(int));
	CoefsIdx_idx = 0;
	for (i = 0; i < n_Images * n_coefs_per_image; i++)
		piCoefsIndexes_all[i] = -1;

	// set control flags, count unknowns, and set valid unknown coefficients' indexes in piCoefsIndexes_all(invalid set as -1) 
	// piImages_Control_Flag: -1 = unconnected, 1 = control, 0 = not control
	for (index_image = 0; index_image < n_Images; index_image++)
	{
		if (index_image == iControl_idx)
		{
			piImages_Control_Flag[index_image] = 1;
			continue;
		}
		if (iControl_idx < 0)
		{
			pImage1 = GetImage(index_image);
			if (pImage1->sensor == enum_L8 || pImage1->sensor == enum_L9)
			{
				piImages_Control_Flag[index_image] = 1;
				continue;
			}
		}

		// not control, to be solved
		// 12 2nd-order polynomial coefficients
		piImages_Control_Flag[index_image] = 0;
		for (i = 0; i < 12; i++)
		{
			piCoefsIndexes_all[index_image * n_coefs_per_image + i] = CoefsIdx_idx;
			CoefsIdx_idx += 1;
			n_unknowns += 1;
		}
		// RBF centers per image
		for (m = 0; m < K_valid_nImages[index_image]; m++)
		{
			piCoefsIndexes_all[index_image * n_coefs_per_image + (m * 2 + 12)] = CoefsIdx_idx;
			CoefsIdx_idx += 1;
			piCoefsIndexes_all[index_image * n_coefs_per_image + (m * 2 + 12 + 1)] = CoefsIdx_idx;
			CoefsIdx_idx += 1;
			n_unknowns += 2;
		}
	}

	if (LargerOrEqualTo_Num_INT(piImages_Control_Flag, n_Images, 1) <= 0) // added 8/30/2019
	{
		printf("\nError in Adjustment_v1(): no controls.\n");
		scanf_s(" %d", &n);
		exit(1);
	}

	iColMax = m_width;
	iRowMax = m_height;
	iMatchingColMax = iColMax / this->m_iDenseMatchingInterval;
	iMatchingRowMax = iRowMax / this->m_iDenseMatchingInterval;

	n = n_unknowns;
	iMatchedPairsNum = Gettxtlines(pacMatchingSummaryFileName);	// number of matched tiles
	// allocate memory
	if (!(N = (double*)calloc(n * n, sizeof(double)))
		|| !(U = (double*)calloc(n, sizeof(double)))
		|| !(A = (double*)calloc(n, sizeof(double)))
		|| !(A_flag = (bool*)calloc(n, sizeof(bool)))
		|| !(dCoefs = (double*)calloc(n, sizeof(double)))
		|| !(dCoefs_all = (double*)calloc(n_Images * n_coefs_per_image, sizeof(double)))
		|| !(dCoefs_all_RBFs = (double*)calloc(n_Images * n_coefs_per_image_RBFs, sizeof(double)))
		|| !(pfShifts_x = (float*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(float)))
		|| !(pfShifts_y = (float*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(float)))
		|| !(piMatchingCols = (int*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(int)))
		|| !(piMatchingRows = (int*)calloc(iMatchingRowMax * iMatchingColMax, sizeof(int)))
		|| !(pdResiduals = (double*)calloc(2 * iMatchedPairsNum, sizeof(double))))
	{
		printf("\nError in Adjustment_v1(): insufficient memory.\n");
		scanf_s(" %d", &n);
		exit(1);
	}

	// set initial affine transformation coefficients
	memset(dCoefs_all, 0, n_Images * n_coefs_per_image*sizeof(double));
	for (i = 0; i < n_Images; i++)
	{
		dCoefs_all[i * n_coefs_per_image + 1] = 1;
		dCoefs_all[i * n_coefs_per_image + 8] = 1;
	}

	// set unknown coefficients' approximates
	CoefsIdx_idx = 0;
	for (i = 0; i < n_Images * n_coefs_per_image; i++)
	{
		if (piCoefsIndexes_all[i] >= 0)
		{
			dCoefs[CoefsIdx_idx] = dCoefs_all[i];
			CoefsIdx_idx += 1;
		}
	}

	// Start adjustment
	iIterNum = 0;
	while (iIterNum < 1)
	{
		// list observation equations
		file_matching = Readtxt(pacMatchingSummaryFileName);

		memset(U, 0, n * sizeof(double));
		memset(N, 0, n * n * sizeof(double));
		for (iMatch = 0; iMatch < iMatchedPairsNum; iMatch++)
		{
			printf("matched pair %d / %d\n", iMatch + 1, iMatchedPairsNum);
			// note: dMeanShift_x and dMeanShift_y are fitted translation tranformation coefficients a0 and b0
			fscanf_s(file_matching, "%d%d%lf%lf%lf%d", &index_image1, &index_image2, &dMeanShift_x, &dMeanShift_y, &RMSE, &ntie);
			fscanf_s(file_matching, "%s", pacMatchingFileDir, STRLEN);

			index_image1 = index_image1 - 1; // L1TP_matching.txt is indexed starting from 1; L1TP_L1GS_matching is indexed from 0
			index_image2 = index_image2 - 1;

			// v2_2_5
			if (piImages_Control_Flag[index_image1] < 0 || piImages_Control_Flag[index_image2] < 0)
				continue; // not connected, skip

			pImage1 = GetImage(index_image1);
			pImage2 = GetImage(index_image2);

			bImage1_control = piImages_Control_Flag[index_image1] == 1;
			bImage2_control = piImages_Control_Flag[index_image2] == 1;

			if (bImage1_control && bImage2_control)
				continue;

			// v2_2_5, use the same criteria as those used in FindUnconnectedSAFEs()
			if (IsGoodMatching_MSS(true, RMSE, dMeanShift_x, dMeanShift_y) == false)
				continue;

			// get matched point coordinates in the tile image (piMatchingCols, piMatchingRows) and correspoinding shift x, shift y
			sprintf_s(pacMatchingFileName, STRLEN, "%s/matching/layer0_dense_matching_rsp%d", pacMatchingFileDir, this->m_iDenseMatchingInterval);
			ntie = GetShifts_v1(pacMatchingFileName, iMatchingColMax, iMatchingRowMax, piMatchingCols, piMatchingRows, this->m_iDenseMatchingInterval, 0.95f, pfShifts_x, pfShifts_y);

			K1 = K_valid_nImages[index_image1];
			K2 = K_valid_nImages[index_image2];

			// list observation equations per matched point
			for (i = 0; i < ntie; i++)
			{
				// translation-fitting RMSE based check, added 2/21/2017, see FitTranslationTransform()
//				if (sqrt((dMeanShift_x + pfShifts_x[i]) * (dMeanShift_x + pfShifts_x[i]) + (dMeanShift_y + pfShifts_y[i]) * (dMeanShift_y + pfShifts_y[i])) > 2 * RMSE)
//					continue;

				// image 1
				x1_0 = piMatchingCols[i] + pfShifts_x[i]; // v1: apply shifts here
				y1_0 = piMatchingRows[i] + pfShifts_y[i];

				if (bImage1_control == false)
				{
					memmove(dCoefs_i1, dCoefs_all + index_image1 * n_coefs_per_image, 12 * sizeof(double));
					for (m = 0; m < K1; m++)
					{
						dCoefs_i1[12 + 4 * m] = dCoefs_all[index_image1 * n_coefs_per_image + 12 + m * 2]; // wx_k[m]
						dCoefs_i1[12 + 4 * m + 1] = dCoefs_all[index_image1 * n_coefs_per_image + 12 + m * 2 + 1]; // wy_k[m]
						dCoefs_i1[12 + 4 * m + 2] = xs_k_valid_nImages[index_image1 * RBFs_K + m]; // xs_k[m]
						dCoefs_i1[12 + 4 * m + 3] = ys_k_valid_nImages[index_image1 * RBFs_K + m]; // ys_k[m]
					}
					GetTransformedCoords(x1_0, y1_0, 4, dCoefs_i1, &x1, &y1, K1);

					// calcualte squared distances to kernal points and RBFs function values
					xs_k_i = xs_k_valid_nImages + index_image1 * RBFs_K;
					ys_k_i = ys_k_valid_nImages + index_image1 * RBFs_K;
					for (m = 0; m < K1; m++)
					{
						D_k = (xs_k_i[m] - x1_0) * (xs_k_i[m] - x1_0) + (ys_k_i[m] - y1_0) * (ys_k_i[m] - y1_0);
						RBFs_k_1[m] = D_k * log(sqrt(D_k));// exp(-D_k);
					}
				}
				else
				{
					x1 = x1_0;
					y1 = y1_0;
				}

				// image 2
				x2_0 = piMatchingCols[i] + pfShifts_x[i]; // v1: apply shifts here
				y2_0 = piMatchingRows[i] + pfShifts_y[i];

				if (bImage2_control == false)
				{
					memmove(dCoefs_i2, dCoefs_all + index_image2 * n_coefs_per_image, 12 * sizeof(double));
					for (m = 0; m < K2; m++)
					{
						dCoefs_i2[12 + 4 * m] = dCoefs_all[index_image2 * n_coefs_per_image + 12 + m * 2]; // wx_k[m]
						dCoefs_i2[12 + 4 * m + 1] = dCoefs_all[index_image2 * n_coefs_per_image + 12 + m * 2 + 1]; // wy_k[m]
						dCoefs_i2[12 + 4 * m + 2] = xs_k_valid_nImages[index_image2 * RBFs_K + m]; // xs_k[m]
						dCoefs_i2[12 + 4 * m + 3] = ys_k_valid_nImages[index_image2 * RBFs_K + m]; // ys_k[m]
					}
					GetTransformedCoords(x2_0, y2_0, 4, dCoefs_i2, &x2, &y2, K2);

					// calcualte squared distances to kernal points and RBFs function values
					xs_k_i = xs_k_valid_nImages + index_image2 * RBFs_K;
					ys_k_i = ys_k_valid_nImages + index_image2 * RBFs_K;
					for (m = 0; m < K2; m++)
					{
						D_k = (xs_k_i[m] - x2_0) * (xs_k_i[m] - x2_0) + (ys_k_i[m] - y2_0) * (ys_k_i[m] - y2_0);
						RBFs_k_2[m] = D_k * log(sqrt(D_k));// exp(-D_k);
					}
				}
				else
				{
					x2 = x2_0;
					y2 = y2_0;
				}

				L = pfShifts_x[i] - (x1 - x2);
				// note: L is shift x; when the mean shift x is close to MAX_SENTINEL2A_SHIFT_THRESHOLD, individual shift x value can be larger than MAX_SENTINEL2A_SHIFT_THRESHOLD;
				//  this condition checking on pioint i's shift value is not valid because of similar check is already done above (2*RMSE check); here it mainly works to check whether the 
				//  calculations of coord_h1_UTM and coord_h2_UTM are correct
				memset(A, 0, n * sizeof(double));
				memset(A_flag, 0, n * sizeof(bool));
				if (bImage1_control == false)
				{
					piCoefsIdx = piCoefsIndexes_all + index_image1 * n_coefs_per_image;

					// x2 = a0 + a1*x1 + a2*y1 + wx1*exp(-Di1^2) + wx2*exp(-Di2^2) + wx3*exp(-Di3^2) + ...
					// Di1^2 = (x1 - xk1)^2 + (y1 - yk1)^2

					// v1: change (x1, y1) to (x1_0, y1_0)
					if (piCoefsIdx[0] >= 0)
						A[piCoefsIdx[0]] = 1;
					if (piCoefsIdx[1] >= 0)
						A[piCoefsIdx[1]] = x1_0;
					if (piCoefsIdx[2] >= 0)
						A[piCoefsIdx[2]] = y1_0;
					if (piCoefsIdx[3] >= 0)
						A[piCoefsIdx[3]] = x1_0 * x1_0;
					if (piCoefsIdx[4] >= 0)
						A[piCoefsIdx[4]] = x1_0 * y1_0;
					if (piCoefsIdx[5] >= 0)
						A[piCoefsIdx[5]] = y1_0 * y1_0;
					// A[6] ~ A[11] = 0
					for (m = 0; m < 6; m++)
						A_flag[piCoefsIdx[m]] = true;  // added 7/25/2022

					for (m = 0; m < K1; m++)
					{
						if (piCoefsIdx[12 + 2 * m] >= 0)
						{
							A[piCoefsIdx[12 + 2 * m]] = RBFs_k_1[m]; // df / d(wxk)
							A_flag[piCoefsIdx[12 + 2 * m]] = true; // added 7/25/2022
						}
					}
				}
				if (bImage2_control == false)
				{
					// v1: change (x2, y2) to (x2_0, y2_0)
					piCoefsIdx = piCoefsIndexes_all + index_image2 * n_coefs_per_image;
					if (piCoefsIdx[0] >= 0)
						A[piCoefsIdx[0]] = -1;
					if (piCoefsIdx[1] >= 0)
						A[piCoefsIdx[1]] = -x2_0;
					if (piCoefsIdx[2] >= 0)
						A[piCoefsIdx[2]] = -y2_0;
					if (piCoefsIdx[3] >= 0)
						A[piCoefsIdx[3]] = -x2_0 * x2_0;
					if (piCoefsIdx[4] >= 0)
						A[piCoefsIdx[4]] = -x2_0 * y2_0;
					if (piCoefsIdx[5] >= 0)
						A[piCoefsIdx[5]] = -y2_0 * y2_0;
					// A[6] ~ A[11] = 0
					for (m = 0; m < 6; m++)
						A_flag[piCoefsIdx[m]] = true;  // added 7/25/2022

					for (m = 0; m < K2; m++)
					{
						if (piCoefsIdx[12 + 2 * m] >= 0)
						{
							A[piCoefsIdx[12 + 2 * m]] = -RBFs_k_2[m]; // df / d(wxk)
							A_flag[piCoefsIdx[12 + 2 * m]] = true; // added 7/25/2022
						}
					}
				}

				for (k = 0; k < n; k++)
				{
					//if (ABS(A[k]) <= 1e-10)
					//	continue;
					if (!A_flag[k])
						continue;

					for (r = 0; r < n; r++)
					{
						//N[k * n + r] += (ABS(A[r]) > 1e-10) ? A[k] * A[r] : 0;
						//if (A[r] > 1e-10 || A[r] < -(1e-10))
						if (A_flag[r])
							N[k * n + r] += A[k] * A[r];
					}
					U[k] += A[k] * L;
				}

				L = pfShifts_y[i] - (y1 - y2);
				memset(A, 0, n * sizeof(double));
				memset(A_flag, 0, n * sizeof(bool));
				if (bImage1_control == false)
				{
					piCoefsIdx = piCoefsIndexes_all + index_image1 * n_coefs_per_image;
					if (piCoefsIdx[6] >= 0)
						A[piCoefsIdx[6]] = 1;
					if (piCoefsIdx[7] >= 0)
						A[piCoefsIdx[7]] = x1_0;
					if (piCoefsIdx[8] >= 0)
						A[piCoefsIdx[8]] = y1_0;
					if (piCoefsIdx[9] >= 0)
						A[piCoefsIdx[9]] = x1_0 * x1_0;
					if (piCoefsIdx[10] >= 0)
						A[piCoefsIdx[10]] = x1_0 * y1_0;
					if (piCoefsIdx[11] >= 0)
						A[piCoefsIdx[11]] = y1_0 * y1_0;
					// A[0] ~ A[5] = 0
					for (m = 6; m < 12; m++)
						A_flag[piCoefsIdx[m]] = true;  // added 7/25/2022

					for (m = 0; m < K1; m++)
					{
						if (piCoefsIdx[12 + 2 * m + 1] >= 0)
						{
							A[piCoefsIdx[12 + 2 * m + 1]] = RBFs_k_1[m]; // df / d(wxk)
							A_flag[piCoefsIdx[12 + 2 * m + 1]] = true; // added 7/25/2022
						}
					}
				}
				if (bImage2_control == false)
				{
					piCoefsIdx = piCoefsIndexes_all + index_image2 * n_coefs_per_image;
					if (piCoefsIdx[6] >= 0)
						A[piCoefsIdx[6]] = -1;
					if (piCoefsIdx[7] >= 0)
						A[piCoefsIdx[7]] = -x2_0;
					if (piCoefsIdx[8] >= 0)
						A[piCoefsIdx[8]] = -y2_0;
					if (piCoefsIdx[9] >= 0)
						A[piCoefsIdx[9]] = -x2_0 * x2_0;
					if (piCoefsIdx[10] >= 0)
						A[piCoefsIdx[10]] = -x2_0 * y2_0;
					if (piCoefsIdx[11] >= 0)
						A[piCoefsIdx[11]] = -y2_0 * y2_0;
					// A[0] ~ A[5] = 0
					for (m = 6; m < 12; m++)
						A_flag[piCoefsIdx[m]] = true;  // added 7/25/2022

					for (m = 0; m < K2; m++)
					{
						if (piCoefsIdx[12 + 2 * m + 1] >= 0)
						{
							A[piCoefsIdx[12 + 2 * m + 1]] = -RBFs_k_2[m]; // df / d(wxk)
							A_flag[piCoefsIdx[12 + 2 * m + 1]] = true; // added 7/25/2022
						}
					}
				}
				for (k = 0; k < n; k++)
				{
					//if (ABS(A[k]) <= 1e-10)
					//	continue;
					if (!A_flag[k])
						continue;

					for (r = 0; r < n; r++)
					{
						//	N[k * n + r] += (ABS(A[r]) > 1e-10) ? A[k] * A[r] : 0;
						//if (A[r] > 1e-10 || A[r] < -(1e-10))
						if (A_flag[r])
							N[k * n + r] += A[k] * A[r];
					}
					U[k] += A[k] * L;
				}
			}
		} // end for (iMatch = 0; iMatch < iMatchedPairsNum; iMatch++)
		fclose(file_matching);

		if (!INVSQR2(N, U, n, pbSolvedCoefsFlag))
		{
			printf(" Error in Landsat_registration::Adjustment_v1(): normal matrix is singular\n");
			scanf_s(" %d", &n);
		}

		// update coefficients
		for (k = 0; k < n; k++)
			dCoefs[k] += U[k];

		bSuc = true;
		for (k = 0; k < n; k++)
		{
			if (ABS(U[k]) > DELTA_LIMIT)
			{
				bSuc = false;
				break;
			}
		}
		if (bSuc == true)
			break;	// solution is accuracy enough

		// update unknown coefficients in dCoefs_all
		CoefsIdx_idx = 0;
		for (i = 0; i < n_Images * n_coefs_per_image; i++)
		{
			if (piCoefsIndexes_all[i] >= 0)
			{
				dCoefs_all[i] = (pbSolvedCoefsFlag[CoefsIdx_idx]) ? dCoefs[CoefsIdx_idx] : -99999; // v2_2_5, only update solved coefficients
				if (pbSolvedCoefsFlag[CoefsIdx_idx] == false)
				{
					if (piImages_Control_Flag[i / n_coefs_per_image] == 1)
						printf("\nUnsolvable coefficients: Image %d\n", i / n_coefs_per_image); // print out unsolved image
					piImages_Control_Flag[i / n_coefs_per_image] = -1; // mark image (i/12) as unconnected
				}
				CoefsIdx_idx += 1;
			}
		}

		iIterNum += 1;
	}

	// put all coefficients into dCoefs_all
	CoefsIdx_idx = 0;
	for (i = 0; i < n_Images * n_coefs_per_image; i++)
	{
		index_image = i / n_coefs_per_image;
		if (piCoefsIndexes_all[i] >= 0)
		{
			dCoefs_all[i] = (piImages_Control_Flag[index_image] == 0) ? dCoefs[CoefsIdx_idx] : 0; // piImages_Control_Flag[] == 0: valid; == -1, unconnected or unsolved
			CoefsIdx_idx += 1;
		}
	}

	// put all coefficients into dCoefs_all_RBFs
	memset(dCoefs_all_RBFs, 0, n_Images * n_coefs_per_image_RBFs * sizeof(double));
	for (index_image = 0; index_image < n_Images; index_image++)
	{
		memcpy(dCoefs_all_RBFs + n_coefs_per_image_RBFs * index_image, dCoefs_all + n_coefs_per_image * index_image, 12 * sizeof(double));
		for (m = 0; m < K_valid_nImages[index_image]; m++)
		{
			dCoefs_all_RBFs[index_image * n_coefs_per_image_RBFs + 12 + 4 * m] = dCoefs_all[index_image * n_coefs_per_image + 12 + m * 2]; // wx_k[m]
			dCoefs_all_RBFs[index_image * n_coefs_per_image_RBFs + 12 + 4 * m + 1] = dCoefs_all[index_image * n_coefs_per_image + 12 + m * 2 + 1]; // wy_k[m]
			dCoefs_all_RBFs[index_image * n_coefs_per_image_RBFs + 12 + 4 * m + 2] = xs_k_valid_nImages[index_image * RBFs_K + m]; // xs_k[m]
			dCoefs_all_RBFs[index_image * n_coefs_per_image_RBFs + 12 + 4 * m + 3] = ys_k_valid_nImages[index_image * RBFs_K + m]; // ys_k[m]
		}
	}

	// calculate residuals per matched tile and overall RMSE
	this->CalculateAdjustmentResiduals_v2(dCoefs_all_RBFs, piImages_Control_Flag, 4, K_valid_nImages, RBFs_K, pdResiduals, dRMSE, pdRMSEsPerImage, pdMisregistrationRes_stats, pdPre_MisregistrationPerImage, pdPost_MisregistrationPerImage, false);
	printf("Mean pairwise misregistration = %.5f; std = %.5f\n\n", pdMisregistrationRes_stats[0], pdMisregistrationRes_stats[1]);

	// output adjusted coordinates
	sprintf_s(pacAdjustmentFileName, STRLEN, "%s/adjustment_RBFs_K%d_%dimages_ctrl_%d.txt", this->m_cOutputDir, RBFs_K, n_Images, iControl_idx);

	fout = Writetxt(pacAdjustmentFileName);
	for (i = 0; i < n_Images; i++)
	{
		fprintf(fout, "%d\t", i);
		if (piImages_Control_Flag[i] < 0)
		{
			// unregistered
			for (k = 0; k < n_coefs_per_image_RBFs; k++)
				fprintf(fout, "-99999 ");
			fprintf(fout, "\n");
			continue;
		}

		for (k = 0; k < n_coefs_per_image_RBFs; k++)
			fprintf(fout, "%.16f ", dCoefs_all_RBFs[i * n_coefs_per_image_RBFs + k]);

		fprintf(fout, "\n");
	}

	// output residuals and RMSE
	fprintf(fout, "\n");
	file_matching = Readtxt(pacMatchingSummaryFileName);
	for (iMatch = 0; iMatch < iMatchedPairsNum; iMatch++)
	{
		fscanf_s(file_matching, "%d%d%lf%lf%lf%d", &index_image1, &index_image2, &dMeanShift_x, &dMeanShift_y, &RMSE, &ntie);
		fscanf_s(file_matching, "%s", pacMatchingFileDir, STRLEN);

		fprintf(fout, "%.6f %.6f\t%d\t%d\t%s\n", pdResiduals[2 * iMatch], pdResiduals[2 * iMatch + 1], index_image1, index_image2, pacMatchingFileDir);
	}
	fprintf(fout, "\nRMSE: %.6f\n\n", dRMSE);
	fclose(file_matching);

	for (i = 0; i < this->m_ImagesNum; i++)
		fprintf(fout, "%d\t%.6f\t%.6f\t%.6f\n", this->GetImage(i)->ID_date, pdRMSEsPerImage[i], pdPre_MisregistrationPerImage[i], pdPost_MisregistrationPerImage[i]);
	fclose(fout);

	free(N);
	free(U);
	free(A);
	free(A_flag);
	free(dCoefs);
	free(dCoefs_all);
	free(pdResiduals);
	free(pfShifts_x);
	free(pfShifts_y);
	free(piMatchingCols);
	free(piMatchingRows);
	free(piCoefsIndexes_all);
	free(pbSolvedCoefsFlag);
	free(piRBFs_ntie_Map);
	free(xs_k_valid_nImages);
	free(ys_k_valid_nImages);
	free(K_valid_nImages);
	free(dCoefs_i1);
	free(dCoefs_i2);
	free(dCoefs_all_RBFs);

	return dRMSE;
}

bool Landsat_registration::InputAdjustedCoefficients(double* ptr_coefs)
{
	char pacAdjustmentFileName[STRLEN];
	FILE* fin = NULL;
	int i, idx, index_image;
	int n_Images;

	n_Images = this->m_ImagesNum;

	sprintf_s(pacAdjustmentFileName, STRLEN, "%s/adjustment_AF_v2_2_5.txt", this->m_cOutputDir);
	
	CheckFileExist(pacAdjustmentFileName);

	fin = Readtxt(pacAdjustmentFileName);
	for (index_image = 0; index_image < n_Images; index_image++)
	{
		fscanf_s(fin, "%d", &idx);
		
		for (i = 0; i < 12; i++)
			fscanf_s(fin, "%lf", ptr_coefs + idx * 12 + i);
	}
	fclose(fin);

	return true;
}

bool Landsat_registration::InputAdjustedCoefficients_RBFs(double* ptr_coefs, int RBFs_K, int iControl_idx)
{
	char pacAdjustmentFileName[STRLEN];
	FILE* fin = NULL;
	int i, idx, index_image;
	int n_Images;
	int n_coefs_per_image_RBFs;

	n_Images = this->m_ImagesNum;

	sprintf_s(pacAdjustmentFileName, STRLEN, "%s/adjustment_RBFs_K%d_%dimages_ctrl_%d.txt", this->m_cOutputDir, RBFs_K, n_Images, iControl_idx);

	CheckFileExist(pacAdjustmentFileName);

	fin = Readtxt(pacAdjustmentFileName);
	n_coefs_per_image_RBFs = 12 + 4 * RBFs_K;
	for (index_image = 0; index_image < n_Images; index_image++)
	{
		fscanf_s(fin, "%d", &idx);

		for (i = 0; i < n_coefs_per_image_RBFs; i++)
			fscanf_s(fin, "%lf", ptr_coefs + idx * n_coefs_per_image_RBFs + i);
	}
	fclose(fin);

	return true;
}

/***
Revised from Sentinel_UTM_Registration::OutputRegisteredTileImagesStack_Polynomial()
Note: the input is 16-bit short and output is 8-bit unsigned char
***/
bool Landsat_registration::OutputRegisteredImagesStack(bool bOutputOriginal)
{
	double pdCoefs_all[MAX_PROC_LANDSAT_IMAGES * 12], pdCoefs[12];
	short int* pshtImageData = NULL;
	unsigned char *pucImageData_reg = NULL;
	int iColMax, iRowMax;
	int i;
	FILE* fout = NULL, *fout_original = NULL;
	char cBandNames[STRLEN_LONG], cBandNames_org[STRLEN_LONG];
	int iLayersNum;
	int n_Images;
	char cOutputFileName_org[STRLEN];
	Landsat_image* pImage = NULL;
	char cLandsatInfo[STRLEN];
	char cOutputFileName[STRLEN];
	char pac_null[2] = "";

	n_Images = this->m_ImagesNum;

	GetImageSize(&iColMax, &iRowMax);

	InputAdjustedCoefficients(pdCoefs_all);

	// set output file name
	pImage = this->GetImage(0);
	sprintf_s(cOutputFileName, STRLEN, "%s/%03d%03d_%04d_reg", this->m_cOutputDir, pImage->path, pImage->row, pImage->year);
	fout = WriteBinary(cOutputFileName);

	// open optional original-data output file
	if (bOutputOriginal)
	{
		cOutputFileName_org[0] = '\0';
		sprintf_s(cOutputFileName_org, STRLEN, "%s/%03d%03d_%04d_org", this->m_cOutputDir, pImage->path, pImage->row, pImage->year);
		fout_original = WriteBinary(cOutputFileName_org);
	}

	if (!(pucImageData_reg = (unsigned char*)calloc(iRowMax * iColMax, sizeof(unsigned char))))
	{
		printf("\nError in Sentinel_UTM_Registration::OutputRegisteredTileImagesStack(): insufficient memory.\n");
		scanf_s(" %d", &i);
		return false;
	}

	iLayersNum = 0;
	cBandNames[0] = '\0';
	cBandNames_org[0] = '\0';
	for (i = 0; i < n_Images; i++)
	{
		pImage = this->GetImage(i);
		if (pImage == NULL)
			continue;

		// input image data
		pshtImageData = (short int*)ReadInputImage(pImage->cImageFilePathName, iColMax, iRowMax, 1);

		// set affine transformation coefficients for the tile image
		memset(pdCoefs, 0, 12 * sizeof(double));
		memmove(pdCoefs, pdCoefs_all + i * 12, 12 * sizeof(double));

		sprintf_s(cLandsatInfo, STRLEN, "%s", pImage->cImageFileName);
		cLandsatInfo[25] = '\0';

		// do image transformation
		if ((ABS(pdCoefs[0]) > 1e-10 || ABS(pdCoefs[6]) > 1e-10) && ABS(pdCoefs[0]) < 9990)
		{
			// registered L1GS or L1TP
			TransformImage_sht2char(pshtImageData, iColMax, iRowMax, 3, pdCoefs, pucImageData_reg);
			sprintf_s(cLandsatInfo, STRLEN, "%s_reg", cLandsatInfo);
		}
		else if (pImage->proc_level == enum_L1TP)
		{
			// output original L1TP
			CopySht2UChar(pshtImageData, iRowMax * iColMax, pucImageData_reg);
			if (pdCoefs[0] < -9990)
				sprintf_s(cLandsatInfo, STRLEN, "%s_unconnected", cLandsatInfo); // unconnected L1TP, still output
		}
		else
			continue; // unregistered L1GS

		// add output band name
		sprintf_s(cBandNames, STRLEN_LONG, "%s,%d: %s", cBandNames, i, cLandsatInfo);
		
		// write data to output file
		fwrite(pucImageData_reg, sizeof(unsigned char), iRowMax * iColMax, fout);
		if (bOutputOriginal)
		{
			fwrite(pshtImageData, sizeof(short), iRowMax * iColMax, fout_original);

			sprintf_s(cLandsatInfo, STRLEN, "%s", pImage->cImageFileName);
			cLandsatInfo[25] = '\0';
			sprintf_s(cBandNames_org, STRLEN_LONG, "%s,%d: %s", cBandNames_org, i, cLandsatInfo);
		}

		iLayersNum += 1;

		free(pshtImageData);
		pshtImageData = NULL;
	}
	fclose(fout);

	// output ENVI hdr file
	OutputEnviHDRFile(cOutputFileName, pac_null, iColMax, iRowMax, iLayersNum, 1, cBandNames);

	if (bOutputOriginal)
	{
		fclose(fout_original);
		OutputEnviHDRFile(cOutputFileName_org, pac_null, iColMax, iRowMax, iLayersNum, 2, cBandNames_org);
	}

	free(pucImageData_reg);

	return true;
}

/***
Revised from Sentinel_UTM_Registration::OutputRegisteredTileImagesStack_Polynomial()
Note: the input is 16-bit short and output is 8-bit unsigned char
Created 8/10/2021
RBFs_K can be 0
***/
bool Landsat_registration::OutputRegisteredImagesStack_RBFs(bool bOutputOriginal, int RBFs_K, int iControl_idx)
{
	unsigned char* pucImageData = NULL;
	unsigned char* pucImageData_reg = NULL;
	int iColMax, iRowMax;
	int i;
	FILE* fout = NULL, * fout_original = NULL;
	char cBandNames[STRLEN_LONG], cBandNames_org[STRLEN_LONG];
	int iLayersNum;
	int n_Images;
	char cOutputFileName_org[STRLEN];
	Landsat_image* pImage = NULL;
	char cLandsatInfo[STRLEN];
	char cOutputFileName[STRLEN];
	char pac_null[2] = "";
	double* dCoefs_all_RBFs = NULL; // per image: 12 polynomial, 2 x K center weights, 2 x K center coordinates
	double* dCoefs_RBFs = NULL;

	int n_coefs_per_image_RBFs;

	n_Images = this->m_ImagesNum;
	n_coefs_per_image_RBFs = 12 + 4 * RBFs_K;
	GetImageSize(&iColMax, &iRowMax);

	if (!(dCoefs_RBFs = (double*)calloc(n_coefs_per_image_RBFs, sizeof(double)))
		|| !(dCoefs_all_RBFs = (double*)calloc(n_Images * n_coefs_per_image_RBFs, sizeof(double))))
	{
		printf("\nError in OutputRegisteredImagesStack_RBFs(): insufficient memory.\n");
		scanf_s(" %d", &i);
		exit(1);
	}

	InputAdjustedCoefficients_RBFs(dCoefs_all_RBFs, RBFs_K, iControl_idx);

	// set output file name
	pImage = this->GetImage(0);
	sprintf_s(cOutputFileName, STRLEN, "%s/%03d%03d_reg_BRFs_K%d_%dimages_ctrl_%d", this->m_cOutputDir, pImage->path, pImage->row, RBFs_K, n_Images, iControl_idx);
	fout = WriteBinary(cOutputFileName);

	// open optional original-data output file
	if (bOutputOriginal)
	{
		cOutputFileName_org[0] = '\0';
		sprintf_s(cOutputFileName_org, STRLEN, "%s/%03d%03d_org_%dimages", this->m_cOutputDir, pImage->path, pImage->row, n_Images);
		fout_original = WriteBinary(cOutputFileName_org);
	}

	if (!(pucImageData_reg = (unsigned char*)calloc(iRowMax * iColMax, sizeof(unsigned char))))
	{
		printf("\nError in Sentinel_UTM_Registration::OutputRegisteredTileImagesStack(): insufficient memory.\n");
		scanf_s(" %d", &i);
		return false;
	}

	iLayersNum = 0;
	cBandNames[0] = '\0';
	cBandNames_org[0] = '\0';
	for (i = 0; i < n_Images; i++)
	{
		pImage = this->GetImage(i);
		if (pImage == NULL)
			continue;

		printf("output image %d\n", i + 1);

		// input image data
		pucImageData = (unsigned char*)ReadInputImage(pImage->cImageFilePathName, iColMax, iRowMax, 0);

		// set affine transformation coefficients for the tile image
		memset(dCoefs_RBFs, 0, n_coefs_per_image_RBFs * sizeof(double));
		memmove(dCoefs_RBFs, dCoefs_all_RBFs + i * n_coefs_per_image_RBFs, n_coefs_per_image_RBFs * sizeof(double));

		sprintf_s(cLandsatInfo, STRLEN, "%s", pImage->cImageFileName);
		cLandsatInfo[25] = '\0';

		// do image transformation
		if ((ABS(dCoefs_RBFs[0]) > 1e-10 || ABS(dCoefs_RBFs[6]) > 1e-10) && ABS(dCoefs_RBFs[0]) < 9990)
		{
			// registered L1GS or L1TP
			TransformImage_UC(pucImageData, iColMax, iRowMax, 4, dCoefs_RBFs, pucImageData_reg, RBFs_K);
			sprintf_s(cLandsatInfo, STRLEN, "%s_reg", cLandsatInfo);
		}
		else if (pImage->proc_level == enum_L1TP)
		{
			// output original L1TP
			memmove(pucImageData_reg, pucImageData, iRowMax * iColMax*sizeof(unsigned char));
			if (dCoefs_RBFs[0] < -9990)
				sprintf_s(cLandsatInfo, STRLEN, "%s_unconnected", cLandsatInfo); // unconnected L1TP, still output
		}
		else
			continue; // unregistered L1GS

		// add output band name
		sprintf_s(cBandNames, STRLEN_LONG, "%s,%d: %s", cBandNames, i, cLandsatInfo);

		// write data to output file
		fwrite(pucImageData_reg, sizeof(unsigned char), iRowMax * iColMax, fout);
		if (bOutputOriginal)
		{
			fwrite(pucImageData, sizeof(unsigned char), iRowMax * iColMax, fout_original);

			sprintf_s(cLandsatInfo, STRLEN, "%s", pImage->cImageFileName);
			cLandsatInfo[25] = '\0';
			sprintf_s(cBandNames_org, STRLEN_LONG, "%s,%d: %s", cBandNames_org, i, cLandsatInfo);
		}

		iLayersNum += 1;

		free(pucImageData);
		pucImageData = NULL;
	}
	fclose(fout);
	printf("output file: %s\n", cOutputFileName);

	// output ENVI hdr file
	OutputEnviHDRFile(cOutputFileName, pac_null, iColMax, iRowMax, iLayersNum, 1, cBandNames);

	if (bOutputOriginal)
	{
		fclose(fout_original);
		OutputEnviHDRFile(cOutputFileName_org, pac_null, iColMax, iRowMax, iLayersNum, 1, cBandNames_org);
		printf("output file: %s\n", cOutputFileName_org);
	}

	free(pucImageData_reg);
	free(dCoefs_all_RBFs);
	free(dCoefs_RBFs);

	return true;
}