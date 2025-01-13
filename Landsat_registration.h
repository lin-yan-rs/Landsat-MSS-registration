#ifndef LANDSAT_REGISTRATION_H
#define LANDSAT_REGISTRATION_H

#include "functions.h"

#define MAX_PROC_LANDSAT_IMAGES (60) 
#define LANDSAT_MSS_MISALIGN_THRESHOLD (1)
#define DENSE_MATCHING_SAMPLING_INTERVAL_MSS (50)
#define MAX_SENTINEL2A_SHIFT_THRESHOLD (6)		// experience value, used to filter out bad matching in adjustment, maximum observed value is about 3.1; 2/26/2017, changed to 6 because over 5 pixels shifts were found in matching D:\Africa_S2\UTM_33\outputs\matching\20160829090552048-20160816090022008_LTD
#define MAX_MSS_SHIFT_THRESHOLD (15)		// experience value, used to filter out bad matching in adjustment
#define MAX_MSS_SHIFT_THRESHOLD_L1TP (4)		// experience value, used to filter out bad matching in adjustment
#define DENSE_MATCHING_SAM_THRESHOLD_MSS (0.995f)
#define MAX_DOY_DIFF (180)

typedef enum { enum_unknown = 0, enum_L1TP = 1, enum_L1GT = 2, enum_L1GS = 3} enum_Landsat_proc_level;
typedef enum { enum_unknown_sensor = 0, enum_L1MSS = 1, enum_L2MSS = 2, enum_L3MSS = 3, enum_L4MSS = 4, enum_L5MSS = 5, enum_L8 = 8, enum_L9 = 9 } enum_Landsat_sensor;
typedef struct {
	char cImageFilePathName[STRLEN]; // binary 16-bit short
	char cImageFileName[STRLEN];
	char cImageDir[STRLEN];

	int path;
	int row;
	int resolution;
	int width;
	int height;
	int DOY;
	int year;
	int month;
	int day;
	int ID;
	int ID_date;
	enum_Landsat_sensor sensor;
	enum_Landsat_proc_level proc_level;
}Landsat_image;

void SetNewLandsatImage(Landsat_image* pLandsatImage);
void InitializeLandsatImage(Landsat_image* pLandsatImage, char * cLandsatImagePathName);
bool IsLandsatImageValid(Landsat_image* pLandsatImage);
enum_Landsat_sensor GetLandsatSensorFromFileName(char* cFileName);
enum_Landsat_proc_level GetLandsatProcLevelFromFileName(char* cFileName);
int GetLandsatResolution(enum_Landsat_sensor sensor, int band);

bool IsGoodMatching_MSS(bool bSamePathRow, double RMSE, double dMeanShift_x, double dMeanShift_y);
int GetShifts(char* pacDenseMatchingFileName, int iColMax, int iRowMax, int* pdCols, int* pdRows, int iSamplingInterval, float* ptr_fShifts_x, float* ptr_fShifts_y);
int GetShifts_v1(char* pacDenseMatchingFileName, int iColMax, int iRowMax, int* pdCols, int* pdRows, int iSamplingInterval, float fSAMThreshold, float* ptr_fShifts_x, float* ptr_fShifts_y);

class Landsat_registration
{
public:
	Landsat_registration();
	Landsat_registration(char* cFileList_L1TP, char* cFileList_L1GS, char* cOutputDir, int width, int height, int iDenseMatchingInterval);
	~Landsat_registration();

	int L1TPScanning_v1(bool bOutputMatchingSummaryFile = true, int iSensorCombine = 1, int iHalfSearchWindow = 1);

	int InputImagesList(char* cFilesList, enum_Landsat_proc_level proc_level);
	
	bool OutputRegisteredImagesStack(bool bOutputOriginal);

	double Adjustment_RBFs_v1(int RBFs_K, int iControl_idx = 0);
	bool OutputRegisteredImagesStack_RBFs(bool bOutputOriginal, int RBFs_K, int iControl_idx = 0);

	void GetImageSize(int* ptr_width, int* ptr_height) { *ptr_width = m_width; *ptr_height = m_height; }
	Landsat_image* GetImage(int index);

	int GetL1TPMisregistrationsFromScanningResults();
	int GetImagesNum() { return m_ImagesNum; }

private:
	Landsat_image m_pImages_L1TP[MAX_PROC_LANDSAT_IMAGES];
	int m_ImagesNum_L1TP;
	double m_misregistrations_L1TP[MAX_PROC_LANDSAT_IMAGES];

	Landsat_image m_pImages_L1GS[MAX_PROC_LANDSAT_IMAGES];
	int m_ImagesNum_L1GS;

	Landsat_image* m_pImages[MAX_PROC_LANDSAT_IMAGES];
	int m_ImagesNum;

	int m_width;
	int m_height;

	char m_cOutputDir[STRLEN];
	int m_iDenseMatchingInterval;

	// functions
	int GetImageIndex(int imageID);
	int GetL1TPImageIndex(int imageID);

	int GetConnectionGraphFromMatchingFile(char* pacMatchingFileName, bool* ptr_bConnectionMap, bool bIndexFrom_0 = true);
	int FindValid_RBF_Centers_PerImage_v1(char* pacMatchingSummaryFileName, int iRBFs_K, unsigned int* piRBFs_ntie_Map, unsigned int iRBF_ties_num_threshold, double *xs_k_valid_nImages, double *ys_k_valid_nImages, int *K_valid_nImages);

	bool CalculateAdjustmentResiduals_v2(double* pCoefs_all, int* piImages_Control_Flag, int iTransformationType, int* piK_valid_nImages, int RBFs_K, double* ptr_dResiduals, double& ptr_dRMSE, double* ptr_dRMSEsPerImage, double* ptr_dMisregistrationRes_stats, double* ptr_dPre_MisregistrationPerImage, double* ptr_dPost_MisregistrationPerImage, bool bIndexedFrom_0);
	bool InputAdjustedCoefficients(double* ptr_coefs);
	bool InputAdjustedCoefficients_RBFs(double* ptr_coefs, int RBFs_K, int iControl_idx = 0);
};

#endif

