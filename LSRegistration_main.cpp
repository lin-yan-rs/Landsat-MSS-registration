#ifndef LSREGISTRATION_MAIN
#define LSREGISTRATION_MAIN

#include <iostream>
#include <fstream>
#include <stdint.h>
#include <string>
#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#include "LSDepthFirst.h"
#include "LSDepthFirst.h"
#include "Landsat_registration.h"

int main(int argc, char** argv) 
{
	char cImageList_L1TP[STRLEN];
	int width, height;
	char cOutputDir[STRLEN];
	char cTarget_pathrow[STRLEN];
	int iRBFs_K = 64;
	int iControl_idx = 0;
	int iDenseMatching_interval = DENSE_MATCHING_SAMPLING_INTERVAL_MSS;

	iControl_idx = -1;

	sprintf_s(cTarget_pathrow, STRLEN, "042036"); // CA
	width = 4933;
	height = 4109;

	// registration
	sprintf_s(cOutputDir, STRLEN, "D:/Lin/MSS/MSS_L1TP_1-5/LSRegistration_MSS1-3_dense_matching/Public_Release/testdata");
	sprintf_s(cImageList_L1TP, STRLEN, "%s/%s_list_all.txt", cOutputDir, cTarget_pathrow);
	Landsat_registration* pLandsat_reg = new Landsat_registration(cImageList_L1TP, "", cOutputDir, width, height, iDenseMatching_interval);
	pLandsat_reg->L1TPScanning_v1(true, 1, 5);
	pLandsat_reg->Adjustment_RBFs_v1(iRBFs_K, iControl_idx);
	pLandsat_reg->OutputRegisteredImagesStack_RBFs(true, iRBFs_K, iControl_idx);
	delete pLandsat_reg;
	
	printf(" \nDONE.\n");
	scanf_s(" %d", &iRBFs_K);
	return 0;
}

#endif
