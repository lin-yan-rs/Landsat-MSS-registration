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
	int iDenseMatching_interval = 50; // dense matching every 50 pixels

	iControl_idx = -1; // -1 means all Landsat 8/9 images are used as controls (references)

	// Settings
	// image names in list_all.txt should be Landsat Collection 2 scene IDs that are used to determine Landsat sensors (1-9) and acquisition dates
	// images should be sorted ascendingly by acquisition dates
	// images should have the same size
	sprintf_s(cTarget_pathrow, STRLEN, "042036"); // does not have to a path/row name
	width = 4933; // image size
	height = 4109;
	sprintf_s(cOutputDir, STRLEN, "D:/Lin/MSS/MSS_L1TP_1-5/LSRegistration_MSS1-3_dense_matching/Public_Release/testdata");
	sprintf_s(cImageList_L1TP, STRLEN, "%s/%s_list_all.txt", cOutputDir, cTarget_pathrow);
	
	Landsat_registration* pLandsat_reg = new Landsat_registration(cImageList_L1TP, "", cOutputDir, width, height, iDenseMatching_interval);
	
	// Step 1: many-to-many dense least-squares matching
	// parameters (true, 1, 5)
	//  true: output a matching summary file used in least-squares adjustment
	//  1: matching strategy mode; = 1 means apply many-to-many matching for all image pairs
	//  5: half search window size; = 5 means the search window is 11 x 11
	pLandsat_reg->L1TPScanning_v1(true, 1, 5);
	
	// Step 2: least-squares adjustment to calcualte RBF tranformation parameters for target (non-control) images 
	// parameters (iRBFs_K, iControl_idx)
	//  iRBFs_K: number of K for RBF (see paper); = 64 means image space is evenly split to 8 x 8 grid cells
	//  iControl_idx: index of the control (reference) image, starting from 0; = -1 means all Landsat 8/9 images are used as controls
	pLandsat_reg->Adjustment_RBFs_v1(iRBFs_K, iControl_idx);
	
	// Step 3: output registered images
	// parameters (true, iRBFs_K, iControl_idx)
	//  true: whether output a stack of the original images
	//  iRBFs_K, iControl_idx: same as above
	pLandsat_reg->OutputRegisteredImagesStack_RBFs(true, iRBFs_K, iControl_idx);
	delete pLandsat_reg;
	
	printf(" \nDONE.\n");
	scanf_s(" %d", &iRBFs_K);
	return 0;
}

#endif
