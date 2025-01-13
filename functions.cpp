#include "functions.h"

/***
Least squares matching (LSM) using SAM

Parameters
imagel:					left image, which is base
imager:					right image to be resampled to match to left image
(img_rowl, img_coll):	matching window size
(x1, y1):				point position on left image
(x2, y2):				(pass-by-address) initial position on right image corresponding to (x1, y1), and also return as matched position
cr:						(pass-by-address) similarity of image pacthes after final resampling; < 0 if image patches exceed image boundaries, or mean value difference of two image patches is greater than shtMeanDiffThreshold
shtMeanDiffThreshold:	threshold; do not do matching if mean value difference of two image patches is greater than this threshold; default = 9999, meaning the threshold is not used

Notes
1. The resampled image patch on the right image has a 1-pixel buffer compared with the patch on the left image. It is (img_rowl + 2) x (img_coll + 2) pixels rather than img_rowl x img_coll pixels.
The buffer is needed for image gradients calculations. This is the trick in this LSM implementation. Empirically, it improves the matching accuracy
by an order of magnitude compared with the implementations without this buffer trick.
2. If the returend (x2, y2) are the same as the input (x2, y2), and cr > 0, it means the resampled image patch in the first least-squares-fitting iteration has decreased similarity,
and so the fitting is stopped returning the original (x2, y2) as the matched position. In this case, increasing the matching window size can possibly provide successful matching.
See the usage of this function in depth-first registration.
3. SAM is used as the similarity metric rather than correlation. Advantage: insensitive to thresholding (0.995 threshold is safe to use). Disadvantage: yields high (close to 1) similarity value when the matched images contain little variance, like water surface;
this can be handled to by checking the variance levels; in dense matching part of this software, water pixels are excluded simply using a threshold applied to NIR band; see DenseMatching().

7/22/2020: add parameter ptr_paras
7/27/2020: add iMinIter condition
7/29/2020
-add parameter h1_threshold (default = 0.5)
-add check with initial SAM, return -6

_UC (11/3/2020): for GCP chips matching
11/6/2020: fix a bug: if cr_current = -5, cr_current will be set as -6

_precise (12/1/2020): use float type for resampled image
4/8/2021: return real correlation coefficient; make it consistent with latest LSMatching_SAM_UC()
4/9/2021: keep iterations if either SAM or corr increases
4/10/2021: return correlation coefficient before normalization
4/13/2021: fixed bug to return correct h1 and h0
4/19/2021: return iter
4/20/2021: added parameter iMaxIter (default = 20)
7/27/2022: 
-change N[8][8] to N[64]
-change adjustment to float value
***/
int LSMatching_SAM_UC_precise_fast(unsigned char* imagel, int coll, int rowl, unsigned char* imager,
	int colr, int rowr, int img_coll, int img_rowl,
	float x1, float y1, float* ptr_x2, float* ptr_y2, float* ptr_cr, float* ptr_cr_real, float* ptr_cr_real_org, unsigned char ucMeanDiffThreshold, float* ptr_paras, int iMinIter, int iMaxIter, double h1_threshold)
{
	int i, j;
	int x, y, x10, y10, x20, y20;
	float A[8], N[64], l, U[8];
	float gx, gy;
	float h0, h1;					// two radiometric parameters
	double par_x[3], par_y[3];		// six geometric parameters
	float par[8];					// total above eight parameters to be fitted by least squares
	float cr_current, cr_new;		// similarity metric, can be SAM or correlation; SAM is used here
	int	img_rowr, img_colr;			// size of resampled image patch on right image; note it is (img_rowl + 2) x (img_coll + 2) pixels rather than img_rowl x img_coll pixels
	float* imgl0, * imgr0, * imgr1;
	int k, r;
	float* imgl = NULL, * imgr = NULL;
	/// int iMaxIter = 100; // added for MSS; 4/8/2021 changed from 10 to 20
	int iter;
	float cr_initial;

	float cr_real_new, cr_real_current, cr_real_init;

	// added 4/10/2021
	float cr_real_org_new, cr_real_org_current, cr_real_org_init;
	float* imgl_org = NULL;
	float cr_org_new, cr_org_current, cr_org_init; // added 5/1/2021, original SAM with no normalization

	*ptr_cr_real = 0;
	*ptr_cr = 0;

	// added 4/8/2021
	bool bNormalize = true;

	iter = 0;

	// for debugging
	bool bTest = false;
	FILE* ftmp = NULL;
	char* filename = NULL, * cBandNames = NULL, * cBandInfo = NULL;
	int iBandsNum = 0;
	if (bTest)
	{
		filename = (char*)calloc(STRLEN, sizeof(char));
		cBandNames = (char*)calloc(STRLEN_LONG, sizeof(char));
		cBandInfo = (char*)calloc(STRLEN, sizeof(char));
		sprintf_s(filename, STRLEN, "D:/Lin/MSS/MSS_L1TP_1-5/LSM_%d_%d", int(x1 + 0.00001) + 1, int(y1 + 0.000001) + 1);
		ftmp = WriteBinary(filename);
	}

	if (ptr_paras != NULL)
		memset(ptr_paras, 0, 8 * sizeof(float));

	img_rowr = img_rowl + 2;
	img_colr = img_coll + 2;

	imgl = (float*)malloc(img_rowl * img_coll * sizeof(float));
	imgr = (float*)malloc(img_rowr * img_colr * sizeof(float));
	imgl_org = (float*)malloc(img_rowl * img_coll * sizeof(float));

	// get left image patch
	x10 = (int)(x1 + 1e-10);
	par_x[0] = (double)(x10);
	par_x[1] = 1;
	par_x[2] = 0;
	y10 = (int)(y1 + 0.000001);
	par_y[0] = (double)(y10);
	par_y[1] = 0;
	par_y[2] = 1;
	if (!resampling_UC_precise(imagel, coll, rowl, imgl, img_coll, img_rowl, x10, y10, par_x, par_y))
	{
		*ptr_cr = -1;
		if (bTest)
			fclose(ftmp);
		goto freedata;
	}
	memmove(imgl_org, imgl, img_coll * img_rowl * sizeof(float));

	if (bNormalize)
	{
		if (!normalized_FLT_simple(imgl, img_coll * img_rowl, 0, 127))
		{
			*ptr_cr = -7;
			if (bTest)
				fclose(ftmp);
			goto freedata;
		}
	}

	if (bTest)
	{
		sprintf_s(cBandInfo, STRLEN, "0: (%d %d)", x10 + 1, y10 + 1);
		sprintf_s(cBandNames, STRLEN_LONG, "%s, %s", cBandNames, cBandInfo);
		memset(imgr, 0, img_colr * img_rowr * sizeof(unsigned char));
		fwrite(imgr, sizeof(float), img_colr, ftmp); // empty line
		for (j = 0; j < img_rowl; j++)
		{
			fwrite(imgr, sizeof(float), 1, ftmp);
			fwrite(imgl + j * img_coll, sizeof(float), img_coll, ftmp);
			fwrite(imgr, sizeof(float), 1, ftmp);
		}
		fwrite(imgr, sizeof(float), img_colr, ftmp); // empty line
		iBandsNum += 1;
	}

	// get right image patch
	x20 = (int)(*ptr_x2 + 1e-10);
	par_x[0] = (double)(x20);
	par_x[1] = 1;
	par_x[2] = 0;
	y20 = (int)(*ptr_y2 + 1e-10);
	par_y[0] = (double)(y20);
	par_y[1] = 0;
	par_y[2] = 1;
	if (!resampling_UC_precise(imager, colr, rowr, imgr, img_colr, img_rowr, x20, y20, par_x, par_y))
	{
		*ptr_cr = -2;
		if (bTest) // added 4/8/2021
			fclose(ftmp);
		goto freedata;
	}

	// get initial pre-normalization correlation coefficient and SAM
	correlation_coefficient_FLT_v2(imgl_org, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_real_org_current);
	cr_real_org_init = cr_real_org_current;
	SAM_FLT(imgl_org, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_org_current);
	cr_org_init = cr_org_current;

	if (bNormalize)
	{
		if (!normalized_FLT_simple(imgr, img_colr * img_rowr, 0, 127))
		{
			*ptr_cr = -7;
			if (bTest)
				fclose(ftmp);
			goto freedata;
		}
	}

	if (bNormalize == false)
	{
		// check mean difference of two image patches based on shtMeanDiffThreshold
		if (GetMeanDiff_FLT(imgl, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl) >= ucMeanDiffThreshold)
		{
			// mean difference greater than shtMeanDiffThreshold, return
			*ptr_cr = -4;
			if (bTest) // added 4/8/2021
				fclose(ftmp);
			goto freedata;
		}
	}

	// initialize 8 parameters
	par_x[0] = *ptr_x2;
	par_x[1] = 1;
	par_x[2] = 0;
	par_y[0] = *ptr_y2;
	par_y[1] = 0;
	par_y[2] = 1;
	h0 = 0;
	h1 = 1;

	// get initial SAM
	SAM_FLT(imgl, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_current);
	cr_initial = cr_current;

	// get initial correlation coefficient
	correlation_coefficient_FLT_v2(imgl, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_real_current);
	cr_real_init = cr_real_current;

	if (bTest)
	{
		sprintf_s(cBandInfo, STRLEN, "0: (%.2f %.2f) %.4f %.4f", par_x[0] - x20, par_y[0] - y20, cr_current, cr_real_current);
		sprintf_s(cBandNames, STRLEN_LONG, "%s, %s", cBandNames, cBandInfo);
		fwrite(imgr, sizeof(float), img_colr * img_rowr, ftmp);
		iBandsNum += 1;
	}

	// do least squres fitting
	iter = 0;
	while (iter < iMaxIter) // add this condition for MSS
	{
		iter += 1;

		imgl0 = imgl;
		imgr0 = imgr1 = imgr + img_colr + 1;

		memset(U, 0, 8 * sizeof(float));
		memset(N, 0, 64 * sizeof(float));

		// traverse every pixel in img_rowl x img_coll window
		for (i = 0; i < img_rowl; i++)
		{
			y = (int)(-img_rowl / 2.0 + i + 1e-10);
			for (j = 0; j < img_coll; j++)
			{
				x = (int)(-img_coll / 2.0 + j + 1e-10);
				l = *imgl0 - ((*imgr1) * h1 + h0);
				gx = (float)(*(imgr1 + 1) - *(imgr1 - 1)) / 2;
				gy = (float)(*(imgr1 + img_colr) - *(imgr1 - img_colr)) / 2;
				if (ABS(*imgl0) < 1e-6 || ABS(*imgr1) < 1e-6 || ABS(*(imgr1 + 1)) < 1e-6 || ABS(*(imgr1 - 1)) < 1e-6 || ABS(*(imgr1 + img_colr)) < 1e-6 || ABS(*(imgr1 - img_colr)) < 1e-6)
				{
					imgl0++;
					imgr1++;
					continue;
				}

				*A = 1;
				*(A + 1) = *imgr1;
				*(A + 2) = gx;
				*(A + 3) = gx * x;
				*(A + 4) = gx * y;
				*(A + 5) = gy;
				*(A + 6) = gy * x;
				*(A + 7) = gy * y;

				// construct normal equation
				for (k = 0; k < 8; k++)
				{
					for (r = 0; r < 8; r++)
						N[k * 8 + r] += A[k] * A[r];
					U[k] += A[k] * l;
				}

				imgl0++;
				imgr1++;
			}
			imgr0 += img_colr;
			imgr1 = imgr0;
		}

		// solve normal equation Nx = U, x saved in U
		if (!INVSQR1_FLT(N, U, 8))
		{
			if (iter == 1)
				*ptr_cr = -8;
			break;
		}

		// update eight parameters
		memmove(par, U, 8 * sizeof(float));
		h0 = h0 + *par;
		h1 = h1 + *(par + 1);

		*par_x += *(par + 2);
		*(par_x + 1) += *(par + 3);
		*(par_x + 2) += *(par + 4);
		*par_y += *(par + 5);
		*(par_y + 1) += *(par + 6);
		*(par_y + 2) += *(par + 7);

		// resample right image patch using new geometric coefficients
		if (!resampling_UC_precise(imager, colr, rowr, imgr, img_colr, img_rowr, x20, y20, par_x, par_y))
		{
			*ptr_cr = -3;
			if (bTest)
				fclose(ftmp);
			goto freedata;
		}

		// update pre-normorlizationnu correlation coefficient and SAM
		correlation_coefficient_FLT_v2(imgl_org, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_real_org_new);
		SAM_FLT(imgl_org, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_org_new);

		if (bNormalize)
		{
			if (!normalized_FLT_simple(imgr, img_colr * img_rowr, 0, 127))
			{
				*ptr_cr = -7;
				if (bTest)
					fclose(ftmp);
				goto freedata;
			}
		}

		// calculate new similarity (SAM) between left image patch and newly resampled right image patch; conventially correlation_coefficient() is used
		SAM_FLT(imgl, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_new);
		// update correlation coefficient
		correlation_coefficient_FLT_v2(imgl, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_real_new, (float)(h0), (float)(h1));

		if (bTest)
		{
			sprintf_s(cBandInfo, STRLEN, "%d:(%.2f %.2f) %.4f %.4f (%.3f %.3f)  %.2f %.2f", iter, par_x[0] - x20, par_y[0] - y20, cr_new, cr_real_new, par_x[1], par_y[2], h1, h0);
			sprintf_s(cBandNames, STRLEN_LONG, "%s, %s", cBandNames, cBandInfo);
			fwrite(imgr, sizeof(float), img_colr * img_rowr, ftmp);
			iBandsNum += 1;
		}

		// check whether similarity increases
		if (cr_new > cr_current || cr_real_new > cr_real_current || iter <= iMinIter) // 11/3/2017: change from > to >= (the results are consistent with the linux version with the change)
		{
			cr_current = cr_new; // similarity inceasing, go on
			cr_real_current = cr_real_new;
			cr_real_org_current = cr_real_org_new;
			cr_org_current = cr_org_new;
		}
		else
		{
			// similarity deceasing, reverse parameters to values in previous iteration and stop least squares fitting
			*par_x -= *(par + 2);
			*(par_x + 1) -= *(par + 3);
			*(par_x + 2) -= *(par + 4);
			*par_y -= *(par + 5);
			*(par_y + 1) -= *(par + 6);
			*(par_y + 2) -= *(par + 7);

			// added 4/13/2021
			h0 -= *par;
			h1 -= *(par + 1);
			break;
		}
	}

	if (bTest)
	{
		fclose(ftmp);
		OutputEnviHDRFile(filename, "", img_colr, img_rowr, iBandsNum, 4, cBandNames);
	}

	//	if (h1 < h1_threshold || h1 > 1.4)// added 7/22/2020
	if (h1 < h1_threshold)
		*ptr_cr = -5;
	else if (*ptr_cr == 0) // added 11/6/2020
	{
		//	if (cr_current >= cr_initial || cr_real_current >= cr_real_init) // added 7/29/2020; second condition added 5/1/2021
		if (cr_current >= cr_initial || cr_real_current >= cr_real_init) // added 7/29/2020; second condition added 5/1/2021
		{
			// final similarity
			*ptr_cr = MAX(cr_current, cr_org_current); // 5/1/2021, originally equals to cr_current
			*ptr_cr_real = cr_real_current;
			*ptr_cr_real_org = cr_real_org_current;

			if (ptr_paras != NULL)
			{
				ptr_paras[0] = h0;
				ptr_paras[1] = h1;
				ptr_paras[2] = (float)(par_x[0] - *ptr_x2);
				ptr_paras[3] = (float)(par_y[0] - *ptr_y2);
				ptr_paras[4] = (float)par_x[1];
				ptr_paras[5] = (float)par_x[2];
				ptr_paras[6] = (float)par_y[1];
				ptr_paras[7] = (float)par_y[2];
			}

			// final matched position on right image, corresponding to (x1, y1) on left image
			// 11/27/2020, moved below
			*ptr_x2 = (float)(par_x[0]);
			*ptr_y2 = (float)(par_y[0]);
		}
		else
			*ptr_cr = -6;
	}

freedata:
	free(imgl);
	free(imgr);
	free(imgl_org);
	if (bTest)
	{
		free(filename);
		free(cBandNames);
		free(cBandInfo);
	}

	return iter;
}

/***
* Created 9/14/2021
***/
int GetConnectionGraphFromMatchingFile(char* pacMatchingFileName, int n_Images, bool* ptr_bConnectionMap, bool bIndexFrom_0, float fRMSE_threshold, float fShift_threshold, int ntie_threshold)
{
	int iMatchedTilesNum;
	FILE* file_matching = NULL;
	int iMatch;
	int index1, index2;
	int ntie;
	double RMSE, dMeanShift_x, dMeanShift_y;
	char pacMatchingFileDir[STRLEN];
	int iMatch_valid;
	
	int i;

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

		if (RMSE > fRMSE_threshold || sqrt(dMeanShift_x*dMeanShift_x + dMeanShift_y*dMeanShift_y) > fShift_threshold || ntie < ntie_threshold)
			continue;

		ptr_bConnectionMap[index1 * n_Images + index2] = true;
		ptr_bConnectionMap[index2 * n_Images + index1] = true;
		iMatch_valid += 1;
	}

	fclose(file_matching);

	return iMatch_valid;
}

/***
* Created 9/16/2021
* Revised from GetConnectionGraphFromMatchingFile()
* for each dense-matching-pixel, count how many matchings are obtained with other images 
* piImages_Control_Flag: 1 = control, -1 = thrown off image; do not consider in the counting if either image is flagged as -1
***/
int Get_DenseMatching_tie_count_map_PerImage(char* pacMatchingSummaryFileName, int n_Images, int matching_width, int matching_height, int * piImages_Control_Flag, bool bIndexFrom_0, float fRMSE_threshold, float fShift_threshold, int ntie_threshold, float fSAM_threshold, int* pi_tie_CountMapsPerImage)
{
	int iMatchedTilesNum;
	FILE* file_matching = NULL;
	int iMatch;
	int index1, index2;
	int ntie;
	double RMSE, dMeanShift_x, dMeanShift_y;
	char pacMatchingFileName[STRLEN];
	int iMatch_valid;

	int i;
	float* pfMatchingMaps = NULL;
	float* pfSAMMap = NULL;
	int iSize;
	int Row_box, Col_box;

	iSize = matching_width * matching_height;

	if (!(pfMatchingMaps = (float*)calloc(iSize * 4, sizeof(float))))
	{
		printf("\nError in GetMatchingCoverageRatioPerImage(): insufficient memory.\n");
		scanf_s(" %d", &i);
		exit(1);
	}

	memset(pi_tie_CountMapsPerImage, 0, n_Images * iSize * sizeof(int));

	iMatchedTilesNum = Gettxtlines(pacMatchingSummaryFileName);	// number of matched tiles
	file_matching = Readtxt(pacMatchingSummaryFileName);
	iMatch_valid = 0;
	for (iMatch = 0; iMatch < iMatchedTilesNum; iMatch++)
	{
		// note: dMeanShift_x and dMeanShift_y are fitted translation tranformation coefficients a0 and b0	
		fscanf_s(file_matching, "%d%d%lf%lf%lf%d", &index1, &index2, &dMeanShift_x, &dMeanShift_y, &RMSE, &ntie);
		fscanf_s(file_matching, "%s", pacMatchingFileName, STRLEN);
		if (bIndexFrom_0 == false)
		{
			index1 -= 1;
			index2 -= 1;
		}
		if (piImages_Control_Flag[index1] == -1 || piImages_Control_Flag[index2] == -1)
			continue;

		if (RMSE > fRMSE_threshold || sqrt(dMeanShift_x * dMeanShift_x + dMeanShift_y * dMeanShift_y) > fShift_threshold || ntie < ntie_threshold)
			continue;

		pfMatchingMaps = (float*)ReadInputImage(pacMatchingFileName, iSize, 4, 3);
		pfSAMMap = pfMatchingMaps + iSize * 3;
		ntie = LargerOrEqualTo_Num_FLT(pfSAMMap, iSize, fSAM_threshold);
		if (ntie < ntie_threshold)
			continue;

		for (Row_box = 0; Row_box < matching_height; Row_box++)
		{
			for (Col_box = 0; Col_box < matching_width; Col_box++)
			{
				if (pfSAMMap[Row_box * matching_width + Col_box] >= fSAM_threshold)
				{
					pi_tie_CountMapsPerImage[index1 * iSize + Row_box * matching_width + Col_box] += 1;
					pi_tie_CountMapsPerImage[index2 * iSize + Row_box * matching_width + Col_box] += 1;
				}
			}
		}
		iMatch_valid += 1;
	}

	fclose(file_matching);
//	FILE* ftmp = WriteBinary("D:/Lin/MSS/MSS_L1GS_L1-3_WRS1/37_30/L1TP/registration/output_chips_2/360300200/A");
//	fwrite(pi_tie_CountMapsPerImage, sizeof(int), iSize * n_Images, ftmp);
//	fclose(ftmp);
	free(pfMatchingMaps);

	return iMatch_valid;
}

/***
* Created 9/14/2021
***/
int FindUnconnectedImages(char* pacMatchingFileName, int n_Images, int iConnectNumThreshold, bool* ptr_bUnconnectedImagesFlag, bool bIndexFrom_0, float fRMSE_threshold, float fShift_threshold, int ntie_threshold)
{
	bool* pbConnectionMap = NULL;
	int i, k;
	int iUnconnectedImagesNum, iUnconnectedImagesNum_new;
	int iConnectionNum;

	if (!(pbConnectionMap = (bool*)calloc(n_Images * n_Images, sizeof(bool))))
	{
		printf("\nError in FindUnconnectedImages(): insufficient memory.\n");
		scanf_s(" %d", &i);
		exit(1);
	}

	for (i = 0; i < n_Images; i++)
		ptr_bUnconnectedImagesFlag[i] = false;

	iUnconnectedImagesNum = 0;
	while (1)
	{
		iUnconnectedImagesNum_new = 0;

		GetConnectionGraphFromMatchingFile(pacMatchingFileName, n_Images, pbConnectionMap, bIndexFrom_0, fRMSE_threshold, fShift_threshold, ntie_threshold);

		for (i = 0; i < n_Images; i++)
		{
			if (ptr_bUnconnectedImagesFlag[i] == true)
				continue;

			// get number of connections for SAFE i
			iConnectionNum = 0;
			for (k = 0; k < n_Images; k++)
			{
				if (pbConnectionMap[i * n_Images + k] == true)
					iConnectionNum += 1;
			}

			// check whether SAFE is considered connected
			if (iConnectionNum < iConnectNumThreshold)
			{
				ptr_bUnconnectedImagesFlag[i] = true;
				iUnconnectedImagesNum_new += 1;
			}
		}

		iUnconnectedImagesNum += iUnconnectedImagesNum_new;

		if (iUnconnectedImagesNum_new == 0)
			break;
	}

	free(pbConnectionMap);

	return iUnconnectedImagesNum;
}

/***
v1 (9/4/2019):
-add parameter fSAMThreshold;
-get threhsold that is the smaller of fSAMThreshold and SAM stat (max of median and mean), and discard matchings with SAMs smaller than the threshold
v2 (9/14/2021): input dense matching data output by Image_stack_many_to_many_match_in_subset()
***/
int GetShifts_v2(char* pacDenseMatchingFileName, int iColMax, int iRowMax, int* pdCols, int* pdRows, float fSAMThreshold, float* ptr_fShifts_x, float* ptr_fShifts_y)
{
	float* pfParalaxMap_x = NULL, * pfParalaxMap_y = NULL, * pfCorrMap = NULL;
	int n;
	int ntie;
	int Row, Col;
	//	float pfSAM_stats[4];
	//	float median_SAM;
	float SAM_threshold;

	// input densing matching tie points coordinates and coefficients
	pfParalaxMap_x = (float*)ReadInputImage(pacDenseMatchingFileName, iRowMax, iColMax * 4, 3);
	pfParalaxMap_y = pfParalaxMap_x + iRowMax * iColMax;
	pfCorrMap = pfParalaxMap_x + 3 * iRowMax * iColMax;

	/*	if (GetStatistics_FLT(pfCorrMap, iRowMax * iColMax, 0, pfSAM_stats) == false)
		{
			printf("Error in GetShifts_v2(): unable to get SAM statistics.\n");
			scanf_s(" %d", &n);
			exit(1);
		}

		// get SAM threshold
		median_SAM = MAX(pfSAM_stats[2], pfSAM_stats[3]);
		SAM_threshold = MIN(median_SAM, fSAMThreshold);*/
	SAM_threshold = fSAMThreshold;

	// set matched coordinates
	ntie = 0;
	for (Row = 0; Row < iRowMax; Row++)
	{
		for (Col = 0; Col < iColMax; Col++)
		{
			n = Row * iColMax + Col;
			if (pfCorrMap[n] > SAM_threshold)
			{
				pdRows[ntie] = Row;
				pdCols[ntie] = Col;
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
pfAffineCoefs[4]
***/
void GetDecomposedAffineComponents(float *pfAffineCoefs, float *pfScale_x, float *pfScale_y, float *pfShear)
{
	float a1, a2, b1, b2;
	float detA;

	*pfScale_x = 1;
	*pfScale_y = 1;
	*pfShear = 1;

	a1 = pfAffineCoefs[0];
	a2 = pfAffineCoefs[1];
	b1 = pfAffineCoefs[2];
	b2 = pfAffineCoefs[3];

	detA = a1 * b2 - a2 * b1;
	if (ABS(detA) < 1e-8)
		return;

	*pfScale_x = sqrt(a1 * a1 + a2 * a2);
	*pfScale_y = detA / (*pfScale_x);
	*pfShear = (a1 * b1 + a2 * b2) / detA;
	
	return;
}

/***
Created 7/28/2021
Use more strict criteria for L1TP images matching
***/
bool Check_LSM_components(float* pfLSM_paras, float fSAM, float corr_real, float corr_real_org, float* pfScale_distort, float* pfShear_distort, bool bL1TP_matching)
{
	// added 7/28/2021
	float x_distort, y_distort, distort;
	float distort_xy;
	float h1;
	bool bSuc;
	// decomposed affine transformation
	float scale_x, scale_y;
	float shear;
	float scale_distort, shear_distort;

	float sam_threshold = 0.96f;

	x_distort = ABS(pfLSM_paras[4] - 1.0f);
	y_distort = ABS(pfLSM_paras[7] - 1.0f);
	distort = MAX(x_distort, y_distort);
	distort_xy = MAX(ABS(pfLSM_paras[5]), ABS(pfLSM_paras[6]));
	h1 = pfLSM_paras[1];

	// get decomposed affine transformation components
	GetDecomposedAffineComponents(pfLSM_paras + 4, &scale_x, &scale_y, &shear);
	scale_distort = MAX(ABS(scale_x - 1), ABS(scale_y - 1));
	shear_distort = ABS(shear);
	distort = scale_distort;

	bSuc = false;
	if (fSAM >= sam_threshold && corr_real >= 0.5 && distort <= 0.04 && pfLSM_paras[1] > 0.3f && pfLSM_paras[1] < 1.6f) // changed corr_real >= 0.6 to 0.5
	{
		bSuc = true;
		if (distort > 0.02f && (pfLSM_paras[1] >= 1.5f || pfLSM_paras[1] < 0.5f || corr_real_org < 0.5))
			bSuc = false;
		if (distort > 0.03f && (pfLSM_paras[1] >= 1.2f || shear_distort > 0.05 || corr_real_org < 0.6))
			bSuc = false;
		if (shear_distort > 0.055)
			bSuc = false;

		if (bL1TP_matching == false)
		{
			if (corr_real < 0.6) // now it is equivalent to GCP lib matching criteria
				bSuc = false;
		}

		// more strict criteria for L1TP pixel-to-pixel matching
		if (bL1TP_matching == true)
		{
			if (distort > 0.02) // 7/22/2022, changed from 0.02 to 0.025
				bSuc = false;

			if (distort > 0.15f && (pfLSM_paras[1] >= 1.5f || pfLSM_paras[1] < 0.5f || corr_real_org < 0.5 || corr_real < 0.6))
				bSuc = false;
		}
	}

	if (pfScale_distort != NULL)
		*pfScale_distort = distort;
	if (pfShear_distort != NULL)
		*pfShear_distort = shear_distort;

	return bSuc;
}

/***
Calculate mean value difference between two input image patches
***/
float GetMeanDiff_FLT(float* img1, int col1, int row1, float* img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height)
{
	int i, j, n;
	int	height2, width2;
	float* image1, * image10, * image2, * image20;
	float pixel1, pixel2;
	double diff;

	height2 = height / 2;
	width2 = width / 2;

	n = width * height;

	image10 = image1 = img1 + (y1 - height2) * col1 + (x1 - width2);
	image20 = image2 = img2 + (y2 - height2) * col2 + (x2 - width2);
	diff = 0;
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			pixel1 = *image1;
			pixel2 = *image2;
			if (ABS(pixel1) > 1e-6 && ABS(pixel2) > 1e-6)
				diff += (pixel1 > pixel2) ? pixel1 - pixel2 : pixel2 - pixel1;
			image1++;
			image2++;
		}
		image10 += col1;  image1 = image10;
		image20 += col2;  image2 = image20;
	}
	diff /= n;

	return (float)(diff);
}

/***
SAM calculation
***/
void SAM_FLT(float* img1, int col1, int row1, float* img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height, float* ptr_cr)
{
	int i, j, n;
	int	height2, width2;
	double Sxx, Syy, Sxy;
	float* image1, * image10, * image2, * image20;
	float pixel1, pixel2;

	height2 = height / 2;
	width2 = width / 2;

	n = width * height;

	Sxx = 0.0;
	Syy = 0.0;
	Sxy = 0;
	image10 = image1 = img1 + (y1 - height2) * col1 + (x1 - width2);
	image20 = image2 = img2 + (y2 - height2) * col2 + (x2 - width2);

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			pixel1 = *image1;
			pixel2 = *image2;
			if (ABS(pixel1) > 1e-6 && ABS(pixel2) > 1e-6)
			{
				Sxx += pixel1 * pixel1;
				Syy += pixel2 * pixel2;
				Sxy += pixel1 * pixel2;
			}
			image1++;
			image2++;
		}
		image10 += col1;  image1 = image10;
		image20 += col2;  image2 = image20;
	}
	if (Sxx > 0 && Syy > 0)
		*ptr_cr = (float)(Sxy / sqrt(Sxx * Syy));
	else
		*ptr_cr = 0;

	return;
}

int ConvertToDOY(int iYear, int iMonth, int iDay)
{
	int i;
	int iDOY = 0;

	for (i = 0; i < iMonth; i++)
	{
		if (i == 0)
			iDOY += iDay;
		else if (i == 1)
			iDOY += 31;
		else if (i == 2)
			iDOY += (iYear % 4) ? 28 : 29;
		else if (i == 3)
			iDOY += 31;
		else if (i == 4)
			iDOY += 30;
		else if (i == 5)
			iDOY += 31;
		else if (i == 6)
			iDOY += 30;
		else if (i == 7)
			iDOY += 31;
		else if (i == 8)
			iDOY += 31;
		else if (i == 9)
			iDOY += 30;
		else if (i == 10)
			iDOY += 31;
		else if (i == 11)
			iDOY += 30;
	}

	return iDOY;
}

/***
Least squares matching (LSM) using SAM

Parameters
imagel:					left image, which is base
imager:					right image to be resampled to match to left image
(img_rowl, img_coll):	matching window size
(x1, y1):				point position on left image
(x2, y2):				(pass-by-address) initial position on right image corresponding to (x1, y1), and also return as matched position
cr:						(pass-by-address) similarity of image pacthes after final resampling; < 0 if image patches exceed image boundaries, or mean value difference of two image patches is greater than shtMeanDiffThreshold
shtMeanDiffThreshold:	threshold; do not do matching if mean value difference of two image patches is greater than this threshold; default = 9999, meaning the threshold is not used

Notes
1. The resampled image patch on the right image has a 1-pixel buffer compared with the patch on the left image. It is (img_rowl + 2) x (img_coll + 2) pixels rather than img_rowl x img_coll pixels.
The buffer is needed for image gradients calculations. This is the trick in this LSM implementation. Empirically, it improves the matching accuracy
by an order of magnitude compared with the implementations without this buffer trick.
2. If the returend (x2, y2) are the same as the input (x2, y2), and cr > 0, it means the resampled image patch in the first least-squares-fitting iteration has decreased similarity,
and so the fitting is stopped returning the original (x2, y2) as the matched position. In this case, increasing the matching window size can possibly provide successful matching.
See the usage of this function in depth-first registration.
3. SAM is used as the similarity metric rather than correlation. Advantage: insensitive to thresholding (0.995 threshold is safe to use). Disadvantage: yields high (close to 1) similarity value when the matched images contain little variance, like water surface; 
this can be handled to by checking the variance levels; in dense matching part of this software, water pixels are excluded simply using a threshold applied to NIR band; see DenseMatching(). 

7/22/2020: add parameter ptr_paras
7/27/2020: add iMinIter condition
7/29/2020
-add parameter h1_threshold (default = 0.5)
-add check with initial SAM, return -6
11/6/2020: fixed some bugs
***/
void LSMatching_SAM(short int *imagel, int coll, int rowl, short int *imager,
	int colr, int rowr, int img_coll, int img_rowl,
	float x1, float y1, float *ptr_x2, float *ptr_y2, float *ptr_cr, short int shtMeanDiffThreshold, float* ptr_paras, int iMinIter, double h1_threshold)
{
	int i, j;
	int x, y, x10, y10, x20, y20;
	double A[8], N[8][8], l, U[8];
	double gx, gy;
	double h0, h1;					// two radiometric parameters
	double par_x[3], par_y[3];		// six geometric parameters
	double par[8];					// total above eight parameters to be fitted by least squares
	float cr_current, cr_new;		// similarity metric, can be SAM or correlation; SAM is used here
	int	img_rowr, img_colr;			// size of resampled image patch on right image; note it is (img_rowl + 2) x (img_coll + 2) pixels rather than img_rowl x img_coll pixels
	short int *imgl0, *imgr0, *imgr1;
	int k, r;
	short int *imgl = NULL, *imgr = NULL;
	int iMaxIter = 10; // added for MSS
	int iter;
	float cr_initial;

	// for debugging
	bool bTest = false;
	FILE* ftmp = NULL;
	char filename[STRLEN] = "", cBandNames[STRLEN_LONG] = "", cBandInfo[STRLEN] = "";
	int iBandsNum = 0;
	if (bTest)
	{
		sprintf_s(filename, STRLEN, "D:/Lin/MSS/LSM_%d_%d", int(x1 + 0.00001) + 1, int(y1 + 0.000001) + 1);
		ftmp = WriteBinary(filename);
	}

	if (ptr_paras != NULL)
		memset(ptr_paras, 0, 8 * sizeof(float));

	img_rowr = img_rowl + 2;
	img_colr = img_coll + 2;

	imgl = (short int *)malloc(img_rowl*img_coll *sizeof(short int));
	imgr = (short int *)malloc(img_rowr*img_colr *sizeof(short int));

	// get left image patch
	x10 = (int)(x1 + 1e-10);
	par_x[0] = (double)(x10);
	par_x[1] = 1;
	par_x[2] = 0;
	y10 = (int)(y1 + 0.000001);
	par_y[0] = (double)(y10);
	par_y[1] = 0;
	par_y[2] = 1;
	if (!resampling(imagel, coll, rowl, imgl, img_coll, img_rowl, x10, y10, par_x, par_y))
	{
		*ptr_cr = -1;
		free(imgl);
		free(imgr);
		return;
	}

	if (bTest)
	{
		sprintf_s(cBandInfo, STRLEN, "0: (%d %d)", x10 + 1, y10 + 1);
		sprintf_s(cBandNames, STRLEN_LONG, "%s, %s", cBandNames, cBandInfo);
		memset(imgr, 0, img_colr * img_rowr * sizeof(short));
		fwrite(imgr, sizeof(short), img_colr, ftmp); // empty line
		for (j = 0; j < img_rowl; j++)
		{
			fwrite(imgr, sizeof(short int), 1, ftmp);
			fwrite(imgl + j * img_coll, sizeof(short int), img_coll, ftmp);
			fwrite(imgr, sizeof(short int), 1, ftmp);
		}
		fwrite(imgr, sizeof(short), img_colr, ftmp); // empty line
		iBandsNum += 1;
	}

	// get right image patch
	x20 = (int)(*ptr_x2 + 1e-10);
	par_x[0] = (double)(x20);
	par_x[1] = 1;
	par_x[2] = 0;
	y20 = (int)(*ptr_y2 + 1e-10);
	par_y[0] = (double)(y20);
	par_y[1] = 0;
	par_y[2] = 1;
	if (!resampling(imager, colr, rowr, imgr, img_colr, img_rowr, x20, y20, par_x, par_y))
	{
		*ptr_cr = -2;
		free(imgl);
		free(imgr);
		return;
	}

	// check mean difference of two image patches based on shtMeanDiffThreshold
	if (GetMeanDiff(imgl, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl) >= shtMeanDiffThreshold)
	{
		// mean difference greater than shtMeanDiffThreshold, return
		*ptr_cr = -4;
		free(imgl);
		free(imgr);
		return;
	}

	// initialize 8 parameters
	par_x[0] = *ptr_x2;
	par_x[1] = 1;
	par_x[2] = 0;
	par_y[0] = *ptr_y2;
	par_y[1] = 0;
	par_y[2] = 1;
	h0 = 0;
	h1 = 1;

	// get initial SAM
	SAM(imgl, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_current);
	cr_initial = cr_current;

	if (bTest)
	{
		sprintf_s(cBandInfo, STRLEN, "0: (%.2f %.2f) %.6f", par_x[0] - x20, par_y[0] - y20, cr_current);
		sprintf_s(cBandNames, STRLEN_LONG, "%s, %s", cBandNames, cBandInfo);
		fwrite(imgr, sizeof(short), img_colr * img_rowr, ftmp);
		iBandsNum += 1;
	}

	// do least squres fitting
	iter = 0;
	while (iter < iMaxIter) // add this condition for MSS
	{
		iter += 1;

		imgl0 = imgl;
		imgr0 = imgr1 = imgr + img_colr + 1;

		memset(U, 0, 8 * sizeof(double));
		memset(N, 0, 64 * sizeof(double));

		// traverse every pixel in img_rowl x img_coll window
		for (i = 0; i<img_rowl; i++)
		{
			y = (int)(-img_rowl / 2.0 + i + 1e-10);
			for (j = 0; j<img_coll; j++)
			{
				x = (int)(-img_coll / 2.0 + j + 1e-10);
				l = *imgl0 - ((*imgr1)*  h1  + h0);
				gx = (double)(*(imgr1 + 1) - *(imgr1 - 1)) / 2;
				gy = (double)(*(imgr1 + img_colr) - *(imgr1 - img_colr)) / 2;

				*A = 1;
				*(A + 1) = *imgr1;
				*(A + 2) = gx;
				*(A + 3) = gx*x;
				*(A + 4) = gx*y;
				*(A + 5) = gy;
				*(A + 6) = gy*x;
				*(A + 7) = gy*y;

				// construct normal equation
				for (k = 0; k < 8; k++)
				{
					for (r = 0; r < 8; r++)
						N[k][r] += A[k] * A[r];
					U[k] += A[k] * l;
				}

				imgl0++;
				imgr1++;
				x++;
			}
			imgr0 += img_colr;
			imgr1 = imgr0;
			y++;
		}

		// solve normal equation Nx = U, x saved in U
		if (!INVSQR1((double*)N, U, 8))
			break;

		// update eight parameters
		memmove(par, U, 8 * sizeof(double));
		h0 = h0 + *par;
		h1 = h1 + *(par + 1);

		*par_x += *(par + 2);
		*(par_x + 1) += *(par + 3);
		*(par_x + 2) += *(par + 4);
		*par_y += *(par + 5);
		*(par_y + 1) += *(par + 6);
		*(par_y + 2) += *(par + 7);

		// resample right image patch using new geometric coefficients
		if (!resampling(imager, colr, rowr, imgr, img_colr, img_rowr, x20, y20, par_x, par_y))
		{
			*ptr_cr = -3;
			free(imgl);
			free(imgr);
			return;
		}

		// calculate new similarity (SAM) between left image patch and newly resampled right image patch; conventially correlation_coefficient() is used
		SAM(imgl, img_coll, img_rowl, imgr, img_colr, img_rowr, img_coll / 2, img_rowl / 2, img_colr / 2, img_rowr / 2, img_coll, img_rowl, &cr_new);
	
		if (bTest)
		{
			sprintf_s(cBandInfo, STRLEN, "%d:(%.2f %.2f)  %.2f-%.2f  %.6f", iter, par_x[0] - x20, par_y[0] - y20, h0, h1, cr_new);
			sprintf_s(cBandNames, STRLEN_LONG, "%s, %s", cBandNames, cBandInfo);
			fwrite(imgr, sizeof(short), img_colr * img_rowr, ftmp);
			iBandsNum += 1;
		}

		// check whether similarity increases
		if (cr_new >= cr_current || iter <= iMinIter) // 11/3/2017: change from > to >= (the results are consistent with the linux version with the change)
			cr_current = cr_new; // similarity inceasing, go on
		else
		{
			// similarity deceasing, reverse parameters to values in previous iteration and stop least squares fitting
			*par_x -= *(par + 2);
			*(par_x + 1) -= *(par + 3);
			*(par_x + 2) -= *(par + 4);
			*par_y -= *(par + 5);
			*(par_y + 1) -= *(par + 6);
			*(par_y + 2) -= *(par + 7);

			break;
		}
	}

	if (bTest)
	{
		fclose(ftmp);
		OutputEnviHDRFile(filename, "", img_colr, img_rowr, iBandsNum, 2, cBandNames);
	}

	if (h1 < h1_threshold || h1 > 1.4)// added 7/22/2020
		*ptr_cr = -5; // 11/6/2020, fixed cr_current = -5 to *ptr_cr = -5
	else // added 11/6/2020
	{
		if (cr_current >= cr_initial) // added 7/29/2020
		{
			// final similarity
			*ptr_cr = cr_current;

			// final matched position on right image, corresponding to (x1, y1) on left image
			*ptr_x2 = (float)(par_x[0]);
			*ptr_y2 = (float)(par_y[0]);

			if (ptr_paras != NULL)
			{
				ptr_paras[0] = (float)h0;
				ptr_paras[1] = (float)h1;
				ptr_paras[2] = (float)(par_x[0] - *ptr_x2);
				ptr_paras[3] = (float)(par_y[0] - *ptr_y2);
				ptr_paras[4] = (float)par_x[1];
				ptr_paras[5] = (float)par_x[2];
				ptr_paras[6] = (float)par_y[1];
				ptr_paras[7] = (float)par_y[2];
			}
		}
		else
			*ptr_cr = -6;
	}

	free(imgl);
	free(imgr);

	return;
}

/***
Bilinear resampling using geometric affine transformation parameteres
11/6/2020
-change
if (x20 < 0 || x20 >= col - 1 || y20 < 0 || y20 >= row - 1)
TO
if (x20 < 0 || x20 > col - 1 || y20 < 0 || y20 > row - 1)
***/
int resampling(short int *image, int col, int row, short int *img, int img_col, int img_row, int x1, int y1, double *par_x, double *par_y)
{
	int i, j, i1, i2, j1, j2;
	short int g, *image0, *img0;
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
			x2 = par_x[0] + j*par_x[1] + i*par_x[2];		//	x2=x+a0+x*a1+y*a2
			y2 = par_y[0] + j*par_y[1] + i*par_y[2];		//	y2=y+b0+y*b1+x*b2
			x20 = (int)(x2 + 1e-10);
			y20 = (int)(y2 + 1e-10);

			if (x20 < 0 || x20 > col - 1 || y20 < 0 || y20 > row - 1)
				return 0;

			image0 = image + y20*col + x20;
			p = x2 - x20;
			q = y2 - y20;
			g = 0;
			bit = 0;
			if ((1 - p) * (1 - q) > 1e-10) // top left (x20, y20)
				bit += (*image0)*(1 - p) * (1 - q);
			if (p * (1 - q) > 1e-10 && x20 + 1 <= col - 1) // top right (x20+1, y20)
				bit += (*(image0 + 1)) * p * (1 - q);
			if ((1 - p) * q > 1e-10 && y20 + 1 <= row - 1) // bottom left (x20, y20+1)
				bit += (*(image0 + col)) * (1 - p) * q;
			if (p * q > 1e-10 && x20 + 1 <= col - 1 && y20 + 1 <= row - 1) // bottom right (x20+1, y20+1)
				bit += (*(image0 + col + 1)) * p * q;
			g = (short int)(bit + 0.5f);
			*img0++ = g;
		}
	}

	return 1;
}

/***
SAM calculation
***/
void SAM(short int *img1, int col1, int row1, short int *img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height, float *ptr_cr)
{
	int i, j, n;
	int	height2, width2;
	double Sxx, Syy, Sxy; // 7/27/2021, changed from float to double, same as SAM_UC() and SAM_FLT()
	short int *image1, *image10, *image2, *image20;
	short pixel1, pixel2;

	height2 = height / 2;
	width2 = width / 2;

	n = width*height;

	Sxx = 0.0;
	Syy = 0.0;
	Sxy = 0;
	image10 = image1 = img1 + (y1 - height2)*col1 + (x1 - width2);
	image20 = image2 = img2 + (y2 - height2)*col2 + (x2 - width2);

	for (i = 0; i<height; i++)
	{
		for (j = 0; j < width; j++)
		{
			pixel1 = *image1;
			pixel2 = *image2;
			Sxx += pixel1*pixel1;
			Syy += pixel2*pixel2;
			Sxy += pixel1*pixel2;
			image1++;
			image2++;
		}
		image10 += col1;  image1 = image10;
		image20 += col2;  image2 = image20;
	}
	if (Sxx > 0 && Syy > 0)
		*ptr_cr = (float)(Sxy / sqrt(Sxx*Syy));
	else
		*ptr_cr = 0;

	return;
}

/***
Calculate mean value difference between two input image patches
***/
short int GetMeanDiff(short int *img1, int col1, int row1, short int *img2, int col2, int row2, int x1, int y1, int x2, int y2, int width, int height)
{
	int i, j, n;
	int	height2, width2;
	short int *image1, *image10, *image2, *image20;
	short pixel1, pixel2;
	float diff;

	height2 = height / 2;
	width2 = width / 2;

	n = width*height;

	image10 = image1 = img1 + (y1 - height2)*col1 + (x1 - width2);
	image20 = image2 = img2 + (y2 - height2)*col2 + (x2 - width2);
	diff = 0;
	for (i = 0; i<height; i++)
	{
		for (j = 0; j < width; j++)
		{
			pixel1 = *image1;
			pixel2 = *image2;
			diff += (pixel1 > pixel2) ? pixel1 - pixel2 : pixel2 - pixel1;
			image1++;
			image2++;
		}
		image10 += col1;  image1 = image10;
		image20 += col2;  image2 = image20;
	}
	diff /= n;

	return (short int)(diff + 0.5f);
}

/***
Solve equation Ax = B
x saved in B
***/
int INVSQR1(double *A, double *B, int n)
{
	int k, i, j, i0 = 0;
	double C;
	double T;

	for (k = 0; k<n; k++)
	{
		C = 0;
		for (i = k; i<n; i++)
		{
			if (fabs(A[i*n + k]) >= fabs(C))
			{
				C = A[i*n + k];
				i0 = i;
			}
		}
		if (i != k)
		{
			for (j = k; j<n; j++)
			{
				T = A[k*n + j];
				A[k*n + j] = A[i0*n + j];
				A[i0*n + j] = T;
			}
			T = B[k];
			B[k] = B[i0];
			B[i0] = T;
		}
		if (fabs(C) <= 0.0)
			return 0;
		C = 1 / C;
		for (j = k + 1; j<n; j++)
		{
			A[k*n + j] *= C;
			for (i = k + 1; i<n; i++)
			{
				A[i*n + j] = A[i*n + j] - A[i*n + k] * A[k*n + j];
			}
		}
		B[k] *= C;
		for (i = k + 1; i<n; i++)
		{
			B[i] = B[i] - B[k] * A[i*n + k];
		}
	}
	for (i = n - 2; i >= 0; i--)
	{
		for (j = i + 1; j<n; j++)
		{
			B[i] = B[i] - A[i*n + j] * B[j];
		}
	}

	return 1;
}

int INVSQR1_FLT(float* A, float* B, int n)
{
	int k, i, j, i0 = 0;
	float C;
	float T;

	for (k = 0; k < n; k++)
	{
		C = 0;
		for (i = k; i < n; i++)
		{
			if (fabs(A[i * n + k]) >= fabs(C))
			{
				C = A[i * n + k];
				i0 = i;
			}
		}
		if (i != k)
		{
			for (j = k; j < n; j++)
			{
				T = A[k * n + j];
				A[k * n + j] = A[i0 * n + j];
				A[i0 * n + j] = T;
			}
			T = B[k];
			B[k] = B[i0];
			B[i0] = T;
		}
		if (fabs(C) <= 0.0)
			return 0;
		C = 1 / C;
		for (j = k + 1; j < n; j++)
		{
			A[k * n + j] *= C;
			for (i = k + 1; i < n; i++)
			{
				A[i * n + j] = A[i * n + j] - A[i * n + k] * A[k * n + j];
			}
		}
		B[k] *= C;
		for (i = k + 1; i < n; i++)
		{
			B[i] = B[i] - B[k] * A[i * n + k];
		}
	}
	for (i = n - 2; i >= 0; i--)
	{
		for (j = i + 1; j < n; j++)
		{
			B[i] = B[i] - A[i * n + j] * B[j];
		}
	}

	return 1;
}

/***
Solve equation Ax = B
x saved in B
v2 (12/13/2017): flag x whose corresponding B is 0; solve the remaining x
9/3/2019: changed condition from ABS(B[i]) < 1e-10 to ABS(B[i]) < 1e-20 to flag unsolved coefficients
***/
int INVSQR2(double* A, double* B, int n, bool* pbSolvedUnknownsFlags)
{
	double* A_ = NULL, * B_ = NULL;
	int n_;
	int i, k;
	int idx, idx_x, idx_y;

	for (i = 0; i < n; i++)
		pbSolvedUnknownsFlags[i] = true;

	// get number of valid unknowns and flag invalid unknowns in pbSolvedUnknownsFlags
	n_ = n;
	for (i = 0; i < n; i++)
	{
		if (ABS(B[i]) < 1e-20)
		{
			pbSolvedUnknownsFlags[i] = false;
			n_ -= 1;
		}
	}

	if (!(A_ = (double*)calloc(n_ * n_, sizeof(double)))
		|| !(B_ = (double*)calloc(n_, sizeof(double))))
	{
		printf("\nError in INVSQR2(): insufficient memory.\n");
		scanf_s(" %d", &n_);
		exit(1);
	}

	// put valid values into A_ and B_
	idx = 0;
	for (i = 0; i < n; i++)
	{
		if (pbSolvedUnknownsFlags[i] == true)
		{
			B_[idx] = B[i];
			idx += 1;
		}
	}

	idx_y = 0;
	for (i = 0; i < n; i++)
	{
		idx_x = 0;
		if (pbSolvedUnknownsFlags[i] == true)
		{
			for (k = 0; k < n; k++)
			{
				if (pbSolvedUnknownsFlags[k] == true)
				{
					A_[idx_y * n_ + idx_x] = A[i * n + k];
					idx_x += 1;
				}
			}
			idx_y += 1;
		}
	}

	// solve it
	if (INVSQR1(A_, B_, n_) == 0)
		return 0;

	// put solved unknowns to B
	memset(B, 0, n * sizeof(double));
	idx = 0;
	for (i = 0; i < n; i++)
	{
		if (pbSolvedUnknownsFlags[i] == true)
		{
			B[i] = B_[idx];
			idx += 1;
		}
	}

	return 1;
}

/***
Calculate mean and standard deviation
***/
void Std1(double *x, int n, double *ptr_stdv, double *ptr_ave)
{
	int i; 

	*ptr_ave = 0;
	*ptr_stdv = 0;
	*ptr_ave = Mean1(x, n);
	for (i = 0; i < n; i++)
		*ptr_stdv += (*(x + i) - (*ptr_ave))*(*(x + i) - (*ptr_ave));

	*ptr_stdv /= (n - 1);
	*ptr_stdv = sqrt(*ptr_stdv);

	return;
}

void Std1_FLT(float* x, int n, float* ptr_stdv, float* ptr_ave)
{
	int i;

	*ptr_ave = 0;
	*ptr_stdv = 0;
	*ptr_ave = Mean1_FLT(x, n);
	for (i = 0; i < n; i++)
		*ptr_stdv += (*(x + i) - (*ptr_ave)) * (*(x + i) - (*ptr_ave));

	*ptr_stdv /= (n - 1);
	*ptr_stdv = sqrt(*ptr_stdv);

	return;
}

/***
Calculate mean and RMSE
***/
void RMSE1(double *x, int n, int t, double *ptr_rmse, double *ptr_avg)
{
	int i;

	if (n <= t)
	{
		*ptr_rmse = 9999;
		*ptr_avg = 9999;
		return;
	}

	*ptr_avg = 0;
	*ptr_rmse = 0;

	*ptr_avg = Mean1(x, n);
	for (i = 0; i < n; i++)
		*ptr_rmse += (*(x + i))*(*(x + i));

	*ptr_rmse /= (n - t);
	*ptr_rmse = sqrt(*ptr_rmse);

	return;
}

double Mean1(double *x, int n)
{
	double sum = 0;
	int i;

	for (i = 0; i<n; i++)
		sum += x[i];

	return sum / n;
}

float Mean1_FLT(float* x, int n)
{
	double sum = 0;
	int i;

	for (i = 0; i < n; i++)
		sum += x[i];

	return (float)(sum / n);
}

unsigned char Mean1_UC(unsigned char* x, int n)
{
	float sum = 0;
	int i;

	for (i = 0; i < n; i++)
		sum += (float)(x[i]);

	return (unsigned char)(sum / n + 0.5f);
}

int LargerOrEqualTo_Num_INT(int *piData, int n, int iValue)
{
	int num, i;
	int *piPtr = piData;

	num = 0;
	for (i = 0; i < n; i++)
		num += (((*piPtr++) >= iValue) ? 1 : 0);

	return num;
}

int LargerOrEqualTo_Num_FLT(float* pData, int n, float value)
{
	int num, i;
	float* pPtr = pData;

	num = 0;
	for (i = 0; i < n; i++)
		num += (((*pPtr++) >= value) ? 1 : 0);

	return num;
}

/***
Fit six affine transformation coefficients in
x2 = a0 + a1*x1 + a2*y1
y2 = b0 + b1*x1 + b2*y1

Returned parameters
Coefs:			[a0, a1, a2, b0, b1, b2]
Errors:			fitting residuals for n points
errors_mean:	mean residual
fitting_rmse:	fitting RMSE
***/
void FitAffineTransform(double *x1, double *y1, double *x2, double *y2, int n, double Coefs[], double *ptr_Errors, double *ptr_errors_mean, double *ptr_fitting_rmse)
{
	double U[6];
	int iIterNum;
	int i, k, r;
	double L;
	double A[6], N[6][6];
	double xdif, ydif;
	double fitting_rmse_, errors_mean_;

	memset(ptr_Errors, 0, n * sizeof(double));
	*ptr_errors_mean = 0;
	*ptr_fitting_rmse = 0;

	if (n < 6)
		return;

	memset(Coefs, 0, 6 * sizeof(double));
	Coefs[1] = 1;
	Coefs[5] = 1;
	iIterNum = 0;
	while (iIterNum < 10)
	{
		memset(A, 0, 6 * sizeof(double)); // added 8/9/2021
		memset(N, 0, 36 * sizeof(double));
		memset(U, 0, 6 * sizeof(double));

		for (i = 0; i < n; i++)
		{
			A[0] = 1;
			A[1] = x1[i];
			A[2] = y1[i];
			A[3] = 0;
			A[4] = 0;
			A[5] = 0;
			L = x2[i] - (Coefs[0] + Coefs[1] * x1[i] + Coefs[2] * y1[i]);
			for (k = 0; k < 6; k++)
			{
				for (r = 0; r < 6; r++)
					N[k][r] += A[k] * A[r];
				U[k] += A[k] * L;
			}
			A[0] = 0;
			A[1] = 0;
			A[2] = 0;
			A[3] = 1;
			A[4] = x1[i];
			A[5] = y1[i];
			L = y2[i] - (Coefs[3] + Coefs[4] * x1[i] + Coefs[5] * y1[i]);
			for (k = 0; k < 6; k++)
			{
				for (r = 0; r < 6; r++)
					N[k][r] += A[k] * A[r];
				U[k] += A[k] * L;
			}
		}

		if (!INVSQR1((double*)N, U, 6))
			break;

		for (k = 0; k < 6; k++)
			Coefs[k] += U[k];

		if (abs(U[0]) < DELTA_LIMIT && abs(U[3]) < DELTA_LIMIT)
			break;

		iIterNum += 1;
	}

	memmove(U, Coefs, 6 * sizeof(double));
	for (i = 0; i < n; i++)
	{
		xdif = U[0] + U[1] * x1[i] + U[2] * y1[i] - x2[i];
		ydif = U[3] + U[4] * x1[i] + U[5] * y1[i] - y2[i];
		ptr_Errors[i] = sqrt(xdif*xdif + ydif*ydif);
	}

	RMSE1(ptr_Errors, n, 6, &fitting_rmse_, &errors_mean_);

	*ptr_errors_mean = errors_mean_;
	*ptr_fitting_rmse = fitting_rmse_;

	return;
}

/***
Fit two translation transformation coefficients in
x2 = a0 + x1
y2 = b0 + y1

a0 + x1 - x2 = v
->
a0^ + a0_init + x1 - x2 = v

a0^ - L = v

=> L = x2 - (a0_init + x1)

Simply: a0 + x1 - x2 = -L, i.e. L = -(a0 + x1 - x2)

Returned parameters
Coefs:			[a0, b0]
Errors:			fitting residuals for n points
errors_mean:	mean residual
fitting_rmse:	fitting RMSE
***/
void FitTranslationTransform(double *x1, double *y1, double *x2, double *y2, int n, double Coefs[], double *ptr_Errors, double *ptr_errors_mean, double *ptr_fitting_rmse)
{
	double U[2];
	int iIterNum;
	int i, k, r;
	double L;
	double A[2], N[2][2];
	double xdif, ydif;
	double fitting_rmse_, errors_mean_;

	memset(ptr_Errors, 0, n * sizeof(double));
	*ptr_errors_mean = 0;
	*ptr_fitting_rmse = 0;

	if (n < 2)
		return;

	memset(Coefs, 0, 2 * sizeof(double));
	iIterNum = 0;
	while (iIterNum < 10)
	{
		memset(A, 0, 2 * sizeof(double)); // added 8/9/2021
		memset(N, 0, 4 * sizeof(double));
		memset(U, 0, 2 * sizeof(double));

		for (i = 0; i < n; i++)
		{
			A[0] = 1;
			A[1] = 0;
			L = x2[i] - (Coefs[0] + x1[i]);
			for (k = 0; k < 2; k++)
			{
				for (r = 0; r < 2; r++)
					N[k][r] += A[k] * A[r];
				U[k] += A[k] * L;
			}
			A[0] = 0;
			A[1] = 1;
			L = y2[i] - (Coefs[1] + y1[i]);
			for (k = 0; k < 2; k++)
			{
				for (r = 0; r < 2; r++)
					N[k][r] += A[k] * A[r];
				U[k] += A[k] * L;
			}
		}

		if (!INVSQR1((double*)N, U, 2))
			break;

		for (k = 0; k < 2; k++)
			Coefs[k] += U[k];

		if (abs(U[0]) < DELTA_LIMIT && abs(U[1]) < DELTA_LIMIT)
			break;

//		break; // theoretically iterations are not needed as a0 and b0 are independent

		iIterNum += 1;
	}

	memmove(U, Coefs, 2 * sizeof(double));
	for (i = 0; i < n; i++)
	{
		xdif = U[0] + x1[i] - x2[i];
		ydif = U[1] + y1[i] - y2[i];
		ptr_Errors[i] = sqrt(xdif*xdif + ydif*ydif);
	}

	RMSE1(ptr_Errors, n, 2, &fitting_rmse_, &errors_mean_);

	*ptr_errors_mean = errors_mean_;
	*ptr_fitting_rmse = fitting_rmse_;

	return;
}

/***
Fit 12 polymonial transformation coefficients in
x2 = a0 + a1*x1 + a2*y1 + a3*x1*x1 + a4*x1*y1 + a5*y1*y1
y2 = b0 + b1*x1 + b2*y1 + b3*x1*x1 + b4*x1*y1 + b5*y1*y1

Returned parameters
Coefs:			[a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5]
Errors:			fitting residuals for n points
errors_mean:	mean residual
fitting_rmse:	fitting RMSE
***/
void FitPolynomialTransform(double *x1, double *y1, double *x2, double *y2, int n, double Coefs[], double *ptr_Errors, double *ptr_errors_mean, double *ptr_fitting_rmse)
{
	double U[12];
	int iIterNum;
	int i, k, r;
	double L;
	double A[12], N[12][12];
	double xdif, ydif;
	double fitting_rmse_, errors_mean_;

	memset(ptr_Errors, 0, n * sizeof(double));
	*ptr_errors_mean = 0;
	*ptr_fitting_rmse = 0;

	if (n < 12)
		return;

	memset(Coefs, 0, 12 * sizeof(double));
	Coefs[1] = 1;
	Coefs[8] = 1;
	iIterNum = 0;
	while (iIterNum < 10)
	{
		memset(A, 0, 12 * sizeof(double)); // added 8/9/2021
		memset(N, 0, 144 * sizeof(double));
		memset(U, 0, 12 * sizeof(double));
		
		for (i = 0; i < n; i++)
		{

			// x2 = a0 + a1*x1 + a2*y1 + a3*x1*x1 + a4*x1*y1 + a5*y1*y1
			A[0] = 1;
			A[1] = x1[i];
			A[2] = y1[i];
			A[3] = x1[i] * x1[i];
			A[4] = x1[i] * y1[i];
			A[5] = y1[i] * y1[i];
			A[6] = 0;
			A[7] = 0;
			A[8] = 0;
			A[9] = 0;
			A[10] = 0;
			A[11] = 0;
			L = x2[i] - (Coefs[0] + Coefs[1] * x1[i] + Coefs[2] * y1[i] + Coefs[3] * x1[i] * x1[i] + Coefs[4] * x1[i] * y1[i] + Coefs[5] * y1[i] * y1[i]);
			for (k = 0; k < 12; k++)
			{
				for (r = 0; r < 12; r++)
					N[k][r] += A[k] * A[r];
				U[k] += A[k] * L;
			}

			// y2 = b0 + b1*x1 + b2*y1 + b3*x1*x1 + b4*x1*y1 + b5*y1*y1
			A[0] = 0;
			A[1] = 0;
			A[2] = 0;
			A[3] = 0;
			A[4] = 0;
			A[5] = 0;
			A[6] = 1;
			A[7] = x1[i];
			A[8] = y1[i];
			A[9] = x1[i] * x1[i];
			A[10] = x1[i] * y1[i];
			A[11] = y1[i] * y1[i];
			L = y2[i] - (Coefs[6] + Coefs[7] * x1[i] + Coefs[8] * y1[i] + Coefs[9] * x1[i] * x1[i] + Coefs[10] * x1[i] * y1[i] + Coefs[11] * y1[i] * y1[i]);
			for (k = 0; k < 12; k++)
			{
				for (r = 0; r < 12; r++)
					N[k][r] += A[k] * A[r];
				U[k] += A[k] * L;
			}
		}

		if (!INVSQR1((double*)N, U, 12))
			break;

		for (k = 0; k < 12; k++)
			Coefs[k] += U[k];

		if (abs(U[0]) < DELTA_LIMIT && abs(U[6]) < DELTA_LIMIT)
			break;

		iIterNum += 1;
	}

	memmove(U, Coefs, 12 * sizeof(double));
	for (i = 0; i < n; i++)
	{
		xdif = U[0] + U[1] * x1[i] + U[2] * y1[i] + U[3] * x1[i] * x1[i] + U[4] * x1[i] * y1[i] + U[5] * y1[i] * y1[i] - x2[i];
		ydif = U[6] + U[7] * x1[i] + U[8] * y1[i] + U[9] * x1[i] * x1[i] + U[10] * x1[i] * y1[i] + U[11] * y1[i] * y1[i] - y2[i];
		ptr_Errors[i] = sqrt(xdif*xdif + ydif*ydif); // fitting residual
	}

	RMSE1(ptr_Errors, n, 12, &fitting_rmse_, &errors_mean_);

	*ptr_errors_mean = errors_mean_;
	*ptr_fitting_rmse = fitting_rmse_;

	return;
}

int GetGridded_RBF_Centers(int iWidth, int iHeight, int iRBFs_K, double* pd_x_k, double* pd_y_k, int *pi_grid_width, int *pi_grid_height)
{
	int iSplit_K, iPtNum, iPtNum_diff;
	int grid_row, grid_col;
	int i;

	iSplit_K = (int)(sqrt(iRBFs_K * 1.0) + 1e-10);
	iPtNum = 0;
	for (grid_row = 0; grid_row < iSplit_K; grid_row++)
	{
		for (grid_col = 0; grid_col < iSplit_K; grid_col++)
		{
			if (iPtNum >= iRBFs_K)
				break;

			pd_x_k[iPtNum] = iWidth / iSplit_K * grid_col + iWidth / iSplit_K / 2;
			pd_y_k[iPtNum] = iHeight / iSplit_K * grid_row + iHeight / iSplit_K / 2;

			iPtNum += 1;
		}
	}

	iPtNum_diff = iRBFs_K - iPtNum;
	for (i = 0; i < iPtNum_diff; i++)
	{
		pd_x_k[iPtNum] = rand() * 1.0 / RAND_MAX * iWidth;
		pd_y_k[iPtNum] = rand() * 1.0 / RAND_MAX * iHeight;
		iPtNum += 1;
	}

	*pi_grid_width = iWidth / iSplit_K;
	*pi_grid_height = iHeight / iSplit_K;

	return iPtNum;
}
/***
Created 3/8/2019

Fit 12 polymonial transformation coefficients in
x2 = a0 + a1*x1 + a2*y1 + a3*x1*x1 + a4*x1*y1 + a5*y1*y1
y2 = b0 + b1*x1 + b2*y1 + b3*x1*x1 + b4*x1*y1 + b5*y1*y1

Returned parameters
Coefs:			[a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5]
Errors:			fitting residuals for n points
errors_mean:	mean residual
fitting_rmse:	fitting RMSE

Note: 12 + 2K coefficients are solved; additional 2K coefficients of RBFs centers' coordinates are also returned in Coefs
v1 (3/12/2019): gridded centers; input image width and height
v2 (2/12/2020): use grid fill number to remove K gridded centers; return new K
***/
int FitRBFsTransform_TPS_Poly_v2(double* x1, double* y1, double* x2, double* y2, int n, int K, int iWidth, int iHeight, double Coefs[], double* ptr_Errors, double* ptr_errors_mean, double* ptr_fitting_rmse)
{
	int iIterNum;
	int i, k, r, m;
	double L;
	double* U = NULL, * A = NULL, * N = NULL; // for N, only the first (12+2K)*(12+2K) elements are used
	double xdif, ydif;
	double fitting_rmse_, errors_mean_;

	int iCoefsNum;

	double wx_k[MAX_RBFs_K], wy_k[MAX_RBFs_K], x_k[MAX_RBFs_K], y_k[MAX_RBFs_K]; // 4*K coefficients
	double D_k[MAX_RBFs_K], RBFs_k[MAX_RBFs_K], RBFs_sum;

	int iSplit_K;
	int grid_row, grid_col;
	int iPtNum, iPtNum_diff;

	// v2
	unsigned int iTiePointsNum_k[MAX_RBFs_K]; // number of tie points per grid
	int grid_width, grid_height;
	unsigned int tie_grid_threshold = RBF_TIE_NUMBER_THRESHOLD;
	double x_k_new[MAX_RBFs_K], y_k_new[MAX_RBFs_K];
	int K_new;

	if (K > MAX_RBFs_K)
		K = MAX_RBFs_K;

	if (n < (2 * K + 12))
		return 0;

	iCoefsNum = 12 + 2 * K;

	if (!(U = (double*)calloc(iCoefsNum, sizeof(double)))
		|| !(A = (double*)calloc(iCoefsNum, sizeof(double)))
		|| !(N = (double*)calloc(iCoefsNum * iCoefsNum, sizeof(double))))
	{
		printf("\nError in FitRBFsTransform_TPS_Poly_v2(): insufficient memory.\n");
		scanf_s(" %d", &K);
		exit(1);
	}

	memset(ptr_Errors, 0, n * sizeof(double));
	*ptr_errors_mean = 0;
	*ptr_fitting_rmse = 0;

	memset(Coefs, 0, iCoefsNum * sizeof(double));

	// assign initial coefficients values
	// 0 ~ 11: a0, a1, a2, b0, b1, b2
	Coefs[1] = 1; // a1
	Coefs[8] = 1; // b2

	// 12 ~ 5 + 4*K: wxk[1], wyk[1]; wxk[2], wyk[2], ...
	iSplit_K = (int)(sqrt(K * 1.0) + 1e-10);

	iPtNum = 0;
	for (grid_row = 0; grid_row < iSplit_K; grid_row++)
	{
		for (grid_col = 0; grid_col < iSplit_K; grid_col++)
		{
			if (iPtNum >= K)
				break;

			x_k[iPtNum] = iWidth / iSplit_K * grid_col + iWidth / iSplit_K / 2;
			y_k[iPtNum] = iHeight / iSplit_K * grid_row + iHeight / iSplit_K / 2;

			iPtNum += 1;
		}
	}
	iPtNum_diff = K - iPtNum;
	for (i = 0; i < iPtNum_diff; i++)
	{
		x_k[iPtNum] = rand() * 1.0 / RAND_MAX * iWidth;
		y_k[iPtNum] = rand() * 1.0 / RAND_MAX * iHeight;
		iPtNum += 1;
	}

	// v2
	grid_width = iWidth / iSplit_K;
	grid_height = iHeight / iSplit_K;
	// count tie points per grid]
	memset(iTiePointsNum_k, 0, MAX_RBFs_K * sizeof(unsigned int));
	for (i = 0; i < n; i++)
	{
		for (m = 0; m < K; m++)
		{
			if (ABS(x1[i] - x_k[m]) < grid_width / 2 && ABS(y1[i] - y_k[m]) < grid_height / 2)
				iTiePointsNum_k[m] += 1;
		}
	}
	// set new control points
	K_new = 0;
	for (m = 0; m < K; m++)
	{
		if (iTiePointsNum_k[m] >= tie_grid_threshold)
		{
			x_k_new[K_new] = x_k[m];
			y_k_new[K_new] = y_k[m];
			K_new += 1;
		}
	}
	memset(x_k, 0, MAX_RBFs_K * sizeof(double));
	memset(y_k, 0, MAX_RBFs_K * sizeof(double));
	K = K_new;
	memmove(x_k, x_k_new, K_new * sizeof(double));
	memmove(y_k, y_k_new, K_new * sizeof(double));

	iCoefsNum = 12 + 2 * K;

	iIterNum = 0;
	while (iIterNum < 10)
	{
		memset(N, 0, iCoefsNum * iCoefsNum * sizeof(double)); // added 8/9/2021
		memset(U, 0, iCoefsNum * sizeof(double));

		// get kernal points (centers) and weights
		for (m = 0; m < K; m++)
		{
			wx_k[m] = Coefs[12 + 2 * m];
			wy_k[m] = Coefs[12 + 2 * m + 1];
		}
		for (i = 0; i < n; i++)
		{
			// calcualte squared distances to kernal points and RBFs function values
			for (m = 0; m < K; m++)
			{
				D_k[m] = (x_k[m] - x1[i]) * (x_k[m] - x1[i]) + (y_k[m] - y1[i]) * (y_k[m] - y1[i]);
				RBFs_k[m] = (D_k[m] < 1e-20) ? 0 : D_k[m] * log(sqrt(D_k[m]));// exp(-D_k[m]);
			}

			// x2 = a0 + a1*x1 + a2*y1 + wx1*exp(-Di1^2) + wx2*exp(-Di2^2) + wx3*exp(-Di3^2) + ...
			// Di1^2 = (x1 - xk1)^2 + (y1 - yk1)^2
			memset(A, 0, (12 + 2 * K) * sizeof(double));
			A[0] = 1;
			A[1] = x1[i];
			A[2] = y1[i];
			A[3] = x1[i] * x1[i];
			A[4] = x1[i] * y1[i];
			A[5] = y1[i] * y1[i];
			A[6] = 0;
			A[7] = 0;
			A[8] = 0;
			A[9] = 0;
			A[10] = 0;
			A[11] = 0;
			for (m = 0; m < K; m++)
				A[12 + 2 * m] = RBFs_k[m]; // df / d(wxk)

			RBFs_sum = 0;
			for (m = 0; m < K; m++)
				RBFs_sum += wx_k[m] * RBFs_k[m];

			L = x2[i] - (Coefs[0] + Coefs[1] * x1[i] + Coefs[2] * y1[i] + Coefs[3] * x1[i] * x1[i] + Coefs[4] * x1[i] * y1[i] + Coefs[5] * y1[i] * y1[i] + RBFs_sum);

			for (k = 0; k < iCoefsNum; k++)
			{
				for (r = 0; r < iCoefsNum; r++)
					N[k * iCoefsNum + r] += A[k] * A[r];
				U[k] += A[k] * L;
			}

			// y2 = b0 + b1*x1 + b2*y1 + wy1*exp(-Di1^2) + wy2*exp(-Di2^2) + wy3*exp(-Di3^2) + ...
			memset(A, 0, (12 + 2 * K) * sizeof(double));
			A[0] = 0;
			A[1] = 0;
			A[2] = 0;
			A[3] = 0;
			A[4] = 0;
			A[5] = 0;
			A[6] = 1;
			A[7] = x1[i];
			A[8] = y1[i];
			A[9] = x1[i] * x1[i];
			A[10] = x1[i] * y1[i];
			A[11] = y1[i] * y1[i];
			for (m = 0; m < K; m++)
				A[12 + 2 * m + 1] = RBFs_k[m]; // df / d(wyk)

			RBFs_sum = 0;
			for (m = 0; m < K; m++)
				RBFs_sum += wy_k[m] * RBFs_k[m];

			L = y2[i] - (Coefs[6] + Coefs[7] * x1[i] + Coefs[8] * y1[i] + Coefs[9] * x1[i] * x1[i] + Coefs[10] * x1[i] * y1[i] + Coefs[11] * y1[i] * y1[i] + RBFs_sum);

			for (k = 0; k < iCoefsNum; k++)
			{
				for (r = 0; r < iCoefsNum; r++)
					N[k * iCoefsNum + r] += A[k] * A[r];
				U[k] += A[k] * L;
			}
		}

		if (!INVSQR1(N, U, iCoefsNum))
			break;

		for (k = 0; k < iCoefsNum; k++)
			Coefs[k] += U[k];

		if (abs(U[0]) < DELTA_LIMIT && abs(U[6]) < DELTA_LIMIT)
			break;

		iIterNum += 1;
	}

	memmove(U, Coefs, iCoefsNum * sizeof(double));
	for (m = 0; m < K; m++)
	{
		wx_k[m] = Coefs[12 + 2 * m];
		wy_k[m] = Coefs[12 + 2 * m + 1];
	}
	for (i = 0; i < n; i++)
	{
		for (m = 0; m < K; m++)
		{
			D_k[m] = (x_k[m] - x1[i]) * (x_k[m] - x1[i]) + (y_k[m] - y1[i]) * (y_k[m] - y1[i]);
			RBFs_k[m] = (D_k[m] < 1e-20) ? 0 : D_k[m] * log(sqrt(D_k[m]));
		}

		RBFs_sum = 0;
		for (m = 0; m < K; m++)
			RBFs_sum += wx_k[m] * RBFs_k[m];

		xdif = U[0] + U[1] * x1[i] + U[2] * y1[i] + U[3] * x1[i] * x1[i] + U[4] * x1[i] * y1[i] + U[5] * y1[i] * y1[i] + RBFs_sum - x2[i];

		RBFs_sum = 0;
		for (m = 0; m < K; m++)
			RBFs_sum += wy_k[m] * RBFs_k[m];

		ydif = U[6] + U[7] * x1[i] + U[8] * y1[i] + U[9] * x1[i] * x1[i] + U[10] * x1[i] * y1[i] + U[11] * y1[i] * y1[i] + RBFs_sum - y2[i];

		ptr_Errors[i] = sqrt(xdif * xdif + ydif * ydif); // fitting residual
	}

	RMSE1(ptr_Errors, n, iCoefsNum, &fitting_rmse_, &errors_mean_); //

	*ptr_errors_mean = errors_mean_;
	*ptr_fitting_rmse = fitting_rmse_;

	// save RBFs xk and ys to Coefs
	for (m = 0; m < K; m++)
	{
		Coefs[12 + 4 * m] = wx_k[m];
		Coefs[12 + 4 * m + 1] = wy_k[m];
		Coefs[12 + 4 * m + 2] = x_k[m];
		Coefs[12 + 4 * m + 3] = y_k[m];
	}

	free(U);
	free(A);
	free(N);

	return K_new;
}

/***
Get transformation coefficients for f() in
(x2, y2) = f(x1, y1) where f() is the transformation.

Parameters
iTransformationType:	1 - translation, 2 - affine, 3 - polynomial
Coefs[12]:				transformation coefficients, [a0, b0] for translation,
						[a0, a1, a2, b0, b1, b2] for affine, and
						[a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5] for polynomial
***/
void GetTransformedCoords(double x1, double y1, int iTransformationType, double *Coefs, double *ptr_x2, double *ptr_y2, int iRBFs_K)
{
	double x2_, y2_;
	double D_k, RBFs_k, RBFs_x_sum, RBFs_y_sum;
	int m;
	double x_k, y_k, wx_k, wy_k;

	*ptr_x2 = 0;
	*ptr_y2 = 0;

	if (iTransformationType == 4 && iRBFs_K == 0)
	{
		printf("\nError in GetTransformedCoords(): RBFs with K = 0.\n");
		scanf_s(" %d", &m);
		exit(2);
	}

	switch (iTransformationType)
	{
	case 1:
		x2_ = Coefs[0] + x1;
		y2_ = Coefs[1] + y1;
		break;
	case 2:
		x2_ = Coefs[0] + Coefs[1] * x1 + Coefs[2] * y1;
		y2_ = Coefs[3] + Coefs[4] * x1 + Coefs[5] * y1;
		break;
	case 3:
		x2_ = Coefs[0] + Coefs[1] * x1 + Coefs[2] * y1 + Coefs[3] * x1*x1 + Coefs[4] * x1*y1 + Coefs[5] * y1*y1;
		y2_ = Coefs[6] + Coefs[7] * x1 + Coefs[8] * y1 + Coefs[9] * x1*x1 + Coefs[10] * x1*y1 + Coefs[11] * y1*y1;
		break;
	case 4:
		RBFs_x_sum = 0;
		RBFs_y_sum = 0;
		for (m = 0; m < iRBFs_K; m++)
		{
			wx_k = Coefs[12 + 4 * m];
			wy_k = Coefs[12 + 4 * m + 1];
			x_k = Coefs[12 + 4 * m + 2];
			y_k = Coefs[12 + 4 * m + 3];

			D_k = (x_k - x1) * (x_k - x1) + (y_k - y1) * (y_k - y1);
			if (D_k < 1e-20)
				continue;
			RBFs_k = D_k * log(sqrt(D_k));
			RBFs_x_sum += wx_k * RBFs_k;
			RBFs_y_sum += wy_k * RBFs_k;
		}
		x2_ = Coefs[0] + Coefs[1] * x1 + Coefs[2] * y1 + Coefs[3] * x1 * x1 + Coefs[4] * x1 * y1 + Coefs[5] * y1 * y1 + RBFs_x_sum;
		y2_ = Coefs[6] + Coefs[7] * x1 + Coefs[8] * y1 + Coefs[9] * x1 * x1 + Coefs[10] * x1 * y1 + Coefs[11] * y1 * y1 + RBFs_y_sum;
		break;
	default:
		x2_ = Coefs[0] + Coefs[1] * x1 + Coefs[2] * y1;
		y2_ = Coefs[3] + Coefs[4] * x1 + Coefs[5] * y1;
		break;
	}

	*ptr_x2 = x2_;
	*ptr_y2 = y2_;

	return;
}

/***
iTransformationType: TRANSLATION = 1, AFFINE = 2, POLYNOMIAL = 3
For MSS processing
***/
void TransformImage_sht2char(short* pshtData, int iColMax, int iRowMax, int iTransformationType, double* pdCoefs, unsigned char* pucData_trsf)
{
	int Row, Col;

	double dTargetRow, dTargetCol;
	int iTargetRow, iTargetCol;

	float fTargetValue;
	float fWeight, fWeightSum;
	float fDX, fDY;
	short int shtTargetValue;

	memset(pucData_trsf, 0, iRowMax * iColMax * sizeof(unsigned char));
	for (Row = 0; Row < iRowMax; Row++)
	{
		for (Col = 0; Col < iColMax; Col++)
		{
			GetTransformedCoords((double)Col, (double)Row, iTransformationType, pdCoefs, &dTargetCol, &dTargetRow);

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
				fTargetValue += pshtData[iTargetRow * iColMax + iTargetCol] * fWeight;
				fWeightSum += fWeight;
			}

			// top right point
			iTargetCol += 1;
			if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
			{
				fWeight = fDX * (1 - fDY);
				fTargetValue += pshtData[iTargetRow * iColMax + iTargetCol] * fWeight;
				fWeightSum += fWeight;
			}

			// bottom right point
			iTargetRow += 1;
			if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
			{
				fWeight = fDX * fDY;
				fTargetValue += pshtData[iTargetRow * iColMax + iTargetCol] * fWeight;
				fWeightSum += fWeight;
			}

			// bottom left point
			iTargetCol -= 1;
			if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
			{
				fWeight = (1 - fDX) * fDY;
				fTargetValue += pshtData[iTargetRow * iColMax + iTargetCol] * fWeight;
				fWeightSum += fWeight;
			}

			if (fWeightSum > 0)
			{
				fTargetValue /= fWeightSum;
				shtTargetValue = (short)(fTargetValue + 0.5f);
				shtTargetValue = (shtTargetValue >= 0) ? shtTargetValue : 0;
				shtTargetValue = (shtTargetValue <= 255) ? shtTargetValue : 255;
				pucData_trsf[Row * iColMax + Col] = (unsigned char)shtTargetValue;
			}
		}
	}

	return;
}

/***
iTransformationType: TRANSLATION = 1, AFFINE = 2, POLYNOMIAL = 3, RBFs = 4
For MSS processing
***/
void TransformImage_UC(unsigned char* pucData, int iColMax, int iRowMax, int iTransformationType, double* pdCoefs, unsigned char* pucData_trsf, int iRBFs_K)
{
	int Row, Col;

	double dTargetRow, dTargetCol;
	int iTargetRow, iTargetCol;

	float fTargetValue;
	float fWeight, fWeightSum;
	float fDX, fDY;
	short int shtTargetValue;

	memset(pucData_trsf, 0, iRowMax * iColMax * sizeof(unsigned char));
	for (Row = 0; Row < iRowMax; Row++)
	{
		for (Col = 0; Col < iColMax; Col++)
		{
			GetTransformedCoords((double)Col, (double)Row, iTransformationType, pdCoefs, &dTargetCol, &dTargetRow, iRBFs_K);

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
				fTargetValue += pucData[iTargetRow * iColMax + iTargetCol] * fWeight;
				fWeightSum += fWeight;
			}

			// top right point
			iTargetCol += 1;
			if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
			{
				fWeight = fDX * (1 - fDY);
				fTargetValue += pucData[iTargetRow * iColMax + iTargetCol] * fWeight;
				fWeightSum += fWeight;
			}

			// bottom right point
			iTargetRow += 1;
			if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
			{
				fWeight = fDX * fDY;
				fTargetValue += pucData[iTargetRow * iColMax + iTargetCol] * fWeight;
				fWeightSum += fWeight;
			}

			// bottom left point
			iTargetCol -= 1;
			if (iTargetCol >= 0 && iTargetCol < iColMax && iTargetRow >= 0 && iTargetRow < iRowMax)
			{
				fWeight = (1 - fDX) * fDY;
				fTargetValue += pucData[iTargetRow * iColMax + iTargetCol] * fWeight;
				fWeightSum += fWeight;
			}

			if (fWeightSum > 0)
			{
				fTargetValue /= fWeightSum;
				shtTargetValue = (short)(fTargetValue + 0.5f);
				shtTargetValue = (shtTargetValue >= 0) ? shtTargetValue : 0;
				shtTargetValue = (shtTargetValue <= 255) ? shtTargetValue : 255;
				pucData_trsf[Row * iColMax + Col] = (unsigned char)shtTargetValue;
			}
		}
	}

	return;
}

void CopySht2UChar(short* pshtData, int iSize, unsigned char* pucData_cpy)
{
	int i;
	unsigned char* puc = NULL;
	short int* psht = NULL;
	short int shtTargetValue;

	psht = pshtData;
	puc = pucData_cpy;
	for (i = 0; i < iSize; i++)
	{
		shtTargetValue = *psht++;
		shtTargetValue = (shtTargetValue >= 0) ? shtTargetValue : 0;
		shtTargetValue = (shtTargetValue <= 255) ? shtTargetValue : 255;
		*puc = (unsigned char)shtTargetValue;
		puc++;
	}
	return;
}

int imsub(short int* im, int ncol, int nrow, int x, int y, int h, short int* sub)
{
	int i, j;
	int w;

	w = h * 2 + 1;
	memset(sub, 0, w * w * sizeof(short int));

	if (x<h || x>ncol - h - 1 || y<h || y>nrow - h - 1)
		return 0;

	for (i = 0; i < w; i++)
	{
		for (j = 0; j < w; j++)
			sub[i * w + j] = im[(y - h + i) * ncol + (x - h + j)];
	}

	return 1;
}

/***
correlation coefficient
***/
double corr2(short int *x, short int *y, int n)
{
	double	dSumL1 = 0, dSumL2 = 0, dSumR1 = 0, dSumR2 = 0, dSumLR = 0;
	short int *pSrc, *pObj;
	int     i1;
	double  dCoefL, dCoef1, dCoef2, dCoef;
	int     MATCH_FULL_SIZE;
	short int src, obj;

	pSrc = x;
	pObj = y;
	MATCH_FULL_SIZE = n;
	for (i1 = 0; i1<MATCH_FULL_SIZE; i1++)
	{
		src = (*pSrc);
		obj = (*pObj);
		dSumL1 += (src);
		dSumL2 += (double)((src)*(src));
		dSumR1 += (obj);
		dSumR2 += (double)((obj)* (obj));
		dSumLR += (double)((src)* (obj));
		pSrc++;
		pObj++;
	}

	dCoefL = dSumL2 - dSumL1*dSumL1 / MATCH_FULL_SIZE;
	dCoef1 = dSumLR - dSumL1*dSumR1 / MATCH_FULL_SIZE;
	dCoef2 = dCoefL*(dSumR2 - dSumR1*dSumR1 / MATCH_FULL_SIZE);
	if (dCoef2<0.0000001 && dCoef2>-0.0000001)
		dCoef = 0.0;
	else
		dCoef = dCoef1 / sqrt(dCoef2);

	return dCoef;
}

/***
correlation coefficient
Created 4/20/2021
Skip fill values
***/
double corr2_UC(unsigned char* x, unsigned char* y, int n, unsigned char fillvalue)
{
	double	dSumL1 = 0, dSumL2 = 0, dSumR1 = 0, dSumR2 = 0, dSumLR = 0;
	unsigned char* pSrc = NULL, * pObj = NULL;
	int     i1;
	double  dCoefL, dCoef1, dCoef2, dCoef;
	int     MATCH_FULL_SIZE;
	unsigned char src, obj;

	pSrc = x;
	pObj = y;
	MATCH_FULL_SIZE = n;
	for (i1 = 0; i1 < MATCH_FULL_SIZE; i1++)
	{
		src = (*pSrc);
		obj = (*pObj);
		if (src != fillvalue && obj != fillvalue)
		{
			dSumL1 += (src);
			dSumL2 += ((double)(src)) * (src);
			dSumR1 += (obj);
			dSumR2 += ((double)(obj)) * (obj);
			dSumLR += ((double)(src)) * (obj);
		}
		pSrc++;
		pObj++;
	}

	dCoefL = dSumL2 - dSumL1 * dSumL1 / MATCH_FULL_SIZE;
	dCoef1 = dSumLR - dSumL1 * dSumR1 / MATCH_FULL_SIZE;
	dCoef2 = dCoefL * (dSumR2 - dSumR1 * dSumR1 / MATCH_FULL_SIZE);
	if (dCoef2<0.0000001 && dCoef2>-0.0000001)
		dCoef = 0.0;
	else
		dCoef = dCoef1 / sqrt(dCoef2);

	return dCoef;
}

void CalcGradient2D(float* pImg, float* pImgX, float* pImgY, int iWidth, int iHeight)
{
	int i1, j1;
	for (i1 = 1; i1 < iHeight - 1; i1++)
	{
		for (j1 = 1; j1 < iWidth - 1; j1++)
		{
			pImgX[i1*iWidth + j1] = (pImg[i1*iWidth + j1 + 1] - pImg[i1*iWidth + j1 - 1]) / 2;
			pImgY[i1*iWidth + j1] = (pImg[(i1 + 1)*iWidth + j1] - pImg[(i1 - 1)*iWidth + j1]) / 2;
		}
	}
}

void CalcMultiply(float* pImg1, float* pImg2, float* pImgNew, int iSize)
{
	int i1;
	for (i1 = 0; i1 < iSize; i1++)
		pImgNew[i1] = pImg1[i1] * pImg2[i1];
}

void CalcDivide(float* pImg1, float* pImg2, float* pImgNew, int iSize)
{
	int i1;
	for (i1 = 0; i1 < iSize; i1++)
		pImgNew[i1] = (ABS(pImg2[i1]) > 0.1e-10f) ? pImg1[i1] / pImg2[i1] : 0;
}

void CalcAdd(float* pImg1, float* pImg2, float* pImgNew, int iSize)
{
	int i1;
	for (i1 = 0; i1 < iSize; i1++)
		pImgNew[i1] = pImg1[i1] + pImg2[i1];
}

void CalcSubtract(float* pImg1, float* pImg2, float* pImgNew, int iSize)
{
	int i1;
	for (i1 = 0; i1 < iSize; i1++)
		pImgNew[i1] = pImg1[i1] - pImg2[i1];
}

void CalcMinMaxMean(float* pImg, float *pdMin, float *pdMax, float *pdMean, int iSize)
{
	int iValidNum;
	int i1;

	*pdMin = 1000000000;
	*pdMax = -1000000000;
	*pdMean = 0;

	iValidNum = 0;

	for (i1 = 0; i1 < iSize; i1++)
	{
		if (ABS(pImg[i1]) < 0.1e-10)
			continue;

		*pdMin = *pdMin < pImg[i1] ? *pdMin : pImg[i1];
		*pdMax = *pdMax > pImg[i1] ? *pdMax : pImg[i1];
		*pdMean += pImg[i1];
		iValidNum += 1;
	}
	if (iValidNum > 0)
		*pdMean /= iValidNum;
	else
		*pdMean = 0;

	return;
}

void CalcMinMaxMeanWithMask(float* pImg, unsigned char *pucMask, float *pdMin, float *pdMax, float *pdMean, int iSize)
{
	int iValidNum;
	int i1;

	*pdMin = 1000000000;
	*pdMax = -1000000000;
	*pdMean = 0;

	iValidNum = 0;

	for (i1 = 0; i1 < iSize; i1++)
	{
		if (ABS(pImg[i1]) < 0.1e-10)
			continue;
		if (pucMask[i1] == 0)
			continue;

		*pdMin = *pdMin < pImg[i1] ? *pdMin : pImg[i1];
		*pdMax = *pdMax > pImg[i1] ? *pdMax : pImg[i1];
		*pdMean += pImg[i1];
		iValidNum += 1;
	}
	if (iValidNum > 0)
		*pdMean /= iValidNum;
	else
		*pdMean = 0;

	return;
}

void ZeroBuffer(float *pBuffer, int iSize)
{
	int i1;
	float* pPntr = pBuffer;
	for (i1 = 0; i1 < iSize; i1++)
		*(pPntr++) = 0;
}

void AddConst(float *pBuffer, float c, int iSize)
{
	int i;
	for (i = 0; i < iSize; i++)
		pBuffer[i] = pBuffer[i] + c;
}

// get a Gaussian Kernel
float* GetGaussian(double dSigma, int *iFilterWidth)
{
	int		i, j;
	int		iWidth = (int)(3 * sqrt(8.0)*dSigma + 1e-10);
	if ((iWidth % 2) == 0)	iWidth = iWidth - 1; // to get the odd size
	*iFilterWidth = iWidth;
	int		iHalfWidth = (iWidth - 1) / 2, iSize = iWidth*iWidth;

	float* dBuffer = (float *)calloc(iSize, sizeof(float));
	float* pdBuffer = dBuffer;
	double dSum = 0;
	double dRef = -0.5 / pow(dSigma, 2);

	for (i = -iHalfWidth; i <= iHalfWidth; i++)
	{
		for (j = -iHalfWidth; j <= iHalfWidth; j++)
		{
			*pdBuffer = (float)(exp(dRef * (i*i + j*j)));
			dSum += *pdBuffer;
			pdBuffer++;
		}
	}

	pdBuffer = dBuffer;
	for (i = 0; i < iSize; i++)
	{
		*pdBuffer = (float)(*pdBuffer / dSum);
		pdBuffer++;
	}
	return dBuffer;
}

void Conv2same(short int *pImg, short int *pImgNew, int iWidth, int iHeight, float* dFilter, int w)
{
	int		h;
	int		iFilterSize;
	int		nCol, nRow;
	int		nSize;
	float	dCurSum, fsum;
	int		i, j, k, l;
	int k1, k2, l1, l2;

	h = (w - 1) / 2;
	nCol = iWidth + 2 * h, nRow = iHeight + 2 * h;
	iFilterSize = w*w;
	nSize = nRow*nCol;

	memset(pImgNew, 0, iWidth*iHeight*sizeof(short int));

	for (i = h; i < iHeight - h; i++) // old: for (i = 0; i <= iHeight - w; i++)
	{
		for (j = h; j < iWidth - h; j++) // old: for (j = 0; j <= iWidth - w; j++)
		{
			if (ISFILLVALUE(pImg[i*iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			for (k = 0; k<w; k++)
			{
				for (l = 0; l < w; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h)*iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h)*iWidth + (j + l - h)] * dFilter[k*w + l];	// old: dCurSum += pImg[(i + k)*iWidth + (j + l)] * dFilter[k*w + l];
					fsum += dFilter[k*w + l];
				}
			}

			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i*iWidth + j] = (short int)(dCurSum + 1e-10);
		}
	}

	// Upper part
	for (i = 0; i<h; i++)
	{
		for (j = 0; j<iWidth; j++)
		{
			if (ISFILLVALUE(pImg[i*iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			k1 = (h - i)>0 ? (h - i) : 0;
			k2 = (iHeight - i)<w ? (iHeight - i) : w;
			l1 = (h - j)>0 ? (h - j) : 0;
			l2 = (iWidth - j)<w ? (iWidth - j) : w;
			for (k = k1; k<k2; k++)
			{
				for (l = l1; l<l2; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h)*iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h)*iWidth + (j + l - h)] * dFilter[k*w + l];
					fsum += dFilter[k*w + l];
				}
			}
			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i*iWidth + j] = (short int)(dCurSum + 1e-10);
		}
	}

	// lower part
	for (i = iHeight - h; i<iHeight; i++)
	{
		for (j = 0; j<iWidth; j++)
		{
			if (ISFILLVALUE(pImg[i*iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			k1 = (h - i)>0 ? (h - i) : 0;
			k2 = (iHeight - i)<w ? (iHeight - i) : w;
			l1 = (h - j)>0 ? (h - j) : 0;
			l2 = (iWidth - j)<w ? (iWidth - j) : w;
			for (k = k1; k<k2; k++)
			{
				for (l = l1; l<l2; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h)*iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h)*iWidth + (j + l - h)] * dFilter[k*w + l];
					fsum += dFilter[k*w + l];
				}
			}
			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i*iWidth + j] = (short int)(dCurSum + 1e-10);
		}
	}
	// left part
	for (i = 0; i<iHeight; i++)
	{
		for (j = 0; j<h; j++)
		{
			if (ISFILLVALUE(pImg[i*iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			k1 = (h - i)>0 ? (h - i) : 0;
			k2 = (iHeight - i)<w ? (iHeight - i) : w;
			l1 = (h - j)>0 ? (h - j) : 0;
			l2 = (iWidth - j)<w ? (iWidth - j) : w;
			for (k = k1; k<k2; k++)
			{
				for (l = l1; l<l2; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h)*iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h)*iWidth + (j + l - h)] * dFilter[k*w + l];
					fsum += dFilter[k*w + l];
				}
			}
			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i*iWidth + j] = (short int)(dCurSum + 1e-10);
		}
	}
	// right part
	for (i = 0; i<iHeight; i++)
	{
		for (j = iWidth - h; j<iWidth; j++)
		{
			if (ISFILLVALUE(pImg[i*iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			k1 = (h - i)>0 ? (h - i) : 0;
			k2 = (iHeight - i)<w ? (iHeight - i) : w;
			l1 = (h - j)>0 ? (h - j) : 0;
			l2 = (iWidth - j)<w ? (iWidth - j) : w;
			for (k = k1; k<k2; k++)
			{
				for (l = l1; l<l2; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h)*iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h)*iWidth + (j + l - h)] * dFilter[k*w + l];
					fsum += dFilter[k*w + l];
				}
			}
			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i*iWidth + j] = (short int)(dCurSum + 1e-10);
		}
	}

	return;
}

void Conv2same_UC(unsigned char* pImg, unsigned char* pImgNew, int iWidth, int iHeight, float* dFilter, int w)
{
	int		h;
	int		iFilterSize;
	int		nCol, nRow;
	int		nSize;
	float	dCurSum, fsum;
	int		i, j, k, l;
	int k1, k2, l1, l2;

	h = (w - 1) / 2;
	nCol = iWidth + 2 * h, nRow = iHeight + 2 * h;
	iFilterSize = w * w;
	nSize = nRow * nCol;

	memset(pImgNew, 0, iWidth * iHeight * sizeof(unsigned char));

	for (i = h; i < iHeight - h; i++) // old: for (i = 0; i <= iHeight - w; i++)
	{
		for (j = h; j < iWidth - h; j++) // old: for (j = 0; j <= iWidth - w; j++)
		{
			if (ISFILLVALUE(pImg[i * iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			for (k = 0; k < w; k++)
			{
				for (l = 0; l < w; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h) * iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h) * iWidth + (j + l - h)] * dFilter[k * w + l];	// old: dCurSum += pImg[(i + k)*iWidth + (j + l)] * dFilter[k*w + l];
					fsum += dFilter[k * w + l];
				}
			}

			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i * iWidth + j] = (unsigned char)(dCurSum + 1e-10);
		}
	}

	// Upper part
	for (i = 0; i < h; i++)
	{
		for (j = 0; j < iWidth; j++)
		{
			if (ISFILLVALUE(pImg[i * iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			k1 = (h - i) > 0 ? (h - i) : 0;
			k2 = (iHeight - i) < w ? (iHeight - i) : w;
			l1 = (h - j) > 0 ? (h - j) : 0;
			l2 = (iWidth - j) < w ? (iWidth - j) : w;
			for (k = k1; k < k2; k++)
			{
				for (l = l1; l < l2; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h) * iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h) * iWidth + (j + l - h)] * dFilter[k * w + l];
					fsum += dFilter[k * w + l];
				}
			}
			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i * iWidth + j] = (unsigned char)(dCurSum + 1e-10);
		}
	}

	// lower part
	for (i = iHeight - h; i < iHeight; i++)
	{
		for (j = 0; j < iWidth; j++)
		{
			if (ISFILLVALUE(pImg[i * iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			k1 = (h - i) > 0 ? (h - i) : 0;
			k2 = (iHeight - i) < w ? (iHeight - i) : w;
			l1 = (h - j) > 0 ? (h - j) : 0;
			l2 = (iWidth - j) < w ? (iWidth - j) : w;
			for (k = k1; k < k2; k++)
			{
				for (l = l1; l < l2; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h) * iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h) * iWidth + (j + l - h)] * dFilter[k * w + l];
					fsum += dFilter[k * w + l];
				}
			}
			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i * iWidth + j] = (unsigned char)(dCurSum + 1e-10);
		}
	}
	// left part
	for (i = 0; i < iHeight; i++)
	{
		for (j = 0; j < h; j++)
		{
			if (ISFILLVALUE(pImg[i * iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			k1 = (h - i) > 0 ? (h - i) : 0;
			k2 = (iHeight - i) < w ? (iHeight - i) : w;
			l1 = (h - j) > 0 ? (h - j) : 0;
			l2 = (iWidth - j) < w ? (iWidth - j) : w;
			for (k = k1; k < k2; k++)
			{
				for (l = l1; l < l2; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h) * iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h) * iWidth + (j + l - h)] * dFilter[k * w + l];
					fsum += dFilter[k * w + l];
				}
			}
			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i * iWidth + j] = (unsigned char)(dCurSum + 1e-10);
		}
	}
	// right part
	for (i = 0; i < iHeight; i++)
	{
		for (j = iWidth - h; j < iWidth; j++)
		{
			if (ISFILLVALUE(pImg[i * iWidth + j]))
				continue;

			dCurSum = 0;
			fsum = 0;
			k1 = (h - i) > 0 ? (h - i) : 0;
			k2 = (iHeight - i) < w ? (iHeight - i) : w;
			l1 = (h - j) > 0 ? (h - j) : 0;
			l2 = (iWidth - j) < w ? (iWidth - j) : w;
			for (k = k1; k < k2; k++)
			{
				for (l = l1; l < l2; l++)
				{
					if (ISFILLVALUE(pImg[(i + k - h) * iWidth + (j + l - h)]))
						continue;

					dCurSum += pImg[(i + k - h) * iWidth + (j + l - h)] * dFilter[k * w + l];
					fsum += dFilter[k * w + l];
				}
			}
			if (ABS(fsum) < 1e-10)
				continue;

			dCurSum /= fsum;
			pImgNew[i * iWidth + j] = (unsigned char)(dCurSum + 1e-10);
		}
	}

	return;
}

/***
Read image data.
Rreturn NULL if cannot open file
Parameters
iDataType:	= 0 - unsigned char, 1 - short int, 2 - int, 3 - float, 4 - int, 5 - double
11/3/2020: use unsigned long data type for memory size
***/
void *ReadInputImage(char *pacImageFileName, int iColMax, int iRowMax, int iDataType)
{
	char cChar;
	int iDataSize;
	void *pImageData = NULL;
	FILE *fin = NULL;

	fopen_s(&fin, pacImageFileName, "rb");
	if (fin == NULL)
	{
		printf("\nError in ReadInputImage(): cannot open %s.\n", pacImageFileName);
		scanf_s(" %c", &cChar, 1);
		exit(1);
		return NULL;
	}

	switch (iDataType)
	{
	case 0:
		iDataSize = sizeof(unsigned char);
		break;
	case 1:
		iDataSize = sizeof(short int);
		break;
	case 3:
		iDataSize = sizeof(float);
		break;
	case 2:
	case 4:
		iDataSize = sizeof(int);
		break;
	default:
		iDataSize = sizeof(double);
	}

	if (!(pImageData = (void *)calloc(((unsigned long)iRowMax) * iColMax, iDataSize)))
	{
		printf("\nError in ReadInputImage: insufficient memory.\n");
		scanf_s(" %s", &cChar, 1);
		exit(1);
	}

	fread((char *)pImageData, iDataSize, ((unsigned long)iRowMax) * iColMax, fin);

	fclose(fin);

	return (void*)(pImageData);
}

/***
Read image data.
Rreturn NULL if cannot open file
Parameters
iDataType:	= 0 - unsigned char, 1 - short int, 2 - int, 3 - float, 4 - int, 5 - double
11/3/2020: use unsigned long data type for memory size
v1 (12/28/2020): return null if file not exist
***/
void* ReadInputImage_v1(char* pacImageFileName, int iColMax, int iRowMax, int iDataType)
{
	char cChar;
	int iDataSize;
	void* pImageData = NULL;
	FILE* fin = NULL;

	fopen_s(&fin, pacImageFileName, "rb");
	if (fin == NULL)
		return NULL;

	switch (iDataType)
	{
	case 0:
		iDataSize = sizeof(unsigned char);
		break;
	case 1:
		iDataSize = sizeof(short int);
		break;
	case 3:
		iDataSize = sizeof(float);
		break;
	case 2:
	case 4:
		iDataSize = sizeof(int);
		break;
	default:
		iDataSize = sizeof(double);
	}

	if (!(pImageData = (void*)calloc(((unsigned long)iRowMax) * iColMax, iDataSize)))
	{
		printf("\nError in ReadInputImage: insufficient memory.\n");
		scanf_s(" %s", &cChar, 1);
		exit(1);
	}

	fread((char*)pImageData, iDataSize, ((unsigned long)iRowMax) * iColMax, fin);

	fclose(fin);

	return (void*)(pImageData);
}

/***
Read image data.
Rreturn NULL if cannot open file
Parameters
iDataType:	= 0 - unsigned char, 1 - short int, 2 - int, 3 - float, 4 - int, 5 - double
v2 (11/2/2020): add parameter iOffset
***/
void* ReadInputImage_v2(char* pacImageFileName, int iColMax, int iRowMax, int iDataType, unsigned long ulOffsetElementsNum)
{
	char cChar;
	int iDataSize;
	void* pImageData = NULL;
	FILE* fin = NULL;

	fopen_s(&fin, pacImageFileName, "rb");
	if (fin == NULL)
	{
		printf("\nError in ReadInputImage(): cannot open %s.\n", pacImageFileName);
		scanf_s(" %c", &cChar, 1);
		exit(1);
		return NULL;
	}

	switch (iDataType)
	{
	case 0:
		iDataSize = sizeof(unsigned char);
		break;
	case 1:
		iDataSize = sizeof(short int);
		break;
	case 3:
		iDataSize = sizeof(float);
		break;
	case 2:
	case 4:
		iDataSize = sizeof(int);
		break;
	default:
		iDataSize = sizeof(double);
	}

	if (!(pImageData = (void*)calloc(((unsigned long)iRowMax) * iColMax, iDataSize)))
	{
		printf("\nError in ReadInputImage: insufficient memory.\n");
		scanf_s(" %s", &cChar, 1);
		exit(1);
	}

	// v2
	_fseeki64(fin, ulOffsetElementsNum * iDataSize, SEEK_SET);

	fread((char*)pImageData, iDataSize, ((unsigned long)iRowMax) * iColMax, fin);

	fclose(fin);

	return (void*)(pImageData);
}

bool FindTargetValueInWindow_UC(unsigned char* pImg, int iWidth, int iHeight, int iTargetCol, int iTargetRow, int w, unsigned char targetValue)
{
	int Row, Col;

	for (Row = iTargetRow - w; Row <= iTargetRow + w; Row++)
	{
		for (Col = iTargetCol - w; Col <= iTargetCol + w; Col++)
		{
			if (Row < 0 || Row >= iHeight || Col < 0 || Col >= iWidth)
				continue;

			if (pImg[Row * iWidth + Col] == targetValue)
				return true;
		}
	}

	return false;
}

// get number of lines in a txt file
int Gettxtlines(char* file)
{
	FILE *fp = NULL;
	int i;
	char cLine[2000];
	
	fopen_s(&fp, file, "r");
	if (fp == NULL)
	{
		printf("%s not exist\n", file);
		scanf_s(" %d", &i);
		exit(1);
	}

	i = 0;
	while (fgets(cLine, 2000, fp) != NULL)
		i += 1;

	fclose(fp);

	return i;
}

// open a file for reading
FILE* Readtxt(char* file)
{
	FILE *fp;
	int i;
	errno_t	err;
	if ((err = fopen_s(&fp, file, "r")) != 0)
	{
		printf_s("%s not exist\n", file);
		scanf_s(" %d", &i);
		exit(1);
	}
	return fp;
}

// open a file for writing
FILE* Writetxt(char* file)
{
	FILE *fp;
	errno_t	err;
	int i;
	if ((err = fopen_s(&fp, file, "w")) != 0)
	{
		printf_s("%s not exist\n", file);
		scanf_s(" %d", &i);
		exit(1);
	}
	return fp;
}

// open a file for writing
FILE* WriteBinary(char* file)
{
	FILE *fp;
	int i;
	fopen_s(&fp, file, "wb");
	if (fp == NULL)
	{
		printf("%s cannot be openned for writing\n", file);
		scanf_s(" %d", &i);
		exit(1);
	}
	return fp;
}

bool FileExist(char* file)
{
	FILE *fp;
	fopen_s(&fp, file, "r");
	if (fp != NULL)
	{
		// file exists
		fclose(fp);
		return true;
	}
	else
		return false;
}

/***
v1 (8/26/2019): print error message and pause if file not exist
***/
bool CheckFileExist(char* file)
{
	FILE* fp;
	int i;
	fopen_s(&fp, file, "r");
	if (fp != NULL)
	{
		// file exists
		fclose(fp);
		return true;
	}
	else
	{
		printf("File not exist: %s\n", file);
		scanf_s(" %d", &i);
		exit(1);
		return false;
	}
}

// Output Envi HDR header file
void OutputEnviHDRFile(char *pacPathName, char *pacDescription, int iSamples, int iLines, int iBands, int iDataType, char *pacBandNames)
{
	char pacPathName_hdr[STRLEN] = "";

	sprintf_s(pacPathName_hdr, STRLEN, "%s.hdr", pacPathName);

	FILE *fout = NULL;
	fout = Writetxt(pacPathName_hdr);
	fprintf(fout, "ENVI\n");
	fprintf(fout, "description = {\n  %s}\n", pacDescription);
	fprintf(fout, "samples\t= %d\n", iSamples);
	fprintf(fout, "lines\t= %d\n", iLines);
	fprintf(fout, "bands   = %d\n", iBands);
	fprintf(fout, "header offset = 0\n");
	fprintf(fout, "file type = ENVI Standard\n");
	fprintf(fout, "data type = %d\n", iDataType);
	fprintf(fout, "interleave = bsq\nsensor type = Unknown\nbyte order = 0\nwavelength units = Unknown\n");
	if (strlen(pacBandNames) > 0)
		fprintf(fout, "band names = {\n %s}", pacBandNames);

	fclose(fout);

	return;
}

bool InputEnviHDRFile(char* pacPathName, int *ptr_iSamples, int *ptr_iLines, int *ptr_iBands, int *ptr_iDataType)
{
	FILE* fin = NULL;
	char cLine[STRLEN] = "";
	char* pacPos = NULL;
	int i;

	// check if file exists
	if (FileExist(pacPathName) == false)
	{
		printf("%s does not exist.\n", pacPathName);
		scanf_s(" %d", &i);
		exit(1);
		return false;
	}

	fopen_s(&fin, pacPathName, "r");
	// read registration file
	while (fgets(cLine, STRLEN, fin) != NULL)
	{
		if (strstr(cLine, "samples	=") != NULL)
			sscanf_s(cLine + strlen("samples	="), "%d", ptr_iSamples);
		else if (strstr(cLine, "lines	=") != NULL)
			sscanf_s(cLine + strlen("lines	="), "%d", ptr_iLines);
		else if (strstr(cLine, "bands   =") != NULL)
			sscanf_s(cLine + strlen("bands   ="), "%d", ptr_iBands);
		else if (strstr(cLine, "data type =") != NULL)
			sscanf_s(cLine + strlen("data type ="), "%d", ptr_iDataType);
		else
			continue;
	}
	fclose(fin);

	return true;
}

bool GetFileDir(char *pacFilePathName, char *pacFileDir)
{
	char * pacPos = strrchr(pacFilePathName, '/');

	if (pacPos == NULL)
		return false;

	pacFileDir[0] = '\0';
	strncpy_s(pacFileDir, STRLEN, pacFilePathName, strlen(pacFilePathName) - strlen(pacPos));
	pacFileDir[strlen(pacFilePathName) - strlen(pacPos)] = '\0';

	return true;
}

bool GetUpperDir(char *pacFilePathName, char *pacFileDir)
{
	char * pacPos = strrchr(pacFilePathName, '/');
	char pacFilePathName_cpy[2000];

	if (pacPos == NULL)
		return false;

	if (strlen(pacPos) == 1)
	{
		// pacFilePathName ends with '/'
		strncpy_s(pacFilePathName_cpy, STRLEN, pacFilePathName, strlen(pacFilePathName) - 1);
		pacFilePathName_cpy[strlen(pacFilePathName) - 1] = '\0';
		pacFileDir[0] = '\0';

		return GetFileDir(pacFilePathName_cpy, pacFileDir);
	}
	else
		return GetFileDir(pacFilePathName, pacFileDir);

	return true;
}

/***
Get file name from a full path name.
Note there should be at least one '/' in pacFilePathName
***/
bool GetFileName(char *pacFilePathName, char *pacFileName)
{
	char * pacPos = strrchr(pacFilePathName, '/');

	pacFileName[0] = '\0';

	if (pacPos == NULL)
		return false;

	strncpy_s(pacFileName, STRLEN, pacPos + 1, strlen(pacPos) - 1);
	pacFileName[strlen(pacPos)-1] = '\0';

	return true;
}

// Get file name from a full path name and remove the extension
bool GetFileName_NoExtension(char *pacFilePathName, char *pacFileName)
{
	char * pacPos = strrchr(pacFilePathName, '/');

	pacFileName[0] = '\0';

	if (pacPos == NULL)
		return false;

	strncpy_s(pacFileName, STRLEN, pacPos + 1, strlen(pacPos) - 1);
	pacFileName[strlen(pacPos) - 1] = '\0';

	// remove extension
	pacPos = strrchr(pacFileName, '.');
	if (pacPos != NULL)
		pacFileName[strlen(pacFileName) - strlen(pacPos)] = '\0';

	return true;
}

float GetStd_UC(unsigned char* pucData, int n, unsigned char* pucMask)
{
	int i;
	float fSum, fMean, fStd;
	int n_valid;

	fSum = 0;
	n_valid = 0;
	for (i = 0; i < n; i++)
	{
		if (pucMask[i] > 0)
		{
			fSum += pucData[i];
			n_valid += 1;
		}
	}

	if (n_valid == 0)
		return -1;

	fMean = fSum / n_valid;
	fSum = 0;
	for (i = 0; i < n; i++)
	{
		if (pucMask[i] > 0)
			fSum += (pucData[i] - fMean) * (pucData[i] - fMean);
	}

	fStd = sqrt(fSum / n_valid);

	return fStd;
}

void GetMinMax(double* pData, int n, bool bExludeZero, double* pMin, double* pMax)
{
	int i;
	*pMax = -DBL_MAX;
	*pMin = DBL_MAX;
	for (i = 0; i < n; i++)
	{
		if (bExludeZero && ABS(pData[i]) < 1e-20)
			continue;

		if (pData[i] > * pMax)
			* pMax = pData[i];
		if (pData[i] < *pMin)
			* pMin = pData[i];
	}

	return;
}
