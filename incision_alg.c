#include "incision_alg.h"
#include "basecfg.h"

#include "showUtils.h"

// Global H array for storing Hough transform counts
#pragma DATA_ALIGN (H_CIRCLE_LEARN_ARRAY, 32);
static Uint16 H_CIRCLE_LEARN_ARRAY[H_X_SIZE][H_Y_SIZE][H_R_SIZE] = {0};
#pragma DATA_ALIGN (H_LINE_ARRAY, 32);
static Uint16 H_LINE_ARRAY[H_SITA_SIZE][H_DIST_SIZE] = {{0}};

// Tables
extern Int16 sinValue[H_SITA_SIZE];
extern Int16 cosValue[H_SITA_SIZE];

// Global Hough parameter for big & small circle
static circleHp  bigCleLearn;
static circleHp	 smallCleLearn;

static hParamCle hParamBigLearn;
static hParamCle hParamSmallLearn;
static hParamCle hParamBigDetect;
static hParamCle hParamSmallDetect;

void sobel (Uint8* restrict srcImgData, Uint8* restrict dstImgData, Int16 width, Int16 height) {
    Int32 i;

    Int16 H, O, V;
    Uint8  i00,  i01,  i02;
    Uint8  i10,        i12;
    Uint8  i20,  i21,  i22;

    Int32 tranSize = width * (height - 2) - 2;

    // Bond I/O 
	Uint8* restrict in  = srcImgData;
    Uint8* restrict out = dstImgData;

    // Iterate over entire image as a single, continuous raster line.
	#pragma MUST_ITERATE(307200, 360960, 2) // 640*480 ~ 752*480
	#pragma UNROLL(2);
    for (i = 0; i < tranSize; ++i)
    {
        // Read in the required 3x3 region from the input.
		i00=in[i        ]; 		i01=in[i    	+1]; 	i02=in[i    	+2];
        i10=in[i+  width];                  			i12=in[i+  width+2];
        i20=in[i+2*width]; 		i21=in[i+2*width+1]; 	i22=in[i+2*width+2];

        // Apply horizontal and vertical filter masks.  The final filter
        // output is the sum of the absolute values of these filters.
        H = -   i00 - 2*i01 -   i02 +
            +   i20 + 2*i21 +   i22;

        V = -   i00         +   i02
            - 2*i10         + 2*i12
            -   i20         +   i22;

        O = abs (H) + abs (V);

        // Clamp to 8-bit range.  The output is always positive due to
        // the absolute value, so we only need to check for overflow.
        if (O > 255) O = 255;

        // Store it.
        out[i + 1 + width] = O;
    }
}


void meanBinary (Uint8* restrict srcImgData, Uint8* restrict dstImgData, Int16 width, Int16 height) {
    Int32 i;

    Int16 sum;
    Uint8 i00, i01, i02;
    Uint8 i10, i11, i12;
    Uint8 i20, i21, i22;

    Int32 tranSize = width * (height - 2) - 2;

    // Bond I/O 
	Uint8* restrict in  = srcImgData;
    Uint8* restrict out = dstImgData;

    // Iterate over entire image as a single, continuous raster line.
	#pragma MUST_ITERATE(307200, 360960, 2) // 640*480 ~ 752*480
	#pragma UNROLL(2);
    for (i = 0; i < tranSize; ++i)
    {
        // Read in the required 3x3 region from the input.
		i00=in[i        ]; 		i01=in[i    	+1]; 	i02=in[i    	+2];
        i10=in[i+  width];      i11=in[i+width	+1];    i12=in[i+  width+2];
        i20=in[i+2*width]; 		i21=in[i+2*width+1]; 	i22=in[i+2*width+2];

        // Calculate the 8-neighbour pixel value sum, include itself
        sum = i00 + i01 + i02 +
              i10 + i11 + i12 +
              i20 + i21 + i22;

        // Clamp to 8-bit range.  The output is always positive due to
        // the absolute value, so we only need to check for overflow.
        if (sum < 5 * 255)
            sum = 0;
        else
            sum = 255;

        // Store it.
        out[i + 1 + width] = sum;
    }
}


void edgeExtract (Uint8* restrict imgData, Uint8* restrict tmpData, Int16 width, Int16 height, circleHp bigCle) {
	Int32 i,j;
	Int16 cleX = bigCle.x;
	Int16 cleY = bigCle.y;

	Bool edgeFind;
	Uint8* restrict imgPtr = imgData;
	Uint8* restrict tmpPtr = tmpData;
	
	Int32 fullSize = width * height;
	Int16 sum;

    memcpy (tmpPtr, imgPtr, sizeof(Uint8) * fullSize); 

	for (j = 0; j < height; ++j) {
		edgeFind = FALSE;
		imgPtr = imgData + j * width + cleX;
		for (i = cleX; i < width; ++i) {
			if (edgeFind == TRUE) {
				*imgPtr++ = 0;
				continue;
			}
			if (*imgPtr == 255) {
				edgeFind = TRUE;
			}
			imgPtr ++;
		}
	}
    
	for (j = 0; j < height; ++j) {
		edgeFind = FALSE;
		imgPtr = imgData + j * width + cleX;
		for (i = cleX; i > 0; --i) {
			if (edgeFind == TRUE) {
				*imgPtr-- = 0;
				continue;
			}
			if (*imgPtr == 255) {
				edgeFind = TRUE;
			}
			imgPtr --;
		}
	}	

	for (i = 0; i < width; ++i) {
		edgeFind = FALSE;
		tmpPtr = tmpData + cleY * width + i;
		for (j = cleY; j < height; ++j) {
			if (edgeFind == TRUE) {
				*tmpPtr = 0;
				tmpPtr += width;
				continue;
			}
			if (*tmpPtr == 255) {
				edgeFind = TRUE;
			}
			tmpPtr += width;
		}
	}
    
	for (i = 0; i < width; ++i) {
		edgeFind = FALSE;
		tmpPtr = tmpData + cleY * width + i;
		for (j = cleY; j > 0; --j) {
			if (edgeFind == TRUE) {
				*tmpPtr = 0;
				tmpPtr -= width;
				continue;
			}
			if (*tmpPtr == 255) {
				edgeFind = TRUE;
			}
			tmpPtr -= width;
		}
	}	
	
	for (i = 0; i < fullSize; ++i) {
		sum = imgData[i] + tmpData[i];
		if (sum > 255) {
			sum = 255;
		}
		imgData[i] = sum;	
	}						
}


void drawCross (Uint8* imgData, Int16 x, Int16 y, Uint8 len, Int16 width, Int16 height) {
	Int16 i,j;

	for (i = x - len; i < x + len; ++i) {
		*(imgData + y * width + i) = 255; 
	}

	for (j = y - len; j < y + len; ++j) {
		*(imgData + j * width + x) = 255; 
	}
}


void initHPForBigCleLearn (Int16 width, Int16 height) {
	Uint32 div = pow (2, BIG_CIRCLE_N);
	Uint32 rMax1 = (height < width ? height : width) / 2;
	Uint32 rMax2 = (Int16)sqrt((float)((H_R_SIZE << BIG_CIRCLE_N) - (div / 2)));
	hParamBigLearn.rMin = (Int16)(sqrt ((float)div));
	hParamBigLearn.rMax = rMax1 < rMax2 ? rMax1 : rMax2;
    hParamBigLearn.n = BIG_CIRCLE_N;
    hParamBigLearn.searchStep = BIG_LEARN_SEARCH_STEP; 

	hParamBigLearn.div2 = div / 2;
	hParamBigLearn.r2DivMax = (hParamBigLearn.rMax * hParamBigLearn.rMax + hParamBigLearn.div2) >> BIG_CIRCLE_N;
	hParamBigLearn.r2DivMin = (hParamBigLearn.rMin * hParamBigLearn.rMin + hParamBigLearn.div2) >> BIG_CIRCLE_N;

	hParamBigLearn.hArray = (Uint16*)H_CIRCLE_LEARN_ARRAY;
	hParamBigLearn.used = FALSE;
}


hParamCle* getHParamBigLearn () {
    return &hParamBigLearn;
}

circleHp* getBigCleLearn () {
	return &bigCleLearn;
}


int hTForBigCleLearn (Uint8* imgData, Int16 width, Int16 height, circleHp* bigCle, hParamCle* hParam) {
   // Iterators
    Int16 i,j,x,y;   // i -- x, j -- y

    // Localize those of the Hough parameter
    Int16 rMin = hParam->rMin;
    Int16 rMax = hParam->rMax;
    Int32 r2Min = rMin * rMin;
    Int32 r2Max = rMax * rMax;
    Int8 n = hParam->n;
	Int16 div2 = pow (2, n-1);
    Int32 r2;
    Int16 r2Div;
	Int16 r2DivMin = hParam->r2DivMin;
	Int16 r2DivMax = hParam->r2DivMax;
    Int8 step = hParam->searchStep;
	Uint16* hArray = hParam->hArray;

	// Ranges for circle center
    Int16 xStart = rMin;
    Int16 xEnd = width - rMin;
    Int16 yStart = rMin;
    Int16 yEnd = height - rMin;
    Int16 xStartTmp, xEndTmp, yStartTmp, yEndTmp;
    Int16 xMin, xMax, yMin, yMax;
    
	// Distance from circle center
    Int16 xDist, yDist;

	// image pointer
	Uint8* imgPtr;

    // Counts
    Int16 tmpNeighCount = 0;
    Int16 maxNeighCount = 0;
    Int16 maxCount = 0;

	// Clear H array
	memset (hArray, 0, sizeof (Uint16) * H_LEARN_SIZE);
	bigCle->hasValue = FALSE;
	bigCle->r = 0;
	bigCle->r2Div = 0;
	bigCle->x = 0;
	bigCle->y = 0;

	// Indicate the learning process has been conducted
	hParam->used = TRUE;

    // Start algorithm
    for (j = 0; j < height; j += step) {
		imgPtr = imgData + j * width;
        for (i = 0; i < width; ++i) {
            if (*imgPtr++ == 255) {
                // current cicle center range
                xStartTmp = i - rMax;
                xEndTmp = i + rMax;
                yStartTmp = j - rMax;
                yEndTmp = j + rMax;

                // compare with the predefined range to choose the shortest one
                xMin = xStartTmp > xStart ? xStartTmp : xStart;
                xMax = xEndTmp < xEnd ? xEndTmp : xEnd;
                yMin = yStartTmp > yStart ? yStartTmp : yStart;
                yMax = yEndTmp < yEnd ? yEndTmp : yEnd;

                // count only those within the r range
                for (x = xMin; x < xMax; ++x) {
                    for (y = yMin; y < yMax; ++y) {
                        xDist = i - x;
                        yDist = j - y;
                        r2 = xDist * xDist + yDist * yDist;
                        if (r2 > r2Min && r2 < r2Max) {
                            r2Div = (r2 + div2) >> n;
                        	H_CIRCLE_LEARN_ARRAY[x][y][r2Div] ++;
                        }
                    }
                }
            }
        }
    }

    // Generally, the biggest count in hCount stands for the most significant
    // circle, here, for improvement of accuracy, add up the counts of the
    // 8-neighborhoods of the center of circle, finding the biggest.
    for (x = xStart; x < xEnd; ++x) {
        for (y = yStart; y < yEnd; ++y) {
            for (r2Div = r2DivMin; r2Div < r2DivMax; ++r2Div) {
                tmpNeighCount = H_CIRCLE_LEARN_ARRAY[x-1][y-1][r2Div] +
                                H_CIRCLE_LEARN_ARRAY[x-1][y][r2Div]   +
                                H_CIRCLE_LEARN_ARRAY[x-1][y+1][r2Div] +
                                H_CIRCLE_LEARN_ARRAY[x][y-1][r2Div]   +
                                H_CIRCLE_LEARN_ARRAY[x][y][r2Div]     +
                                H_CIRCLE_LEARN_ARRAY[x][y+1][r2Div]   +
                                H_CIRCLE_LEARN_ARRAY[x+1][y-1][r2Div] +
                                H_CIRCLE_LEARN_ARRAY[x+1][y][r2Div]   +
                                H_CIRCLE_LEARN_ARRAY[x+1][y+1][r2Div];

                if (tmpNeighCount > maxNeighCount) {
                    maxNeighCount = tmpNeighCount;
                    maxCount = H_CIRCLE_LEARN_ARRAY[x][y][r2Div];
                    bigCle->x = x;
                    bigCle->y = y;
                    bigCle->r2Div = r2Div;
                }
            }
        }
    }

	if (maxNeighCount < 54) {
		bigCle->hasValue = FALSE;
		bigCle->r = 0;
		bigCle->r2Div = 0;
		bigCle->x = 0;
		bigCle->y = 0;
		printf ("No chip is found!\n");
	} else {
    	bigCle->hasValue = TRUE;
    	bigCle->r = (Int16)sqrtf ((float)((bigCle->r2Div << n) - div2));
	    printf ("Step is %d\n", step);
	    printf ("Count number is %d\n", maxCount);
    	printf ("Big circle (%d, %d),radius:%d\n", bigCle->x, bigCle->y, bigCle->r);
	}
    printf ("-----------------------------------------------------\n\n");

    return 0;	
}


void initHPForSmallCleLearn () {
	Uint32 div = pow (2, SMALL_CIRCLE_N);
	hParamSmallLearn.rMin = 1;
	hParamSmallLearn.rMax = (Int16)sqrt((float)((H_R_SIZE << SMALL_CIRCLE_N) - (div / 2)));
	hParamSmallLearn.n = SMALL_CIRCLE_N;
	hParamSmallLearn.searchStep = SMALL_DETECT_LEARN_STEP;
	hParamSmallLearn.div2 = div / 2;
	hParamSmallLearn.r2DivMax = (hParamSmallLearn.rMax * hParamSmallLearn.rMax + hParamSmallLearn.div2) >> SMALL_CIRCLE_N;
	hParamSmallLearn.r2DivMin = (hParamSmallLearn.rMin * hParamSmallLearn.rMin + hParamSmallLearn.div2) >> SMALL_CIRCLE_N;
	hParamSmallLearn.hArray = (Uint16*) H_CIRCLE_LEARN_ARRAY;
	hParamSmallLearn.used = FALSE;
}


hParamCle* getHParamSmallLearn () {
    return &hParamSmallLearn;
}


circleHp* getSmallCleLearn () {
	return &smallCleLearn;
}


int hTForSmallCleLearn (Uint8* imgData, Int16 width, Int16 height, circleHp* smallCle, hParamCle* hParam, circleHp bigCle) {
    // Iterators
    Int32 i,j,x,y;   // i -- x, j -- y

    // Localize those of the Hough parameter
    Int8 rMin = hParam->rMin;
    Int8 rMax = hParam->rMax;
    Int16 r2Min = rMin * rMin;
    Int16 r2Max = rMax * rMax;
    Int8 n = hParam->n;
    Int8 step = hParam->searchStep;
    Int16 div2 = hParam->div2;
    Int8 r2DivMin = hParam->r2DivMin;
    Int8 r2DivMax = hParam->r2DivMax;
	Uint16* hArray = hParam->hArray;
    Int32 r2;
    Int8 r2Div;

    // Range of image pixel fetching
    Int16 iMin = bigCle.x - bigCle.r;
    Int16 iMax = bigCle.x + bigCle.r;
    Int16 jMin = bigCle.y - bigCle.r;
    Int16 jMax = bigCle.y + bigCle.r;
	Uint8* imgPtr;

	// Range of small circle center
    Int16 xMin,xMax,yMin,yMax;
    Int16 xStart = iMin + rMin;
    Int16 xEnd = iMax - rMin;
    Int16 yStart = jMin + rMin;
    Int16 yEnd = jMax - rMin;
    Int16 xStartTmp,xEndTmp,yStartTmp,yEndTmp;
    Int16 xDist, yDist;

    // Thresh for pixel distance from big circle center, make sure to
    // exclude the border & near border pixels
    Int32 r2ThreForInsideBigCle = bigCle.r * bigCle.r;

    // Counts
    Int16 tmpNeighCount = 0;
    Int16 maxNeighCount = 0;
    Int16 maxCount = 0;

	// Clear H array
	memset (hArray, 0, sizeof (Uint16) * H_LEARN_SIZE);
	smallCle->hasValue = FALSE;
	smallCle->r = 0;
	smallCle->r2Div = 0;
	smallCle->x = 0;
	smallCle->y = 0;

	// Indicate the learning process has been conducted
	hParam->used = TRUE;

    // Start algorithm
    for (j = jMin; j < jMax; j += step) {
		imgPtr = imgData + j * width + iMin;
        for (i = iMin; i < iMax; ++i) {
            if (*imgPtr++ == 255) {
                // exclude the border & near border pixels
                xDist = i - bigCle.x;
                yDist = j - bigCle.y;
                r2 = xDist * xDist + yDist * yDist;
                if (r2 >= r2ThreForInsideBigCle) {
                    continue;
                }

                // current cicle center range
                xStartTmp = i - rMax;
                xEndTmp = i + rMax;
                yStartTmp = j - rMax;
                yEndTmp = j + rMax;

                // compare with the predefined range to choose the shortest one
                xMin = xStartTmp > xStart ? xStartTmp : xStart;
                xMax = xEndTmp < xEnd ? xEndTmp : xEnd;
                yMin = yStartTmp > yStart ? yStartTmp : yStart;
                yMax = yEndTmp < yEnd ? yEndTmp : yEnd;

                // count only those within the r range
                for (x = xMin; x < xMax; ++x) {
                    for (y = yMin; y < yMax; ++y) {
                        xDist = i - x;
                        yDist = j - y;
                        r2 = xDist * xDist + yDist * yDist;
                        if (r2 > r2Min && r2 < r2Max) {
                            r2Div = (r2 + div2) >> n;
                        	H_CIRCLE_LEARN_ARRAY[x][y][r2Div] ++;
                        }
                    }
                }
            }
        }
    }

    // Generally, the biggest count in hCount stands for the most significant
    // circle, here, for improvement of accuracy, add up the counts of the
    // 8-neighborhoods of the center of circle, finding the biggest.
    for (x = xStart; x < xEnd; ++x) {
        for (y = yStart; y < yEnd; ++y) {
            for (r2Div = r2DivMin; r2Div < r2DivMax; ++r2Div) {
                tmpNeighCount = H_CIRCLE_LEARN_ARRAY[x-1][y-1][r2Div] +
                                H_CIRCLE_LEARN_ARRAY[x-1][y][r2Div]   +
                                H_CIRCLE_LEARN_ARRAY[x-1][y+1][r2Div] +
                                H_CIRCLE_LEARN_ARRAY[x][y-1][r2Div]   +
                                H_CIRCLE_LEARN_ARRAY[x][y][r2Div]     +
                                H_CIRCLE_LEARN_ARRAY[x][y+1][r2Div]   +
                                H_CIRCLE_LEARN_ARRAY[x+1][y-1][r2Div] +
                                H_CIRCLE_LEARN_ARRAY[x+1][y][r2Div]   +
                                H_CIRCLE_LEARN_ARRAY[x+1][y+1][r2Div];
                if (tmpNeighCount > maxNeighCount) {
                    maxNeighCount = tmpNeighCount;
                    maxCount = H_CIRCLE_LEARN_ARRAY[x][y][r2Div];
                    smallCle->x = x;
                    smallCle->y = y;
                    smallCle->r2Div = r2Div;
                }
            }
        }
    }

    // if circle((x,y) as center, r2Div for r) has more than COUNT_THRESH,
    // then it can be regarded as an explicit circle.
    if (maxNeighCount < 54) {
		smallCle->hasValue = FALSE;
		smallCle->r = 0;
		smallCle->r2Div = 0;
		smallCle->x = 0;
		smallCle->y = 0;
    	printf ("It's the back side, please place a front side chip for learning process!\n");   
    } else {
        smallCle->hasValue = TRUE;
        smallCle->r = (Int16)sqrtf ((float)((smallCle->r2Div << n) - div2));
        printf ("Step is %d\n", step);
        printf ("Count number is %d\n", maxCount);
        printf ("Small circle (%d,%d),radius:%d\n", smallCle->x, smallCle->y, smallCle->r); 		
	}
    printf ("-----------------------------------------------------\n\n");

    return 0;
}


void initHPForBigCleDetect (Int16 width, Int16 height, hParamCle hParamCleLearn) {
	hParamBigDetect = hParamCleLearn;

	// overlap these
	hParamBigDetect.searchStep = BIG_DETECT_SEARCH_STEP;
	if (hParamBigDetect.hArray == NULL) {
		hParamBigDetect.hArray = MEM_calloc (DDR2, sizeof (Uint16) * width * height, 32);
		if (hParamBigDetect.hArray == MEM_ILLEGAL) {
			printf ("Memory allocation failed in init hough param for big circle detect!\n");
			exit(-1);
		}
	}
}


hParamCle getHParamBigDetect () {
	return hParamBigDetect;
}

circleHp hTForBigCleDetect (Uint8* imgData, Int16 width, Int16 height, circleHp bigCle, hParamCle hParam) {
    // Iterators
    Int16 i,j,x,y;   // i -- x, j -- y

    // Localize those of the Hough parameter
    Int16 div2 = hParam.div2;
    Int8 n = hParam.n;
    Int8 step = hParam.searchStep;
	Uint16* hArray = hParam.hArray;
	Int16 r = bigCle.r;

	// Search range
    Int16 xMin,xMax,yMin,yMax;
    Int16 xStart = bigCle.r;
    Int16 xEnd = width - bigCle.r;
    Int16 yStart = bigCle.r;
    Int16 yEnd = height - bigCle.r;
    Int16 xStartTmp,xEndTmp,yStartTmp,yEndTmp;

	// Image pointer
	Uint8* imgPtr;
    Int16 xDist, yDist;
    Int8 r2Div;

    // Counts
    Int16 tmpNeighCount = 0;
    Int16 maxNeighCount = 0;
    Int16 maxCount = 0;

    // Initialize circle variable
    circleHp cle = {0,0,0,0,FALSE};
	cle.r = bigCle.r;

    // Initialize H_CIRCLE_LEARN_ARRAY
	memset (hArray, 0, sizeof(Uint16) / sizeof(Uint8) * width * height);

    // Start algorithm
    for (j = 0; j < height; j += step) {
		imgPtr = imgData + j * width;
        for (i = 0; i < width; ++i) {
            if (*imgPtr++ == 255) {
                // current cicle center range
                xStartTmp = i - r;
                xEndTmp = i + r;
                yStartTmp = j - r;
                yEndTmp = j + r;

                // compare with the predefined range to choose the shortest one
                xMin = xStartTmp > xStart ? xStartTmp : xStart;
                xMax = xEndTmp < xEnd ? xEndTmp : xEnd;
                yMin = yStartTmp > yStart ? yStartTmp : yStart;
                yMax = yEndTmp < yEnd ? yEndTmp : yEnd;

                // count only those within the r range
                for (x = xMin; x < xMax; ++x) {
                    for (y = yMin; y < yMax; ++y) {
                        xDist = i - x;
                        yDist = j - y;
                        r2Div = (xDist * xDist + yDist * yDist + div2) >> n;
                        if (r2Div == bigCle.r2Div) {
                            hArray[y * width + x] ++;
                        }
                    }
                }
            }
        }
    }

	//hArray = hParam.hArray;

    // Generally, the biggest count in hCount stands for the most significant
    // circle, here, for improvement of accuracy, add up the counts of the
    // 8-neighborhoods of the center of circle, finding the biggest.
    for (x = xStart; x < xEnd; ++x) {
        for (y = yStart; y < yEnd; ++y) {
            tmpNeighCount = hArray[(y-1)*width + (x-1)] +
                            hArray[(y  )*width + (x-1)] +
                            hArray[(y+1)*width + (x-1)] +
                            hArray[(y-1)*width + (x  )] +
                            hArray[(y  )*width + (x  )] +
                            hArray[(y+1)*width + (x  )] +
                            hArray[(y-1)*width + (x+1)] +
                            hArray[(y  )*width + (x+1)] +
                            hArray[(y+1)*width + (x+1)] ;

            if (tmpNeighCount > maxNeighCount) {
                maxNeighCount = tmpNeighCount;
                maxCount = *(hArray + y * width + x);
                cle.x = x;
                cle.y = y;
            }
        }
    }
	
	if (maxNeighCount < 54) {
		cle.hasValue = FALSE;
		printf ("No chip is found!\n");
	} else {
	    cle.hasValue = TRUE;
	    printf ("Step is %d\n", step);
	    printf ("Count number is %d\n", maxCount);
	    printf ("Big circle (%d, %d),radius:%d\n", cle.x, cle.y, cle.r);
	}
    printf ("-----------------------------------------------------\n\n");

    return cle;
}


void initHPForSmallCleDetect (Int16 width, Int16 height, hParamCle hParamCleLearn) {
	hParamSmallDetect = hParamCleLearn;

	// overlap these
	hParamSmallDetect.searchStep = SMALL_DETECT_SEARCH_STEP;
	if (hParamSmallDetect.hArray == NULL) {
		hParamSmallDetect.hArray = MEM_calloc (DDR2, sizeof (Uint16) * width * height, 32);
		if (hParamSmallDetect.hArray == MEM_ILLEGAL) {
			printf ("Memory allocation failed in init hough param for small circle detect!\n");
			exit(-1);
		}
	}
}


hParamCle getHParamSmallDetect () {
    return hParamSmallDetect;
}


circleHp hTForSmallCleDetect (Uint8* imgData, Int16 width, Int16 height, circleHp smallCle, hParamCle hParam, circleHp bigCle) {
    // Iterators
    Int32 i,j,x,y;   // i -- x, j -- y

    // Localize those of the Hough parameter
    Int8 n = hParam.n;
    Int8 step = hParam.searchStep;
    Int16 div2 = hParam.div2;
    Int16 r = smallCle.r;
    Int32 r2;
    Int8 r2Div;
	Uint16* hArray = hParam.hArray;

    // Range of image pixel fetching
    Int16 iMin = bigCle.x - bigCle.r;
    Int16 iMax = bigCle.x + bigCle.r;
    Int16 jMin = bigCle.y - bigCle.r;
    Int16 jMax = bigCle.y + bigCle.r;
	Uint8* imgPtr;

	// Range of small circle center
    Int16 xMin,xMax,yMin,yMax;
    Int16 xStart = iMin + r;
    Int16 xEnd = iMax - r;
    Int16 yStart = jMin + r;
    Int16 yEnd = jMax - r;
    Int16 xStartTmp,xEndTmp,yStartTmp,yEndTmp;
    Int16 xDist, yDist;

    // Thresh for pixel distance from big circle center, make sure to
    // exclude the border & near border pixels
    Int32 r2ThreForInsideBigCle = bigCle.r * bigCle.r;

    // Counts
    Int16 tmpNeighCount = 0;
    Int16 maxNeighCount = 0;
    Int16 maxCount = 0;

    // Initialize circle variable
    circleHp cle = {0,0,0,0,FALSE};
	cle.r = smallCle.r;
	cle.r2Div = smallCle.r2Div;

	// Clear H array
	memset (hArray, 0, sizeof (Uint16) * width * height);

    // Start algorithm
    for (j = jMin; j < jMax; j += step) {
		imgPtr = imgData + j * width + iMin;
        for (i = iMin; i < iMax; ++i) {
            if (*imgPtr++ == 255) {
                // exclude the border & near border pixels
                xDist = i - bigCle.x;
                yDist = j - bigCle.y;
                r2 = xDist * xDist + yDist * yDist;
                if (r2 >= r2ThreForInsideBigCle) {
                    continue;
                }

                // current cicle center range
                xStartTmp = i - r;
                xEndTmp = i + r;
                yStartTmp = j - r;
                yEndTmp = j + r;

                // compare with the predefined range to choose the shortest one
                xMin = xStartTmp > xStart ? xStartTmp : xStart;
                xMax = xEndTmp < xEnd ? xEndTmp : xEnd;
                yMin = yStartTmp > yStart ? yStartTmp : yStart;
                yMax = yEndTmp < yEnd ? yEndTmp : yEnd;

                // count only those within the r range
                for (x = xMin; x < xMax; ++x) {
                    for (y = yMin; y < yMax; ++y) {
                        xDist = i - x;
                        yDist = j - y;
                        r2Div = (xDist * xDist + yDist * yDist + div2) >> n;
                        if (r2Div == smallCle.r2Div) {
                            hArray[y * width + x] ++;
                        }
                    }
                }
            }
        }
    }

    // Generally, the biggest count in hCount stands for the most significant
    // circle, here, for improvement of accuracy, add up the counts of the
    // 8-neighborhoods of the center of circle, finding the biggest.
    for (x = xStart; x < xEnd; ++x) {
        for (y = yStart; y < yEnd; ++y) {
            tmpNeighCount = hArray[(y-1)*width + (x-1)] +
                            hArray[(y  )*width + (x-1)] +
                            hArray[(y+1)*width + (x-1)] +
                            hArray[(y-1)*width + (x  )] +
                            hArray[(y  )*width + (x  )] +
                            hArray[(y+1)*width + (x  )] +
                            hArray[(y-1)*width + (x+1)] +
                            hArray[(y  )*width + (x+1)] +
                            hArray[(y+1)*width + (x+1)] ;

            if (tmpNeighCount > maxNeighCount) {
                maxNeighCount = tmpNeighCount;
                maxCount = *(hArray + y * width + x);
                cle.x = x;
                cle.y = y;
            }
        }
    }

    // if circle((x,y) as center, r2Div for r) has more than COUNT_THRESH,
    // then it can be regarded as an explicit circle.
    if (maxNeighCount < 54) {
		cle.hasValue = FALSE;
    	printf ("It's the back side!\n");   
    } else {
        cle.hasValue = TRUE;
        cle.r = (Int16)sqrtf ((float)((cle.r2Div << n) - div2));
        printf ("Step is %d\n", step);
        printf ("Count number is %d\n", maxCount);
        printf ("Small circle (%d,%d),radius:%d\n", cle.x, cle.y, cle.r); 		
	}
    printf ("-----------------------------------------------------\n\n");

    return cle;
}


lineHp hTForLineDetect (Uint8* imgData, Int16 width, Int16 height, circleHp bigCle) {
    Int16 i,j;
	Int32 k,kNum;
    Int16 sita,dist;
	Uint16 count;
    Uint16 maxCount = 0;

	Int16 iStart = bigCle.x - bigCle.r;
	Int16 iEnd = bigCle.x + bigCle.r;
	Int16 jStart = bigCle.y - bigCle.r;
	Int16 jEnd = bigCle.y + bigCle.r; 
	Uint8* imgPtr;
    Int16* hLineArray = (Int16*) H_LINE_ARRAY;

    lineHp line = {0, 0, FALSE};

    // Initialize H_LINE_ARRAY
	memset (H_LINE_ARRAY, 0, sizeof(Uint16) / sizeof(Uint8) * H_DISTSITA_SIZE);

    // Hough transform
	for (j = jStart; j < jEnd; ++j) {
		imgPtr = imgData + j * width + iStart;
		for (i = iStart; i < iEnd; ++i) {
			if (*imgPtr++ == 255) {		
				for (sita = 0; sita < H_SITA_SIZE; ++sita) {
					dist = (Int16)((Int32)((Int32)i * (Int32)cosValue[sita] + (Int32)j * (Int32)sinValue[sita]) >> 10);
					
					if (dist > 0 && dist < H_DIST_SIZE) {
						H_LINE_ARRAY[sita][dist]++;
					}
				}
			}	
		}
	}

    // Find the line
    for (k = 0; k < H_DISTSITA_SIZE; ++k) {
        count = *hLineArray++;
        if (count > maxCount) {
            maxCount = count;
            kNum = k;
        }
    }

	//if (maxCount > LINE_COUNT_THRESH) {
    	line.dist = kNum % H_DIST_SIZE;
    	line.sita = kNum / H_DIST_SIZE;	
		line.hasValue = TRUE;
	//}

	// Draw the line
    dist = line.dist;
    sita = line.sita;
    for (i = 0; i < width; ++i) {
        for (j = 0; j < height; ++j) {
            if (((i * cosValue[sita] + j * sinValue[sita]) >> 10) == dist) {
                *(imgData + j*width + i) = 255;
            }
        }
    }


	/*imgPtr = imgData;
	for (j = 0; j < height; ++j) {
		for (i = 0; i < width; ++i) {
            if (((i * cosValue[sita] + j * sinValue[sita]) >> 10) == dist) {
                *imgPtr = 255;
            }			
			imgPtr++;
		}
		imgPtr += width;
	}*/

	//x = (  B*B*x0  -  A*B*y0  -  A*C  ) / ( A*A + B*B );
	//y  =  ( -A*B*x0 + A*A*y0 - B*C  ) / ( A*A + B*B );
	/*lineXEnd = ((Int32)((Int32)sinValue[sita] * (Int32)sinValue[sita] * lineXStart) >> 20)  - 
			   ((Int32)((Int32)cosValue[sita] * (Int32)sinValue[sita] * lineYStart) >> 20)  + 
			   ((Int32)cosValue[sita] * (Int32)dist >> 10);
	lineYEnd = ((-(Int32)cosValue[sita] * (Int32)sinValue[sita] * lineXStart) >> 20) + 
			   (((Int32)cosValue[sita] * (Int32)cosValue[sita] * lineYStart) >> 20)  +
			   ((Int32)sinValue[sita] * (Int32)dist >> 10);
    lineDist = ((Int32)sinValue[sita] * lineXStart - (Int32)cosValue[sita] * lineYStart) >> 10;
	
	lineXMin = lineXStart < lineXEnd ? lineXStart : lineXEnd;
	lineXMax = lineXStart > lineXEnd ? lineXStart : lineXEnd;
	lineYMin = lineYStart < lineYEnd ? lineYStart : lineYEnd;
	lineYMax = lineYStart > lineYEnd ? lineYStart : lineYEnd;

    for (i = lineXMin; i < lineXMax; ++i) {
        for (j = lineYMin; j < lineYMax; ++j) {
            if ((( i * (Int32)sinValue[sita] - j * (Int32)cosValue[sita]) >> 10) == lineDist) {
                *(imgData + j*width + i) = 255;
            }
        }
    }*/

    printf ("Count number is %d\n", maxCount);
    printf ("Line angle is %d, r is %d\n", line.sita, line.dist);
    printf ("-----------------------------------------------------\n\n");

    return line;
}


int incisionLearn (Uint8* imgData, Int16 width, Int16 height) {
	// Initialization
	char strBuf[40];
	Point commentPos;
	GRAY8 imageGray;
	Int32 fullSize = width * height;

	Uint8* tmpData = MEM_calloc (DDR2, sizeof (Uint8) * fullSize, 32);
	if (tmpData == MEM_ILLEGAL) {
		printf ("Memory allocation failed in tmp image in incision detect!\n");
		exit(-1);
	}

	commentPos.x = 0;
	commentPos.y = 0;
	imageGray.gray = imgData;
	imageGray.width = width;
	imageGray.height = height;


    // Sharpen the image
    sobel (imgData, tmpData, width, height);


    // Binary the image, extract the chip profile
    meanBinary (tmpData, imgData, width, height);
 

    // Find the chip circle
    initHPForBigCleLearn (width, height);
    hTForBigCleLearn (imgData, width, height, getBigCleLearn(), getHParamBigLearn ());

  
    // Find the spot circle
    initHPForSmallCleLearn ();
    hTForSmallCleLearn (imgData, width, height, getSmallCleLearn(), getHParamSmallLearn (), *getBigCleLearn ());


	// Mark detect result
	if (getBigCleLearn()->hasValue) {
		sprintf (strBuf, "Chip Center: (%d, %d)", getBigCleLearn()->x, getBigCleLearn()->y); 
		showString (&imageGray, &commentPos, strBuf);
		commentPos.x = 0;
		commentPos.y += 16;
		sprintf (strBuf, "Chip Radius: %d", getBigCleLearn()->r);
		showString (&imageGray, &commentPos, strBuf);
		commentPos.x = 0;
		commentPos.y += 16;

		drawCross (imgData, getBigCleLearn()->x, getBigCleLearn()->y, 4, width, height);
	} else {
		sprintf (strBuf, "No chip is found!"); 
		showString (&imageGray, &commentPos, strBuf);
    	MEM_free (DDR2, tmpData, sizeof (Uint8*) * fullSize);
		return 0;
	}

	if (getSmallCleLearn()->hasValue) {
		sprintf (strBuf, "Spot Radius: %d", getSmallCleLearn()->r); 
		showString (&imageGray, &commentPos, strBuf);	
		commentPos.x = 0;
		commentPos.y += 16;

		drawCross (imgData, getSmallCleLearn()->x, getSmallCleLearn()->y, 4, width, height);		
	} else {
		sprintf (strBuf, "It's not front side!"); 
		showString (&imageGray, &commentPos, strBuf);		
	}
    
    MEM_free (DDR2, tmpData, sizeof (Uint8*) * fullSize);

	return 0;
}


int incisionDetect (Uint8* imgData, Int16 width, Int16 height) {
    // Initialization
    circleHp bCle, sCle;
	lineHp line;

	Int32 fullSize = width * height;
	char strBuf[40];
	Point commentPos;
	GRAY8 imageGray;

	static Uint8* tmpData = NULL;

	commentPos.x = 0;
	commentPos.y = 0;
	imageGray.gray = imgData;
	imageGray.width = width;
	imageGray.height = height;


	// Judge if the learning process has been conducted!
	if (getHParamBigLearn()->used == FALSE) {
		sprintf (strBuf, "Learning process undone!");
		showString (&imageGray, &commentPos, strBuf);
		return -1;		
	}

	if (getBigCleLearn()-> hasValue == FALSE) {			 
		sprintf (strBuf, "No chip in learning proc!");
		showString (&imageGray, &commentPos, strBuf);				
		return -1;		
	}

	/*if (getSmallCleLearn()-> hasValue == FALSE) {
		sprintf (strBuf, "Back side in learning proc!");
		showString (&imageGray, &commentPos, strBuf);	
		return -1;			
	}*/


	// If has learned, allocate memory for detecting process
	if (tmpData == NULL) {
		tmpData = MEM_calloc (DDR2, sizeof(Uint8) * fullSize, 32);
		if (tmpData == MEM_ILLEGAL) {
			printf ("Memory allocation failed in tmp image in incision detect!");
			exit(-1);
		}
	}

    // Sharpen the image
    sobel (imgData, tmpData, width, height);


    // Binary the image, extract the chip profile
    meanBinary (tmpData, imgData, width, height);
 

    // Find the chip circle
    initHPForBigCleDetect (width, height, *getHParamBigLearn());
    bCle = hTForBigCleDetect (imgData, width, height, *getBigCleLearn (), getHParamBigDetect ());


    // Find the spot circle
    initHPForSmallCleDetect (width, height, *getHParamSmallLearn());
    sCle = hTForSmallCleDetect (imgData, width, height, *getSmallCleLearn() ,getHParamSmallDetect (), bCle);


    // Find the string line
	edgeExtract (imgData, tmpData, width, height, bCle);
    line = hTForLineDetect (imgData, width, height, bCle);


	// Mark detect result
	if (bCle.hasValue) {
		sprintf (strBuf, "Chip Center: (%d, %d)", bCle.x, bCle.y); 
		showString (&imageGray, &commentPos, strBuf);
		commentPos.x = 0;
		commentPos.y += 16;

		drawCross (imgData, bCle.x, bCle.y, 4, width, height);
	} else {
		sprintf (strBuf, "No chip is found!"); 
		showString (&imageGray, &commentPos, strBuf);
		return 0;
	}

	if (sCle.hasValue) {
		sprintf (strBuf, "Chip Side: Front"); 
		showString (&imageGray, &commentPos, strBuf);	
		commentPos.x = 0;
		commentPos.y += 16;		

		drawCross (imgData, sCle.x, sCle.y, 4, width, height);		
	} else {
		sprintf (strBuf, "Chip Side: Back"); 
		showString (&imageGray, &commentPos, strBuf);	
		commentPos.x = 0;
		commentPos.y += 16;				
	}

	if (line.hasValue) {
		sprintf (strBuf, "Incision Angle: %d", line.sita); 
		showString (&imageGray, &commentPos, strBuf);
	}

    //MEM_free (DDR2, tmpData, sizeof (Uint8*) * fullSize);

	return 0;
}


/*-------------------------------- End of file -------------------------------*/

