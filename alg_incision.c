#include "alg_incision.h"

// Global H array for storing Hough transform counts
static Int32 H_CIRCLE_ARRAY[H_X_SIZE][H_Y_SIZE][H_R2DIV_SIZE] = {{{0}}};
static Int32 H_LINE_ARRAY[H_SITA_SIZE][H_DIST_SIZE] = {{0}};

// Tables
extern const Int32 sinValue[H_SITA_SIZE];
extern const Int32 cosValue[H_SITA_SIZE];

// Global Hough parameter for big & small circle
static cleHParam  bHParam;
static cleHParam  sHParam;
static lineHParam lHParam;


void sobel (UInt8* imageData, Int32 width, Int32 height) {
    Int32 i;

    Int32 H, O, V;
    UInt8  i00,  i01,  i02;
    UInt8  i10,        i12;
    UInt8  i20,  i21,  i22;
    UInt8* pI00; UInt8* pI01; UInt8* pI02;
    UInt8* pI10;              UInt8* pI12;
    UInt8* pI20; UInt8* pI21; UInt8* pI22;

    Int32 fullSize = width * height;
    Int32 tranSize = width * (height - 2) - 2;

    // Bond output
    UInt8* out = imageData;

    // Copy init data
    UInt8* in = (UInt8*) malloc (fullSize * sizeof(UInt8));
    if (in == NULL) {
        printf ("Melmory allocation failed in sobel!");
        exit(1);
    }
    for (i = 0; i < fullSize; ++i) {
        *in++ = *imageData++;
    }
    in = in - fullSize;

    // Init 8-neighbour pointers
    pI00 = in;              pI01 = in+1;            pI02 = in+2;
    pI10 = in+width;                                pI12 = in+2+width;
    pI20 = in+2*width;      pI21 = in+1+2*width;    pI22 = in+2+2*width;

    // Iterate over entire image as a single, continuous raster line.
    for (i = 0; i < tranSize; ++i)
    {
        // Read in the required 3x3 region from the input.
        i00 = *pI00++;    i01 = *pI01++;    i02 = *pI02++;
        i10 = *pI10++;                      i12 = *pI12++;
        i20 = *pI20++;    i21 = *pI21++;    i22 = *pI22++;


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

    // free intermediate copy image
    free (in);
}


void meanBinary (UInt8* imageData, Int32 width, Int32 height) {
    Int32 i;

    Int32 sum;
    UInt8 i00, i01, i02;
    UInt8 i10, i11, i12;
    UInt8 i20, i21, i22;
    UInt8* pI00; UInt8* pI01; UInt8* pI02;
    UInt8* pI10; UInt8* pI11; UInt8* pI12;
    UInt8* pI20; UInt8* pI21; UInt8* pI22;

    Int32 fullSize = width * height;
    Int32 tranSize = width * (height - 2) - 2;

    // Bond output
    UInt8* out = imageData;

    // Copy init data
    UInt8* in = malloc (width * height * sizeof(UInt8));
    if (in == NULL) {
        printf ("Melmory allocation failed in meanBinary!");
        exit(1);
    }
    for (i = 0; i < fullSize; ++i) {
        *in++ = *imageData++;
        if (i < width || i > (tranSize + width) || i % width > 650) {
            *imageData = 0;
        }
    }
    in = in - fullSize;

    // Init 8-neighbour pointers
    pI00 = in;              pI01 = in+1;            pI02 = in+2;
    pI10 = in+width;        pI11 = in+1+width;      pI12 = in+2+width;
    pI20 = in+2*width;      pI21 = in+1+2*width;    pI22 = in+2+2*width;

    // Iterate over entire image as a single, continuous raster line.
    for (i = 0; i < tranSize; ++i)
    {
        // Read in the required 3x3 region from the input.
        i00 = *pI00++;    i01 = *pI01++;    i02 = *pI02++;
        i10 = *pI10++;    i11 = *pI11++;    i12 = *pI12++;
        i20 = *pI20++;    i21 = *pI21++;    i22 = *pI22++;

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

    // free intermediate copy image
    free (in);
}


void initHPForBigCle (Int32 width, Int32 height){
    bHParam.xStart = BIG_CIRCLE_R_MIN;
    bHParam.xEnd = width - BIG_CIRCLE_R_MIN;
    bHParam.yStart = BIG_CIRCLE_R_MIN;
    bHParam.yEnd = height - BIG_CIRCLE_R_MIN;
    bHParam.rMin = BIG_CIRCLE_R_MIN;
    bHParam.rMax = BIG_CIRCLE_R_MAX;
    bHParam.div2 = BIG_CIRCLE_DIV_2;
    bHParam.n = BIG_CIRCLE_N;
    bHParam.searchStep = BIG_SEARCH_STEP;
    bHParam.countThresh = BIG_CIRCLE_COUNT_THRESH;
}


cleHParam getBHParam () {
    return bHParam;
}


void initHPForSmallCle (circleHp bigCle) {
    if (bigCle.hasValue) {
        sHParam.xStart = bigCle.x - bigCle.r + SMALL_CIRCLE_R_MIN;
        sHParam.xEnd = bigCle.x + bigCle.r - SMALL_CIRCLE_R_MIN;
        sHParam.yStart = bigCle.y - bigCle.r + SMALL_CIRCLE_R_MIN;
        sHParam.yEnd = bigCle.y + bigCle.r - SMALL_CIRCLE_R_MIN;
        sHParam.rMin = SMALL_CIRCLE_R_MIN;
        sHParam.rMax = SMALL_CIRCLE_R_MAX;
        sHParam.div2 = SMALL_CIRCLE_DIV_2;
        sHParam.n = SMALL_CIRCLE_N;
        sHParam.searchStep = SMALL_SEARCH_STEP;
        sHParam.countThresh = SMALL_CIRCLE_COUNT_THRESH;
    } else {
        printf ("Biggest circle hasn't been found yet!");
        exit(1);
    }
}


cleHParam getSHParam () {
    return sHParam;
}


void initHPForLine (circleHp bigCle) {
    lHParam.iStart      = bigCle.x - bigCle.r;
    lHParam.iEnd        = bigCle.x + bigCle.r;
    lHParam.jStart      = bigCle.y - bigCle.r;
    lHParam.jEnd        = bigCle.y + bigCle.r;
    lHParam.searchStep  = LINE_SEARCH_SETP;
}


lineHParam getLHParam () {
    return lHParam;
}


void clearArray (void* array, Int32 size) {
    int i;

    for (i = 0; i < size; ++i) {
        *((Int32*)array + i) = 0;
    }
}


circleHp hTForBigCle(UInt8* imageData,
                     Int32 width, Int32 height,
                     cleHParam hParam) {
    // Iterators
    Int32 i,j,x,y;   // i -- x, j -- y

    // Localize those of the Hough parameter
    Int32 xMin,xMax,yMin,yMax;
    Int32 xStart = xMin = hParam.xStart;
    Int32 xEnd = xMax = hParam.xEnd;
    Int32 yStart = yMin = hParam.yStart;
    Int32 yEnd = yMax = hParam.yEnd;
    Int32 xStartTmp,xEndTmp,yStartTmp,yEndTmp;
    Int32 xDist, yDist;
    Int32 rMin = hParam.rMin;
    Int32 rMax = hParam.rMax;
    Int32 r2Min = rMin * rMin;
    Int32 r2Max = rMax * rMax;
    Int32 div2 = hParam.div2;
    Int32 n = hParam.n;
    Int32 r2DivMin = (r2Min + div2) >> n;
    Int32 r2DivMax = (r2Max + div2) >> n;
    Int32 r2,r2Div;
    Int32 step = hParam.searchStep;
    Int32 countThresh = hParam.countThresh;

    // Range of image pixel fetching
    Int32 iMin = xStart - rMin;
    Int32 iMax = xEnd + rMin;
    Int32 jMin = yStart - rMin;
    Int32 jMax = yEnd + rMin;

    // Counts
    Int32 tmpNeighCount = 0;
    Int32 maxNeighCount = 0;
    Int32 maxCount = 0;

    // Initialize circle variable
    circleHp cle = {0,0,0,0,False};

    // Initialize H_CIRCLE_ARRAY
    clearArray (H_CIRCLE_ARRAY, H_XYR2DIV_SIZE);

    // Start algorithm
    for (i = iMin; i < iMax; i += step) {
        for (j = jMin; j < jMax; ++j) {
            if (*(imageData + j * width + i) == 255) {
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
                            H_CIRCLE_ARRAY[x][y][r2Div] ++;
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
                tmpNeighCount = H_CIRCLE_ARRAY[x-1][y-1][r2Div] +
                                H_CIRCLE_ARRAY[x-1][y][r2Div]   +
                                H_CIRCLE_ARRAY[x-1][y+1][r2Div] +
                                H_CIRCLE_ARRAY[x][y-1][r2Div]   +
                                H_CIRCLE_ARRAY[x][y][r2Div]     +
                                H_CIRCLE_ARRAY[x][y+1][r2Div]   +
                                H_CIRCLE_ARRAY[x+1][y-1][r2Div] +
                                H_CIRCLE_ARRAY[x+1][y][r2Div]   +
                                H_CIRCLE_ARRAY[x+1][y+1][r2Div];

                if (tmpNeighCount > maxNeighCount) {
                    maxNeighCount = tmpNeighCount;
                    maxCount = H_CIRCLE_ARRAY[x][y][r2Div];
                    cle.x = x;
                    cle.y = y;
                    cle.r2Div = r2Div;
                }
            }
        }
    }

    // if circle((x,y) as center, r2Div for r) has more than COUNT_THRESH,
    // then it can be regarded as an explicit circle.
    if (maxNeighCount > countThresh) {
        cle.r = (Int32)sqrtf ((float)((cle.r2Div << n) - div2));
        cle.hasValue = True;
        printf ("Step is %d\n", step);
        printf ("Count number is %d\n", maxCount);
    }

    return cle;
}


circleHp hTForSmallCle(UInt8* imgData,
                       Int32 width, Int32 height,
                       cleHParam hParam,
                       circleHp bigCle) {
    // Iterators
    Int32 i,j,x,y;   // i -- x, j -- y

    // Localize those of the Hough parameter
    Int32 xMin,xMax,yMin,yMax;
    Int32 xStart = xMin = hParam.xStart;
    Int32 xEnd = xMax = hParam.xEnd;
    Int32 yStart = yMin = hParam.yStart;
    Int32 yEnd = yMax = hParam.yEnd;
    Int32 xStartTmp,xEndTmp,yStartTmp,yEndTmp;
    Int32 xDist, yDist;
    Int32 rMin = hParam.rMin;
    Int32 rMax = hParam.rMax;
    Int32 r2Min = rMin * rMin;
    Int32 r2Max = rMax * rMax;
    Int32 div2 = hParam.div2;
    Int32 n = hParam.n;
    Int32 r2DivMin = (r2Min + div2) >> n;
    Int32 r2DivMax = (r2Max + div2) >> n;
    Int32 r2,r2Div;
    Int32 step = hParam.searchStep;
    Int32 countThresh = hParam.countThresh;

    // Thresh for pixel distance from big circle center, make sure to
    // exclude the border & near border pixels
    Int32 r2ThreForInsideBigCle = (bigCle.r - BIG_CIRCLE_WIDTH) *
                                  (bigCle.r - BIG_CIRCLE_WIDTH);

    // Range of image pixel fetching
    Int32 iMin = xStart - rMin;
    Int32 iMax = xEnd + rMin;
    Int32 jMin = yStart - rMin;
    Int32 jMax = yEnd + rMin;

    // Counts
    Int32 tmpNeighCount = 0;
    Int32 maxNeighCount = 0;
    Int32 maxCount = 0;

    // Initialize circle variable
    circleHp cle = {0,0,0,0,False};

    // Initialize H_CIRCLE_ARRAY
    clearArray (H_CIRCLE_ARRAY, H_XYR2DIV_SIZE);

    // Start algorithm
    for (i = iMin; i < iMax; i += step) {
        for (j = jMin; j < jMax; ++j) {
            if (*(imgData + j * width + i)) {
                // exclude the border & near border pixels
                xDist = i - bigCle.x;
                yDist = j - bigCle.y;
                r2 = xDist * xDist + yDist * yDist;
                if (r2 > r2ThreForInsideBigCle) {
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
                            H_CIRCLE_ARRAY[x][y][r2Div] ++;
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
                tmpNeighCount = H_CIRCLE_ARRAY[x-1][y-1][r2Div] +
                                H_CIRCLE_ARRAY[x-1][y][r2Div]   +
                                H_CIRCLE_ARRAY[x-1][y+1][r2Div] +
                                H_CIRCLE_ARRAY[x][y-1][r2Div]   +
                                H_CIRCLE_ARRAY[x][y][r2Div]     +
                                H_CIRCLE_ARRAY[x][y+1][r2Div]   +
                                H_CIRCLE_ARRAY[x+1][y-1][r2Div] +
                                H_CIRCLE_ARRAY[x+1][y][r2Div]   +
                                H_CIRCLE_ARRAY[x+1][y+1][r2Div];
                if (tmpNeighCount > maxNeighCount) {
                    maxNeighCount = tmpNeighCount;
                    maxCount = H_CIRCLE_ARRAY[x][y][r2Div];
                    cle.x = x;
                    cle.y = y;
                    cle.r2Div = r2Div;
                }
            }
        }
    }

    // if circle((x,y) as center, r2Div for r) has more than COUNT_THRESH,
    // then it can be regarded as an explicit circle.
    if (maxNeighCount > countThresh) {
        cle.r = (Int32)sqrtf ((float)((cle.r2Div << n) - div2));
        cle.hasValue = True;
        printf ("Step is %d\n", step);
        printf ("Count number is %d\n", maxCount);
    }

    return cle;
}


lineHp hTForLine (UInt8* imgData, Int32 width, Int32 height, lineHParam hParam) {
    Int32 i,j,x,y;
    Int32 sita,dist,count,num;
    Int32 maxCount = 0;

    Int32 fullSize = width * height;
    Int32* hLineArray = (Int32*) H_LINE_ARRAY;
    Int32 step = hParam.searchStep;
    UInt8* head = imgData;

    lineHp line = {0, 0, False};

    // Initialize H_LINE_ARRAY
    clearArray (H_LINE_ARRAY, H_DISTSITA_SIZE);

    // Hough transform
    for (i = 0; i < fullSize; i+= step) {
        if (!*head++) continue;

        x = i % width;  y = i / width;

        for (sita = 0; sita < H_SITA_SIZE; ++sita) {
            dist = (x * cosValue[sita] + y * sinValue[sita]) >> 10;

            if (dist > 0 && dist < H_DIST_SIZE)
                H_LINE_ARRAY[sita][dist]++;
        }
    }
    head = imgData;


    // Find the line
    for (i = 0; i < H_DISTSITA_SIZE; ++i) {
        count = *hLineArray++;
        if (count > maxCount) {
            maxCount = count;
            num = i;
        }
    }
    line.dist = num % H_DIST_SIZE;
    line.sita = num / H_DIST_SIZE;
    line.hasValue = True;

    dist = line.dist;
    sita = line.sita;
    for (i = 0; i < width; ++i) {
        for (j = 0; j < height; ++j) {
            if (((i*cosValue[sita] + j*sinValue[sita]) >> 10) == dist) {
                *(head + j*width + i) = 255;
            }
        }
    }

    printf ("Step is %d\n", step);
    printf ("Count number is %d\n", maxCount);

    return line;
}




/*-------------------------------- End of file -------------------------------*/

