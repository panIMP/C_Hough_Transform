#ifndef _ALG_INCISION_H
#define _ALG_INCISION_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


/*------------------------------------------------------------------------------
Standard types
------------------------------------------------------------------------------*/
#ifndef Int8
#define Int8 char
#endif

#ifndef UInt8
#define UInt8 unsigned char
#endif

#ifndef Int32
#define Int32 int
#endif

#ifndef UInt32
#define UInt32 unsigned int
#endif

#ifndef BOOL
#define BOOL short
#endif

#ifndef False
#define False 0
#endif

#ifndef True
#define True 1
#endif


/*------------------------------------------------------------------------------
Algorithm types
------------------------------------------------------------------------------*/
// Descriping the circle found by Hough transform
typedef struct _CIRCLE_HP {
    // x,y cordinate start from the top-left border
    Int32 x;
    Int32 y;
    // radius
    Int32 r;
    // r2Div is a transformed value representing its radius length r.
    // r2Div = (r * r + div / 2) / div,
    // if div = pow(2, n), div2 stands for "div / 2",
    // then, r2Div = (r * r + div2) >> n,
    Int32 r2Div;
    // indicate whether the circle has been found
    BOOL  hasValue;
} circleHp;

// Descriping the parameter for Hough transform detecting circle
typedef struct _HOUGH_CIRCLE_PARAM_HP {
    // searching range for the center of the circle
    Int32 xStart;
    Int32 xEnd;
    Int32 yStart;
    Int32 yEnd;
    // searching range for the radius length of the circle
    Int32 rMin;
    Int32 rMax;
    // see explaination in _CIRCLE_HP
    Int32 div2;
    Int32 n;
    // calculate by every "searchStep" image pixels
    Int32 searchStep;
    // counter number threshold to be regarded as a circle
    Int32 countThresh;
} cleHParam;

// Descriping the line found by Hough transform
typedef struct _LINE_HP {
    // line is descriped as x*cos(sita) + y*sin(sita) = dist;
    Int32 dist;
    Int32 sita;
    // indicate whether the circle has been found
    BOOL  hasValue;
} lineHp;

// Descriping the parameter for Hough transform detecting line
typedef struct _HOUGH_LINE_PARAM_HP {
    // searching range for the pixel on the line
    Int32 iStart;
    Int32 iEnd;
    Int32 jStart;
    Int32 jEnd;
    // calculate by every "searchStep" image pixels
    Int32 searchStep;
} lineHParam;


/*------------------------------------------------------------------------------
Algorithm const variables
------------------------------------------------------------------------------*/
// Only search for big circle with radius at range of
// (BIG_CIRCLE_MIN, BIG_CIRCLE_MAX)
#define BIG_CIRCLE_R_MIN 94
#define BIG_CIRCLE_R_MAX 98

// Only search for small circle with radius at range of
// (SMALL_CIRCLE_MIN, SMALL_CIRCLE_MAX)
#define SMALL_CIRCLE_R_MIN 3
#define SMALL_CIRCLE_R_MAX 6

// "div2" & "n" value explained in _CIRCLE_HP
#define BIG_CIRCLE_DIV_2 128
#define BIG_CIRCLE_N 8    // pow(2, R_MOV_NUMBER) = 2 * DIV_2
#define SMALL_CIRCLE_DIV_2 0
#define SMALL_CIRCLE_N 0    // pow(2, R_MOV_NUMBER) = 2 * DIV_2

// Search step for finding big circle and small circle
#define BIG_SEARCH_STEP 1
#define SMALL_SEARCH_STEP 1
#define LINE_SEARCH_SETP 1

// Border width of the big circle
#define BIG_CIRCLE_WIDTH 4

// Counter number threshold to be regarded as a circle
#define BIG_CIRCLE_COUNT_THRESH   100
#define SMALL_CIRCLE_COUNT_THRESH 8

// H array sizes
#define H_X_SIZE 700
#define H_Y_SIZE 500
#define H_R2DIV_SIZE 40
#define H_YR2DIV_SIZE  20000
#define H_XYR2DIV_SIZE 14000000

#define H_DIST_SIZE 1000
#define H_SITA_SIZE 360
#define H_DISTSITA_SIZE 360000


/*------------------------------------------------------------------------------
Algorithm functions
------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------
Description: Sobel transform
             Input image data:
             yyyyyyyyyyyyyyyy
             yxxxxxxxxxxxxxxy
             yxxxxxxxxxxxxxxy
             yxxxxxxxxxxxxxxy
             yxxxxxxxxxxxxxxy
             yyyyyyyyyyyyyyyy

             Output image data:
             tXXXXXXXXXXXXXXz
             zXXXXXXXXXXXXXXz
             zXXXXXXXXXXXXXXz
             zXXXXXXXXXXXXXXt

             where,
             X = sobel(x)
             t = y (input image value, unchanged)
             z = sobel(y)

Input:       imageData -- aligned source image data
             width -- image width
             height -- image height

Output:      imageData -- filtered image

Returns:     NULL
----------------------------------------------------------------------------*/
void sobel(UInt8* imageData, Int32 width, Int32 height);


/*----------------------------------------------------------------------------
Description: Binary the image by mean statistic thresh in 8-neighbour pixels
             Input image data:
             yyyyyyyyyyyyyyyy
             yxxxxxxxxxxxxxxy
             yxxxxxxxxxxxxxxy
             yxxxxxxxxxxxxxxy
             yxxxxxxxxxxxxxxy
             yyyyyyyyyyyyyyyy

             Output image data:
             tXXXXXXXXXXXXXXz
             zXXXXXXXXXXXXXXz
             zXXXXXXXXXXXXXXz
             tXXXXXXXXXXXXXXz

             where,
             X = meanBinary(x)
             t = y (input image value, unchanged)
             z = meanBinary(y)

Input:       imageData -- aligned source image data
             width -- image width
             height -- image height

Output:      imageData -- filtered image

Returns:     NULL
----------------------------------------------------------------------------*/
void meanBinary (UInt8* imageData, Int32 width, Int32 height);


/*----------------------------------------------------------------------------
Description: Initialize the Hough parameters for finding the biggest circle

Input:       width -- image width
             height -- image height

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
void initHPForBigCle (Int32 width, Int32 height);


/*----------------------------------------------------------------------------
Description: Get the static big circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
cleHParam getBHParam ();


/*----------------------------------------------------------------------------
Description: Initialize the Hough parameters for finding the small circle
             inside the biggest circle

Input:       bigCle -- the biggest circle that has been found previously

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
void initHPForSmallCle (circleHp bigCle);


/*----------------------------------------------------------------------------
Description: Get the static big circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
cleHParam getSHParam ();


/*----------------------------------------------------------------------------
Description: Initialize the Hough parameters for finding the small circle
             inside the biggest circle

Input:       bigCle -- the biggest circle that has been found previously

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
void initHPForLine (circleHp bigCle);


/*----------------------------------------------------------------------------
Description: Get the static small circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
lineHParam getLHParam ();


/*----------------------------------------------------------------------------
Description: Clear the data in 3 dimention array

Input:       array -- array to be cleared
             d1, d2, d3 -- array[d1][d2][d3]

Output:      array

Returns:     NULL
----------------------------------------------------------------------------*/
void clearArray (void* array, Int32 size);


/*----------------------------------------------------------------------------
Description: Find the biggest circle in the image by prehanded Hough Transform
             parameters.

Input:       imageData -- unsigned char image gray value array
             width -- image width
             height -- image height
             hParam -- Hough Transform parameters

Output:      circleHp -- information of the biggest circle that may be found

Returns:     same as above
----------------------------------------------------------------------------*/
circleHp hTForBigCle (UInt8* imageData,
                      Int32 width, Int32 height,
                      cleHParam hParam);


/*----------------------------------------------------------------------------
Description: Find the small circle inside the big circle

Input:       imageData -- unsigned char image gray value array
             width -- image width
             height -- image height
             hParam -- Hough Transform parameters
             bigCle -- big circle information

Output:      circleHp -- information of the biggest circle that may be found

Returns:     same as above
----------------------------------------------------------------------------*/
circleHp hTForSmallCle (UInt8* imageData,
                        Int32 width, Int32 height,
                        cleHParam hParam,
                        circleHp bigCle);


/*----------------------------------------------------------------------------
Description: Find the line as the string of the big circle

Input:       imageData -- unsigned char image gray value array
             width -- image width
             height -- image height
             hParam -- Hough Transform parameters
             bigCle -- big circle information

Output:      lineHp -- information of the string

Returns:     same as above
----------------------------------------------------------------------------*/
lineHp hTForLine (UInt8* imageData,
                  Int32 width, Int32 height,
                  lineHParam hParam);




#endif // alg_incision.h

/*-------------------------------- End of file -------------------------------*/

