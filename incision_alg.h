#ifndef _INCISION_ALG_H
#define _INCISION_ALG_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "tistdtypes.h"

/*=======================================================================================*/
/*=================================Total Process=========================================*/
/*=======================================================================================*/
// "n" value explained in _CIRCLE_LEARN
#define BIG_CIRCLE_N						10    	
// "n" value explained in _CIRCLE_LEARN
#define SMALL_CIRCLE_N						1 

#define H_X_SIZE 							800
#define H_Y_SIZE							600
#define H_R_SIZE							60
#define H_YR_SIZE							36000
#define H_LEARN_SIZE 						28800000

// Descriping the circle found by Hough transform learn process
typedef struct _CIRCLE_HP {
    // x,y cordinate start from the top-left border
    Int16 x;
    Int16 y;
    // radius
    Int16 r;
    // r2Div is a transformed value representing its radius length r.
    // r2Div = (r * r + div / 2) / div,
    // if div = pow(2, n), div2 stands for "div / 2",
    // then, r2Div = (r * r + div2) >> n,
    Int16 r2Div;
    // indicate whether the circle has been found
    Bool hasValue;
} circleHp;

// Descriping the parameter for Hough transform detecting circle
typedef struct _HOUGH_PARAM_CIRCLE {
    // searching range for the radius length of the circle
    Int16 rMin;
    Int16 rMax;
    // see explaination in _CIRCLE_HP
    Int8 n;
    // calculate by every "searchStep" image pixels
    Int8 searchStep;
	// H array
	Uint16* hArray; 
	// r2Div max value -- see _CIRCLE_LEARN
	Int16 r2DivMax;
	Int16 r2DivMin;
	Int16 div2;
	// indicate if the param has been used
	Bool used;
} hParamCle;

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

Input:       srcImgData -- aligned source image data
			 dstImgData -- aligned output image data
             width -- image width
             height -- image height

Output:      dstImgData -- aligned output image data

Returns:     NULL
----------------------------------------------------------------------------*/
void sobel(Uint8* restrict srcImgData, Uint8* restrict dstImgData, Int16 width, Int16 height);


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

Input:       srcImgData -- aligned source image data
			 dstImgData -- aligned output image data
             width -- image width
             height -- image height

Output:      dstImgData -- aligned output image data

Returns:     NULL
----------------------------------------------------------------------------*/
void meanBinary (Uint8* restrict srcImgData, Uint8* restrict dstImgData, Int16 width, Int16 height);


/*----------------------------------------------------------------------------
Description: Extract the edge of the circle

Input:       imgData -- aligned image data
			 tmpData -- another buf to store temporary image
			 width -- image width
			 height -- image height
             bigCle -- big circle that has been found

Output:      imgData -- binary image with only edge of the circle as foreground 

Returns:     NULL
----------------------------------------------------------------------------*/
void edgeExtract (Uint8* restrict imgData, Uint8* restrict tmpData, Int16 width, Int16 height, circleHp bigCle);


/*----------------------------------------------------------------------------
Description: Draw cross at certain pos 

Input:       imgData -- aligned image data
			 x -- pos.x
			 y -- pos.y
			 len -- cross length from center
			 width -- image width
			 height -- image height

Output:      imgData -- binary image with only edge of the circle as foreground 

Returns:     NULL
----------------------------------------------------------------------------*/
void drawCross (Uint8* imgData, Int16 x, Int16 y, Uint8 len, Int16 width, Int16 height);



/*=======================================================================================*/
/*=================================Learn Process=========================================*/
/*=======================================================================================*/

// Search step for finding chip circle and spot circle
#define BIG_LEARN_SEARCH_STEP				4
#define SMALL_DETECT_LEARN_STEP				1


/*----------------------------------------------------------------------------
Description: Initialize the Hough parameters for finding the biggest circle

Input:       width -- image width
             height -- image height

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
void initHPForBigCleLearn (Int16 width, Int16 height);


/*----------------------------------------------------------------------------
Description: Get the static chip circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     hough param for searching chip circle in learning process
----------------------------------------------------------------------------*/
hParamCle* getHParamBigLearn ();


/*----------------------------------------------------------------------------
Description: Get the static chip circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     chip circle info in learning process
----------------------------------------------------------------------------*/
circleHp* getBigCleLearn ();


/*----------------------------------------------------------------------------
Description: Find the chip circle in learing process

Input:       imgData -- unsigned char image gray value array
             width -- image width
             height -- image height
             hParam -- Hough Transform parameters

Output:      bigCle -- chip circle info
			 hParam -- modified hParam

Returns:     NULL
----------------------------------------------------------------------------*/
int hTForBigCleLearn (Uint8* imgData, Int16 width, Int16 height, circleHp* bigCle, hParamCle* hParam);


/*----------------------------------------------------------------------------
Description: Initialize the Hough parameters for finding the spot circle

Input:       NULL

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
void initHPForSmallCleLearn ();


/*----------------------------------------------------------------------------
Description: Get the static small circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     hough param for searching big circle in learning process
----------------------------------------------------------------------------*/
hParamCle* getHParamSmallLearn ();


/*----------------------------------------------------------------------------
Description: Get the static big circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     spot circle info in learning process
----------------------------------------------------------------------------*/
circleHp* getSmallCleLearn ();


/*----------------------------------------------------------------------------
Description: Find the spot circle in learing process

Input:       imgData -- unsigned char image gray value array
             width -- image width
             height -- image height
             hParam -- Hough Transform parameters
			 bigCle -- chip circle info

Output:      smallCle -- spot circle info
			 hParam -- modified hParam

Returns:     NULL
----------------------------------------------------------------------------*/
int hTForSmallCleLearn (Uint8* imgData, Int16 width, Int16 height, circleHp* smallCle, hParamCle* hParam, circleHp bigCle);


/*----------------------------------------------------------------------------
Description: main function of incision learning algorithm

Input:       imgData -- unsigned char image gray value array
             width -- image width
             height -- image height

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
int incisionLearn (Uint8* imgData, Int16 width, Int16 height); 


/*=======================================================================================*/
/*=================================Detect Process========================================*/
/*=======================================================================================*/

#define BIG_DETECT_SEARCH_STEP				8
#define SMALL_DETECT_SEARCH_STEP 			1
  

// Only search for small circle with radius at range of (SMALL_CIRCLE_MIN, SMALL_CIRCLE_MAX)
//#define SMALL_CIRCLE_R_MIN 1
//#define SMALL_CIRCLE_R_MAX 10

#define H_DIST_SIZE 						1000
#define H_SITA_SIZE 						360
#define H_DISTSITA_SIZE 					360000

// Descriping the line found by Hough transform
typedef struct _LINE_HP {
    // line is descriped as x*cos(sita) + y*sin(sita) = dist;
    Int16 dist;
    Int16 sita;
    // indicate whether the circle has been found
    Bool  hasValue;
} lineHp;



/*----------------------------------------------------------------------------
Description: Initialize the Hough parameters for finding the biggest circle

Input:       width -- image width
             height -- image height
			 hParamCleLearn -- correspond hough param in learning process

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
void initHPForBigCleDetect (Int16 width, Int16 height, hParamCle hParamCleLearn);


/*----------------------------------------------------------------------------
Description: Get the static big circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
hParamCle getHParamBigDetect ();


/*----------------------------------------------------------------------------
Description: Find the big circle length

Input:       imgData -- unsigned char image gray value array
             width -- image width
             height -- image height
			 bigCle -- chip circle info get from learning process
             hParam -- Hough Transform parameters

Output:      circleHp -- information of the biggest circle that may be found

Returns:     circleHp
----------------------------------------------------------------------------*/
circleHp hTForBigCleDetect (Uint8* imgData, Int16 width, Int16 height, circleHp bigCle, hParamCle hParam);


/*----------------------------------------------------------------------------
Description: Initialize the Hough parameters for finding the small circle
             inside the biggest circle

Input:       width -- image width
			 height -- image height
			 hParamCleLearn --correspond hough param in learning process

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
void initHPForSmallCleDetect (Int16 width, Int16 height, hParamCle hParamCleLearn);


/*----------------------------------------------------------------------------
Description: Get the static big circle Hough parameter

Input:       NULL

Output:      NULL

Returns:     hParamCle -- hough param for finding spot circle in detecting process
----------------------------------------------------------------------------*/
hParamCle getHParamSmallDetect ();


/*----------------------------------------------------------------------------
Description: Find the small circle inside the big circle

Input:       imgData -- unsigned char image gray value array
             width -- image width
             height -- image height
			 smallCle -- spot circle found in learning process
             hParam -- Hough Transform parameters for finding spot circle
             bigCle -- chip circle information found previously

Output:      circleHp -- information of the spot circle that may be found

Returns:     same as above
----------------------------------------------------------------------------*/
circleHp hTForSmallCleDetect (Uint8* imgData, Int16 width, Int16 height, circleHp smallCle, hParamCle hParam, circleHp bigCle);


/*----------------------------------------------------------------------------
Description: Find the line as the string of the big circle

Input:       imgData -- unsigned char image gray value array
             width -- image width
             height -- image height
             bigCle -- big circle information

Output:      lineHp -- information of the string

Returns:     lineHp
----------------------------------------------------------------------------*/
lineHp hTForLineDetect (Uint8* imgData, Int16 width, Int16 height, circleHp bigCle);


/*----------------------------------------------------------------------------
Description: main function of incision detection algorithm

Input:       imgData -- unsigned char image gray value array
             width -- image width
             height -- image height

Output:      NULL

Returns:     NULL
----------------------------------------------------------------------------*/
int incisionDetect (Uint8* imgData, Int16 width, Int16 height);


#endif // alg_incision.h

/*-------------------------------- End of file -------------------------------*/















