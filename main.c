#include "alg_incision.h"


#include <opencv2/highgui/highgui_c.h>
#include <time.h>



int main(void)
{
    IplImage* image;
    UInt8* imageData;
    Int32 width, height;
    circleHp bCle, sCle;
    lineHp line;
    clock_t timeStart, timeEnd;

    // Initialize
    image = cvLoadImage ("G:/operation/f1.bmp", CV_LOAD_IMAGE_GRAYSCALE);
    imageData = (UInt8*)image->imageData;
    width = image->width;
    height = image->height;

    // Time start
    timeStart = clock ();

    // Sharpen the image
    sobel (imageData, width, height);
    timeEnd = clock ();
    printf ("Soble time consuming: %d\n", (Int32)timeEnd - (Int32)timeStart);
    timeStart = timeEnd;
    printf ("-----------------------------------------------------\n\n");

    // Binary the image, extract the chip profile
    meanBinary (imageData, width, height);
    timeEnd = clock ();
    printf ("Mean time consuming: %d\n", (Int32)timeEnd - (Int32)timeStart);
    timeStart = timeEnd;
    printf ("-----------------------------------------------------\n\n");

    // Find the chip circle
    initHPForBigCle (width, height);
    bCle = hTForBigCle (imageData, width, height, getBHParam ());
    if (!bCle.hasValue) {
        printf ("No circle is found!\n");
        return 0;
    }
    printf ("Big circle (%d, %d),radius:%d\n", bCle.x, bCle.y, bCle.r);
    timeEnd = clock ();
    printf ("Big circle time consuming: %d\n", (Int32)timeEnd - (Int32)timeStart);
    timeStart = timeEnd;
    printf ("-----------------------------------------------------\n\n");

    // Find the spot circle
    initHPForSmallCle (bCle);
    sCle = hTForSmallCle (imageData, width, height, getSHParam (), bCle);
    if (sCle.hasValue)
        printf ("Small circle (%d,%d),radius:%d\n", sCle.x, sCle.y, sCle.r);
    else
        printf ("It's the back side, turn it upside down!\n");
    timeEnd = clock ();
    printf ("Small circle time consuming: %d\n", (Int32)timeEnd - (Int32)timeStart);
    timeStart = timeEnd;
    printf ("-----------------------------------------------------\n\n");

    cvNamedWindow ("source", CV_WINDOW_AUTOSIZE);
    cvShowImage ("source", image);
    cvWaitKey (0);

    // Find the string line
    initHPForLine (bCle);
    line = hTForLine (imageData, width, height, getLHParam ());
    printf ("Line angle is %d, r is %d\n", line.sita, line.dist);
    timeEnd = clock ();
    printf ("Line time consuming: %d\n", (Int32)timeEnd - (Int32)timeStart);
    timeStart = timeEnd;
    printf ("-----------------------------------------------------\n\n");

    cvNamedWindow ("result", CV_WINDOW_AUTOSIZE);
    cvShowImage ("result", image);
    cvWaitKey (0);

    return 0;
}


/*-------------------------------- End of file -------------------------------*/


