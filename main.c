#include "alg_incision.h"
#include <opencv2/highgui/highgui_c.h>

#include <time.h>

int main(void)
{
    IplImage* image;
    circleHp bigCle, smallCle;
    clock_t timeStart, timeEnd;

    image = cvLoadImage ("G:/operation/f1.bmp", CV_LOAD_IMAGE_GRAYSCALE);

    // Time start
    timeStart = clock ();

    sobel ((UInt8*)image->imageData, image->width, image->height);

    meanBinary ((UInt8*)image->imageData, image->width, image->height);

    initHoughParamForBigCle (image->width, image->height);

    // Find the big circle
    bigCle = houghTransFormForBigCle ((UInt8*)image->imageData,
                                      image->width, image->height,
                                      getBHParam ());

    // If big circle exits, find the small circle inside the big circle
    if (bigCle.hasValue) {
        // output
        printf ("The big circle center is   (%d, %d), with radius of %d\n",
                bigCle.x, bigCle.y, bigCle.r);

        initHoughParamForSmallCle (bigCle);

        smallCle = houghTransFormForSmallCle ((UInt8*)image->imageData,
                                              image->width, image->height,
                                              getSHParam (),
                                              bigCle);

        // if small circle exits, meaning it's front side
        if (smallCle.hasValue) {
            printf ("The small circle center is (%d, %d), with radius of %d\n",
                    smallCle.x, smallCle.y, smallCle.r);
        } else {
            printf ("It's the back side, turn it upside down!\n");
        }
    } else {
        printf ("No circle is found!\n");
    }

    // Time end
    timeEnd = clock ();
    printf ("Hough Time consuming is %d ms", timeEnd - timeStart);

    return 0;
}

