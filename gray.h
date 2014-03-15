#ifndef _GRAY_H_
#define _GRAY_H_

#include "tistdtypes.h"

typedef struct _GRAY8
{
    Uint8 * gray;
	Int32 width;	//像素宽度
	Int32 height;	//像素高度
}GRAY8;

typedef struct _Point
{
	Uint16 x;	//width方向
	Uint16 y;	//height方向
}Point;

#endif



