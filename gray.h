#ifndef _GRAY_H_
#define _GRAY_H_

#include "tistdtypes.h"

typedef struct _GRAY8
{
    Uint8 * gray;
	Int32 width;	//���ؿ��
	Int32 height;	//���ظ߶�
}GRAY8;

typedef struct _Point
{
	Uint16 x;	//width����
	Uint16 y;	//height����
}Point;

#endif



