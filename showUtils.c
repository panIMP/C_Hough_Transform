#include "stdio.h"
#include "string.h"
#include "../ziku/showUtils.h"
#include "tistdtypes.h"

extern unsigned char OctArray[267616];

/*
函数名：showString
功能：显示字符串
参数：
	bmp 输入的RGB888像素矩阵
	base 字符串位于图像中的位置
	str 输入的字符串
注意：
	基准坐标base必须保证str不会写到图像外边去，否则数组越界；
	函数返回时，base的值会被重定位到字符串结尾处后一个字符位处
*/
void showString(GRAY8* bmp,Point* base, char* str)
{
	Uint16 width = bmp->width;
	//Uint16 height = bmp->height;
	Uint8* rgb = bmp->gray;
	Uint8* ptr;
	Uint16 a;
	int count=0;
	int i,j,k;
	Uint8 AB,CD;
	unsigned long addr;
	Uint8* tempOct;
	char* pstr = str;
	
	while (pstr[0]!='\0')
	{
		if (pstr[0]<=0x80)
		{
			if( pstr[0]>=33 && pstr[0]<=126 )
			{
				a = 0xa3a1 + pstr[0] - 33;
			}
			else
			{
				a = 0xa1a1;	//其它的都设置为空格
			}
			AB = (Uint8)(a>>8) & 0xff;
			CD = (Uint8)a & 0xff;
			pstr++;
		}
		else
		{
			AB = pstr[0];
			CD = pstr[1];
			pstr += 2;
		}
		//定位
		addr = ((AB-0xa1)*94+(CD-0xa1))*32;
		tempOct = OctArray + addr;
		for(i = 0;i < 16;i++)
		{
			for(j = 0;j < 2;j++)
			{
				for(k = 0;k < 8;k++)
				{
					//这里是针对RGB888图像的显示，
					//针对其它格式的图像只需修改下面这段像素填充代码即可
					//ptr = rgb + 3*((base->y+i)*width + base->x+16*count+(8*j+k));
					ptr = rgb + 1*((base->y+i)*width + base->x+16*count+(8*j+k));
					if(tempOct[i*2+j]&(0x80>>k))
					{
						//ptr[2] = 255 ;
						//ptr[1] = 255 ;
						ptr[0] = 255 ;
					}
					else
					{
						//ptr[2] = ptr[1] = ptr[0] = 0;
						ptr[0] = 0;
					}
				}
			}
		}
		count++;
	}
	//注意这个地方，为了方便后续显示，对base进行了修正
	base->x += 16*count;
}


/*
函数名：showCoordinate
功能：在bmp图的base位置显示坐标值point,形式如“(120,100)”
*/
void showCoordinate(GRAY8* bmp, Point* base, Point* point)
{

	char result[16]={0};
	
	sprintf(result,"(%d,%d)",point->x,point->y);
	
	showString(bmp,base,result);
	
}




