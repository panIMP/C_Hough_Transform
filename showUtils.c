#include "stdio.h"
#include "string.h"
#include "../ziku/showUtils.h"
#include "tistdtypes.h"

extern unsigned char OctArray[267616];

/*
��������showString
���ܣ���ʾ�ַ���
������
	bmp �����RGB888���ؾ���
	base �ַ���λ��ͼ���е�λ��
	str ������ַ���
ע�⣺
	��׼����base���뱣֤str����д��ͼ�����ȥ����������Խ�磻
	��������ʱ��base��ֵ�ᱻ�ض�λ���ַ�����β����һ���ַ�λ��
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
				a = 0xa1a1;	//�����Ķ�����Ϊ�ո�
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
		//��λ
		addr = ((AB-0xa1)*94+(CD-0xa1))*32;
		tempOct = OctArray + addr;
		for(i = 0;i < 16;i++)
		{
			for(j = 0;j < 2;j++)
			{
				for(k = 0;k < 8;k++)
				{
					//���������RGB888ͼ�����ʾ��
					//���������ʽ��ͼ��ֻ���޸�����������������뼴��
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
	//ע������ط���Ϊ�˷��������ʾ����base����������
	base->x += 16*count;
}


/*
��������showCoordinate
���ܣ���bmpͼ��baseλ����ʾ����ֵpoint,��ʽ�硰(120,100)��
*/
void showCoordinate(GRAY8* bmp, Point* base, Point* point)
{

	char result[16]={0};
	
	sprintf(result,"(%d,%d)",point->x,point->y);
	
	showString(bmp,base,result);
	
}




