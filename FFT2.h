#pragma once
#include "iostream"
#include <fstream>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define Maskcols 128
#define Maskrows 256
using namespace std;


//�²��� ֱ����С4��
void pyr_downNM(unsigned int rows, unsigned int cols, int* image, int* new_image, unsigned int n, unsigned int m)
{
	unsigned int width2 = rows / m; //630/3= 210
	for (unsigned int i = 0; i < cols; i = i + n)// 480   +4
	{
		unsigned int yy, ys;
		yy = i / n * width2;
		ys = i * rows;
		for (int unsigned j = 0; j < rows; j = j + m)
		{
			new_image[yy + j / m] = image[j + ys];
		}
	}
}


/***************��������*****************/
typedef struct Complex {
	float real;
	float imagin;
};
//���帴�����㣬�����˷����ӷ�������
void Add_Complex(Complex * src1, Complex *src2, Complex *dst) {
	dst->imagin = src1->imagin + src2->imagin;
	dst->real = src1->real + src2->real;
}
void Sub_Complex(Complex * src1, Complex *src2, Complex *dst) {
	dst->imagin = src1->imagin - src2->imagin;
	dst->real = src1->real - src2->real;
}
void Multy_Complex(Complex * src1, Complex *src2, Complex *dst) {
	float r1 = 0.0, r2 = 0.0;
	float i1 = 0.0, i2 = 0.0;
	r1 = src1->real;
	r2 = src2->real;
	i1 = src1->imagin;
	i2 = src2->imagin;
	dst->imagin = r1 * i2 + r2 * i1;
	dst->real = r1 * r2 - i1 * i2;
}
//exp(-j2pi/N)
void getWN(float n, int size_n, Complex * dst) {
	float x = 2.0*M_PI*n / size_n;
	dst->imagin = -sin(x);
	dst->real = cos(x);
}
//����DFT����������DFT����
void DFT(float * src, Complex * dst, int size) {
	for (int m = 0; m < size; m++) {
		float real = 0.0;
		float imagin = 0.0;
		for (int n = 0; n < size; n++) {
			float x = M_PI * 2 * m*n;
			real += src[n] * cos(x / size);
			imagin += src[n] * (-sin(x / size));
		}
		dst[m].imagin = imagin;
		dst[m].real = real;
	}
}
//����IDFT����
void IDFT(Complex *src, Complex *dst, int size) {
	for (int m = 0; m < size; m++) {
		float real = 0.0;
		float imagin = 0.0;
		for (int n = 0; n < size; n++) {
			float x = M_PI * 2 * m*n / size;
			real += src[n].real*cos(x) - src[n].imagin*sin(x);
			imagin += src[n].real*sin(x) + src[n].imagin*cos(x);
		}
		real /= size;
		imagin /= size;
		if (dst != NULL) {
			dst[m].real = real;
			dst[m].imagin = imagin;
		}
	}
}
//��������
int FFT_remap(float * src, int size_n) {
	if (size_n == 1)
		return 0;
	float * temp = (float *)malloc(sizeof(float)*size_n);
	for (int i = 0; i < size_n; i++)
		if (i % 2 == 0)
			temp[i / 2] = src[i];
		else
			temp[(size_n + i) / 2] = src[i];
	for (int i = 0; i < size_n; i++)
		src[i] = temp[i];
	free(temp);
	FFT_remap(src, size_n / 2);
	FFT_remap(src + size_n / 2, size_n / 2);
	return 1;
}
int FFT_remap(Complex * src, int size_n) {
	if (size_n == 1)
		return 0;
	Complex * temp = (Complex *)malloc(sizeof(Complex)*size_n);
	for (int i = 0; i < size_n; i++)
		if (i % 2 == 0)
			temp[i / 2] = src[i];
		else
			temp[(size_n + i) / 2] = src[i];
	for (int i = 0; i < size_n; i++)
		src[i] = temp[i];
	free(temp);
	FFT_remap(src, size_n / 2);
	FFT_remap(src + size_n / 2, size_n / 2);
	return 1;
}
//����FFT
void FFT(float * src, Complex * dst, int size_n) {

	FFT_remap(src, size_n);
	int k = size_n;
	int z = 0;
	while (k /= 2) {
		z++;
	}
	k = z;
	if (size_n != (1 << k))
		exit(0);
	Complex * src_com = (Complex*)malloc(sizeof(Complex)*size_n);
	if (src_com == NULL)
		exit(0);
	for (int i = 0; i < size_n; i++) {
		src_com[i].real = src[i];
		src_com[i].imagin = 0;
	}
	for (int i = 0; i < k; i++) {
		z = 0;
		for (int j = 0; j < size_n; j++) {
			if ((j / (1 << i)) % 2 == 1) {
				Complex wn;
				getWN(z, size_n, &wn);
				Multy_Complex(&src_com[j], &wn, &src_com[j]);
				z += 1 << (k - i - 1);
				Complex temp;
				int neighbour = j - (1 << (i));
				temp.real = src_com[neighbour].real;
				temp.imagin = src_com[neighbour].imagin;
				Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
				Sub_Complex(&temp, &src_com[j], &src_com[j]);
			}
			else
				z = 0;
		}
	}
	for (int i = 0; i < size_n; i++) {
		dst[i].imagin = src_com[i].imagin;
		dst[i].real = src_com[i].real;
	}
}
void FFT(Complex * src, Complex * dst, int size_n) {

	FFT_remap(src, size_n);
	int k = size_n;
	int z = 0;
	while (k /= 2) {
		z++;
	}
	k = z;
	if (size_n != (1 << k))
		exit(0);
	Complex * src_com = (Complex*)malloc(sizeof(Complex)*size_n);
	if (src_com == NULL)
		exit(0);
	for (int i = 0; i < size_n; i++) {
		src_com[i].real = src[i].real;
		src_com[i].imagin = src[i].imagin;
	}
	for (int i = 0; i < k; i++) {
		z = 0;
		for (int j = 0; j < size_n; j++) {
			if ((j / (1 << i)) % 2 == 1) {
				Complex wn;
				getWN(z, size_n, &wn);
				Multy_Complex(&src_com[j], &wn, &src_com[j]);
				z += 1 << (k - i - 1);
				Complex temp;
				int neighbour = j - (1 << (i));
				temp.real = src_com[neighbour].real;
				temp.imagin = src_com[neighbour].imagin;
				Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
				Sub_Complex(&temp, &src_com[j], &src_com[j]);
			}
			else
				z = 0;
		}
	}
	for (int i = 0; i < size_n; i++) {
		dst[i].imagin = src_com[i].imagin;
		dst[i].real = src_com[i].real;
	}
}
void IFFT(Complex * src, Complex * dst, int size_n) {

	FFT_remap(src, size_n);
	int k = size_n;
	int z = 0;
	while (k /= 2) {
		z++;
	}
	k = z;
	if (size_n != (1 << k))
		exit(0);
	Complex * src_com = (Complex*)malloc(sizeof(Complex)*size_n);
	if (src_com == NULL)
		exit(0);
	for (int i = 0; i < size_n; i++) {
		src_com[i].real = src[i].real;
		src_com[i].imagin = src[i].imagin;
	}
	for (int i = 0; i < k; i++) {
		z = 0;
		for (int j = 0; j < size_n; j++) {
			if ((j / (1 << i)) % 2 == 1) {
				Complex wn;
				getWN(-1.0*z, size_n, &wn);
				Multy_Complex(&src_com[j], &wn, &src_com[j]);
				z += 1 << (k - i - 1);
				Complex temp;
				int neighbour = j - (1 << (i));
				temp.real = src_com[neighbour].real;
				temp.imagin = src_com[neighbour].imagin;
				Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
				Sub_Complex(&temp, &src_com[j], &src_com[j]);
			}
			else
				z = 0;
		}
	}
	for (int i = 0; i < size_n; i++) {
		dst[i].imagin = src_com[i].imagin / size_n;
		dst[i].real = src_com[i].real / size_n;
	}
}


void FFT_2D(float** temp, Complex **get, int size_n)
{
	//  float temp[SIZE][SIZE] = { 0.0 };//�洢ͼ������
	float tempreal[Maskcols][Maskrows] = { 0.0 };//�洢ͼ������任���ʵ��
	float tempimage[Maskcols][Maskrows] = { 0.0 };//�洢ͼ������任����鲿
	Complex src[Maskrows];
	Complex dst[Maskrows];
	///��ͼ����и���Ҷ�任
	//����ÿһ�е�FFT�任
	for (int i = 0; i < Maskcols; i++)
	{
		FFT(temp[i], dst, Maskrows);//�õ�256���ȸ�������dst
		for (int j = 0; j < Maskrows; j++)
		{
			tempreal[i][j] = dst[j].real;//�洢ʵ��
			tempimage[i][j] = dst[j].imagin;//�洢�鲿
		}
	}
	//���Ѿ������б任��ĸ��������ٽ�����fft�任
	for (int j = 0; j < Maskrows; j++)
	{
		for (int i = 0; i < Maskcols; i++)//�и���
		{
			src[i].real = tempreal[i][j];
			src[i].imagin = tempimage[i][j];
		}//���ƺõ�j�е�src[SIZE]
		FFT(src, dst, Maskcols);
		for (int i = 0; i < Maskcols; i++)
		{
			get[i][j].real = dst[i].real;//�ѶԵ�j�н��е�fft�ٴ��ԭ��������
			get[i][j].imagin = dst[i].imagin;
		}
	}//�õ�����j�е�fft�任��Ҳ����ͼ���fft�任tempreal+j*tempimage = F(u,v)
}

void IFFT_2D(float** tempreal, float** tempimage, Complex **get, int size_n)
{
	/********************************************
	****************���з�����Ҷ�任*************
	*********************************************/
	//��i�У�ÿһ�н���һά����Ҷ�任
	Complex src[Maskrows];
	Complex dst[Maskrows];

	for (int i = 0; i < Maskcols; i++)
	{
		for (int j = 0; j < Maskrows; j++)
		{
			//��F(u,v)ÿһ�У������ǵ�i�д洢��Complex һά����src
			src[j].real = tempreal[i][j];
			src[j].imagin = tempimage[i][j];
		}
		IFFT(src, dst, Maskrows);
		//FFT_remap(dst, SIZE);
		for (int j = 0; j < Maskrows; j++)
		{
			tempreal[i][j] = dst[j].real;
			tempimage[i][j] = dst[j].imagin;
		}
	}
	//���ڽ����б仯
	for (int j = 0; j < Maskrows; j++)
	{
		for (int i = 0; i < Maskcols; i++)
		{
			src[i].real = tempreal[i][j];
			src[i].imagin = tempimage[i][j];
		}//���ƺ�ÿһ��
		IFFT(src, dst, Maskcols);
		//FFT_remap(dst, SIZE);
		for (int i = 0; i < Maskcols; i++)
		{
			get[i][j].real = dst[i].real;           //tempreal[i][j] = dst[i].real;
			get[i][j].imagin = dst[i].imagin;       //tempimage[i][j] = dst[i].imagin;
		}
	}
}
