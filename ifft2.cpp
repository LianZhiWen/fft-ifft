#include"FFT2.h"


//�²��� ֱ����С4��
void pyr_down(unsigned int rows, unsigned int cols, int* image, int* new_image, unsigned int n, unsigned int m)
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

//void IFFT_2D(float** tempreal, float** tempimage, Complex **get, int size_n)
//{
//	/********************************************
//	****************���з�����Ҷ�任*************
//	*********************************************/
//	//��i�У�ÿһ�н���һά����Ҷ�任
//	Complex src[Maskrows];
//	Complex dst[Maskrows];
//
//	for (int i = 0; i < Maskcols; i++)
//	{
//		for (int j = 0; j < Maskrows; j++)
//		{
//			//��F(u,v)ÿһ�У������ǵ�i�д洢��Complex һά����src
//			src[j].real = tempreal[i][j];
//			src[j].imagin = tempimage[i][j];
//		}
//		IFFT(src, dst, Maskrows);
//		//FFT_remap(dst, SIZE);
//		for (int j = 0; j < Maskrows; j++)
//		{
//			tempreal[i][j] = dst[j].real;
//			tempimage[i][j] = dst[j].imagin;
//		}
//	}
//	//���ڽ����б仯
//	for (int j = 0; j < Maskrows; j++)
//	{
//		for (int i = 0; i < Maskcols; i++)
//		{
//			src[i].real = tempreal[i][j];
//			src[i].imagin = tempimage[i][j];
//		}//���ƺ�ÿһ��
//		IFFT(src, dst, Maskcols);
//		//FFT_remap(dst, SIZE);
//		for (int i = 0; i < Maskcols; i++)
//		{
//			get[i][j].real = dst[i].real;           //tempreal[i][j] = dst[i].real;
//			get[i][j].imagin = dst[i].imagin;       //tempimage[i][j] = dst[i].imagin;
//		}
//	}
//}



int main(int argc, const char * argv[]) {


	//const char* imagename = "D://Asef//����//1_33.bmp";//�˴�Ϊ���Լ���ͼƬ·��

	////���ļ��ж���ͼ��

	//ifstream infileR, infileX;
	//infileR.open("RR.txt");//���ļ� �˲���ʵ��
	//infileX.open("XX.txt");//���ļ� �˲����鲿

	//Maskimage = Readfilter(infileR, infileX, , maskrows, maskcols);
	FILE* fpp;
	errno_t errr;
	errr = fopen_s(&fpp, "D:\\Asef_eyeC\\ceshi\\4.bmp", "rb");
	fseek(fpp, 1078, 0);
	unsigned char* pData = new unsigned char[480 * 640];          //����һ��������
	fread(pData, sizeof(unsigned char) * (480 * 640), 1, fpp);    //��ȡͼƬ���� fp ->pData


	int width = 480;
	int length = 640;

	int* P = new int[width*length];
	memset(P, 0, width*length);
	int* PP = new int[width*(length - 5)];
	memset(PP, 0, 480 * 635);
	//˳����
	for (int i = width - 1; i >= 0; i--)
	{
		int m = (width - 1 - i) * length;
		int n = i * length;
		for (int j = 0; j < length; j++)
		{
			P[m + j] = pData[n + j];
		}
	}


	for (int i = 0; i < 480; i++)
	{
		for (int j = 5; j < 635; j++)
		{
			PP[i * 630 + j - 5] = P[i * 640 + j];  //��ʾ���

		}

	}


	int* NewPic = new int[120 * 210];
	memset(NewPic, 0, 120 * 210);
	pyr_down(630, 480, PP, NewPic, 4, 3);

	float tempp[120][210] = { 0.0 };//����ͼ��

	for (int i = 0; i < 120; i++)
	{
		for (int j = 0; j < 210; j++)
		{
			tempp[i][j] = NewPic[i * 210 + j];  //��ʾ���
		}

	}

	float temp[128][256] = { 0.0 };//����ͼ��

	for (int i = 4; i < 124; i++)
	{
		for (int j = 23; j < 233; j++)
		{
			temp[i][j] = tempp[i - 4][j - 23];
		}

	}


	float tempreal[Maskcols][Maskrows] = { 0.0 };//�洢ͼ������任���ʵ��
	float tempimage[Maskcols][Maskrows] = { 0.0 };//�洢ͼ������任����鲿
	Complex src[Maskrows];
	Complex dst[Maskrows];

	/****************************************/
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
			tempreal[i][j] = dst[i].real;//�ѶԵ�j�н��е�fft�ٴ��ԭ��������
			tempimage[i][j] = dst[i].imagin;
		}
	}//�õ�����j�е�fft�任��Ҳ����ͼ���fft�任tempreal+j*tempimage = F(u,v)


	//cout << "\n\n**************************************\n\n";
	//for (int i = 0; i < 128; i++)
	//{
	//	for (int j = 0; j < 256; j++)
	//	{
	//		cout  << '\t' << tempreal[i][j] << '\t' << tempimage[i][j] ;

	//	}
	//	cout << " ----------"<<endl;
	//}
	float **IfftR = new float *[Maskcols];
	for (int i = 0; i < Maskrows; i++) {
		IfftR[i] = new float[Maskrows];
	}
	for (int i = 0; i < Maskcols; i++)
	{
		for (int j = 0; j < Maskrows; j++)
		{
			IfftR[i][j] = tempreal[i][j];  //��ʾ���
		}

	}
	float **IfftX = new float *[Maskcols];
	for (int i = 0; i < Maskrows; i++) {
		IfftX[i] = new float[Maskrows];
	}
	for (int i = 0; i < Maskcols; i++)
	{
		for (int j = 0; j < Maskrows; j++)
		{
			IfftX[i][j] = tempimage[i][j];  //��ʾ���
		}

	}
	Complex **MaskIFFT = new Complex*[Maskcols];
	//��ʼ����ά�ṹ��  ��任�Ľ��
	for (int i = 0; i < Maskcols; i++) {
		MaskIFFT[i] = new Complex[Maskrows];
	}
	for (int i = 0; i < Maskcols; i++)
	{

		for (int j = 0; j < Maskrows; j++)
		{
			MaskIFFT[i][j] = { 0.0 };  //��ʾ���
		}

	}

	IFFT_2D(IfftR, IfftX, MaskIFFT, 256);
	for (int i = 4; i < 124; i++)
	{	
		for (int j =23 ; j < 233; j++)
		{
			cout<<MaskIFFT[i][j].real<<" ,";
		}
		cout << endl;
	}

	/*******************����һ�����任���㷨*******************
	****************���з�����Ҷ�任*************
	*********************************************/
	//��i�У�ÿһ�н���һά����Ҷ�任
	//for (int i = 0; i < Maskcols; i++)
	//{
	//	for (int j = 0; j < Maskrows; j++)
	//	{
	//		//��F(u,v)ÿһ�У������ǵ�i�д洢��Complex һά����src
	//		src[j].real = tempreal[i][j];
	//		src[j].imagin = tempimage[i][j];
	//	}
	//	IFFT(src, dst, Maskrows);
	//	//FFT_remap(dst, Maskrows);
	//	for (int j = 0; j < Maskrows; j++)
	//	{
	//		tempreal[i][j] = dst[j].real;
	//		tempimage[i][j] = dst[j].imagin;
	//	}
	//}
	////���ڽ����б仯
	//for (int j = 0; j < Maskrows; j++)
	//{
	//	for (int i = 0; i < Maskcols; i++)
	//	{
	//		src[i].real = tempreal[i][j];
	//		src[i].imagin = tempimage[i][j];
	//	}//���ƺ�ÿһ��
	//	IFFT(src, dst, Maskcols);
	//	//FFT_remap(dst, Maskcols);
	//	for (int i = 0; i < Maskcols; i++)
	//	{
	//		tempreal[i][j] = dst[i].real;
	//		tempimage[i][j] = dst[i].imagin;
	//	}
	//}
	///***************************************/


	return 0;
}