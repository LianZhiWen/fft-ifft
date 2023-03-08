#include"FFT2.h"


//下采样 直接缩小4倍
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
//	****************进行反傅里叶变换*************
//	*********************************************/
//	//对i行，每一行进行一维傅里叶变换
//	Complex src[Maskrows];
//	Complex dst[Maskrows];
//
//	for (int i = 0; i < Maskcols; i++)
//	{
//		for (int j = 0; j < Maskrows; j++)
//		{
//			//把F(u,v)每一行，这里是第i行存储到Complex 一维数组src
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
//	//现在进行列变化
//	for (int j = 0; j < Maskrows; j++)
//	{
//		for (int i = 0; i < Maskcols; i++)
//		{
//			src[i].real = tempreal[i][j];
//			src[i].imagin = tempimage[i][j];
//		}//复制好每一列
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


	//const char* imagename = "D://Asef//测试//1_33.bmp";//此处为你自己的图片路径

	////从文件中读入图像

	//ifstream infileR, infileX;
	//infileR.open("RR.txt");//打开文件 滤波器实部
	//infileX.open("XX.txt");//打开文件 滤波器虚部

	//Maskimage = Readfilter(infileR, infileX, , maskrows, maskcols);
	FILE* fpp;
	errno_t errr;
	errr = fopen_s(&fpp, "D:\\Asef_eyeC\\ceshi\\4.bmp", "rb");
	fseek(fpp, 1078, 0);
	unsigned char* pData = new unsigned char[480 * 640];          //创建一个空数组
	fread(pData, sizeof(unsigned char) * (480 * 640), 1, fpp);    //读取图片数据 fp ->pData


	int width = 480;
	int length = 640;

	int* P = new int[width*length];
	memset(P, 0, width*length);
	int* PP = new int[width*(length - 5)];
	memset(PP, 0, 480 * 635);
	//顺序存放
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
			PP[i * 630 + j - 5] = P[i * 640 + j];  //显示结果

		}

	}


	int* NewPic = new int[120 * 210];
	memset(NewPic, 0, 120 * 210);
	pyr_down(630, 480, PP, NewPic, 4, 3);

	float tempp[120][210] = { 0.0 };//输入图像

	for (int i = 0; i < 120; i++)
	{
		for (int j = 0; j < 210; j++)
		{
			tempp[i][j] = NewPic[i * 210 + j];  //显示结果
		}

	}

	float temp[128][256] = { 0.0 };//输入图像

	for (int i = 4; i < 124; i++)
	{
		for (int j = 23; j < 233; j++)
		{
			temp[i][j] = tempp[i - 4][j - 23];
		}

	}


	float tempreal[Maskcols][Maskrows] = { 0.0 };//存储图像数组变换后的实部
	float tempimage[Maskcols][Maskrows] = { 0.0 };//存储图像数组变换后的虚部
	Complex src[Maskrows];
	Complex dst[Maskrows];

	/****************************************/
	///对图像进行傅里叶变换

	//进行每一行的FFT变换
	for (int i = 0; i < Maskcols; i++)
	{

		FFT(temp[i], dst, Maskrows);//得到256长度复数数组dst

		for (int j = 0; j < Maskrows; j++)
		{
			tempreal[i][j] = dst[j].real;//存储实部
			tempimage[i][j] = dst[j].imagin;//存储虚部
		}
	}
	//对已经进行行变换后的复数数组再进行列fft变换
	for (int j = 0; j < Maskrows; j++)
	{
		for (int i = 0; i < Maskcols; i++)//列复制
		{
			src[i].real = tempreal[i][j];
			src[i].imagin = tempimage[i][j];
		}//复制好第j列到src[SIZE]
		FFT(src, dst, Maskcols);
		for (int i = 0; i < Maskcols; i++)
		{
			tempreal[i][j] = dst[i].real;//把对第j列进行的fft再次填到原来的数组
			tempimage[i][j] = dst[i].imagin;
		}
	}//得到各个j列的fft变换，也就是图像的fft变换tempreal+j*tempimage = F(u,v)


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
			IfftR[i][j] = tempreal[i][j];  //显示结果
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
			IfftX[i][j] = tempimage[i][j];  //显示结果
		}

	}
	Complex **MaskIFFT = new Complex*[Maskcols];
	//初始化二维结构体  逆变换的结果
	for (int i = 0; i < Maskcols; i++) {
		MaskIFFT[i] = new Complex[Maskrows];
	}
	for (int i = 0; i < Maskcols; i++)
	{

		for (int j = 0; j < Maskrows; j++)
		{
			MaskIFFT[i][j] = { 0.0 };  //显示结果
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

	/*******************附加一个反变换的算法*******************
	****************进行反傅里叶变换*************
	*********************************************/
	//对i行，每一行进行一维傅里叶变换
	//for (int i = 0; i < Maskcols; i++)
	//{
	//	for (int j = 0; j < Maskrows; j++)
	//	{
	//		//把F(u,v)每一行，这里是第i行存储到Complex 一维数组src
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
	////现在进行列变化
	//for (int j = 0; j < Maskrows; j++)
	//{
	//	for (int i = 0; i < Maskcols; i++)
	//	{
	//		src[i].real = tempreal[i][j];
	//		src[i].imagin = tempimage[i][j];
	//	}//复制好每一列
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