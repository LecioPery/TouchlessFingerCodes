#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <math.h>
#include <string.h>

#include <matrix.h>
#include <mat.h>

#include <opencv2/highgui/highgui.hpp>

#include "..\Header Files\trickyLib.hpp"
#include "..\Header Files\commonConstants.hpp"

using namespace cv;
using namespace std;

void matParse( cv::Mat &Target, const char *varname )
{

	/** Accept only char type matrices. */
	CV_Assert( Target.depth() != sizeof(uchar) );

	int channels = Target.channels();
	int nRows = Target.rows;
	int nCols = Target.cols;

	cout << "Started the parse '" << varname << "'..." << endl;
	cout << "channels: " << channels << " nRows: " << nRows << " nCols: " << nCols << endl;

	switch ( channels )
	{

		case 1:
			cout << Target;
		break;

		case 3:
			int i, j, k;
			float *p;
			for (k = 0; k < channels; k++)
			{
				cout << "channel: " << k << endl;
				for (i = 0; i < nRows; ++i)
				{
					p = Target.ptr<float>(i);
					for (j = k * nCols; j < (k + 1) * nCols; ++j)
					{
						cout << p[j] << ", ";
					}
					cout << endl;
				}
			}

			/** Alternate version: (SIMILAR - not equal - result)
			unsigned char b;
			unsigned char g;
			unsigned char r;

			for (int j = 0; j < nRows; j++){
				for (int i = 0; i < nCols; i++){
					b = Target.data[Target.step * j + i];
					g = Target.data[Target.step * j + i + 1];
					r = Target.data[Target.step * j + i + 2];
				}
			}
			cout << b << ", " << g << ", " << r << endl;
			*/
		break;

		case 4:
			cout << "TODO." << endl;
		break;
		
		default:
			cout << "Error! Number of channels is: " << channels << endl;
			return;
		break;

	}

	cout << endl << "Ended the parse of '" << varname << "'!" << endl;

}

cv::Mat setMatValuesFloat( cv::Mat img, cv::Range heightInterval, cv::Range widthInterval, float value )
{
	
	/** Goal: do anything fancy with the image, so we know it's matrix has changed! */
	int channels = img.channels();
	int nCols = img.cols * channels;
	int r;
	int s;
	float *t;
	
	for ( r = heightInterval.start - 1; r < heightInterval.end; ++r )
	{

		t = img.ptr<float>(r);
		for ( s = widthInterval.start - 1; s < widthInterval.end; ++s )
		{
			
			t[s] = value;

		}

	}

	return img;

}

cv::string getImageType(int number)
{
	// find type
	int imgTypeInt = number % 8;
	cv::string imgTypeString;

	switch (imgTypeInt)
	{
	case 0:
		imgTypeString = "8U";
		break;
	case 1:
		imgTypeString = "8S";
		break;
	case 2:
		imgTypeString = "16U";
		break;
	case 3:
		imgTypeString = "16S";
		break;
	case 4:
		imgTypeString = "32S";
		break;
	case 5:
		imgTypeString = "32F";
		break;
	case 6:
		imgTypeString = "64F";
		break;
	default:
		break;
	}

	// find channel
	int channel = (number / 8) + 1;

	std::stringstream type;
	type << "CV_" << imgTypeString << "C" << channel;

	return type.str();
}

float mean(float *data, int n)
{
	float mean = 0.0;
	int i;
	for (i = 0; i<n; ++i)
	{
		mean += data[i];
	}
	mean = mean / n;
	return mean;
}

float standard_deviation(float *data, int n)
{
	float mean = 0.0, sum_deviation = 0.0;
	int i;
	for (i = 0; i<n; ++i)
	{
		mean += data[i];
	}
	
	mean /= n;
	
	for (i = 0; i<n; ++i)
		sum_deviation += (data[i] - mean)*(data[i] - mean);

	return sqrt(sum_deviation / (i - 1));
}


float standard_deviation_vector(vector<float> data)
{
	float mean = 0.0, sum_deviation = 0.0;
	int i;
	for (i = 0; i<data.size(); ++i)
	{
		mean += data[i];
	}

	mean /= data.size();

	for (i = 0; i<data.size(); ++i)
		sum_deviation += (data[i] - mean)*(data[i] - mean);

	return sqrt(sum_deviation / (i-1));
}


float std2(Mat imagem)
{
	float result;

	int nRows;
	int nCols;
	int channels;
	int conta=0;

	nRows = imagem.rows;
	nCols = imagem.cols;
	channels = imagem.channels();

	
	vector<float> bloco;

	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			bloco.push_back(imagem.at<float>(i,j));
		}		
	}

	result = standard_deviation_vector(bloco);
	

	//result = standard_deviation(imagem.ptr<float>(0), nCols*nRows*channels);


	return result;

}

/** Usado por stddvim(backup). Descontinuado.
Mat extractIndexes(Mat target, int *lineIndex, int hLenght, int *columnIndex, int cLenght)
{
	Mat result;

	/** Retira da matriz apenas os índices correspondentes ao produto cartesiano:
	linexIndex X columnIndex*//*
	int channels = target.channels();
	int r, s;
	uchar *t, *u;

	u = result.ptr<uchar>(0);

	for (r = lineIndex[0]; r < hLenght; ++r)
	{

		t = target.ptr<uchar>(r);
		for (s = columnIndex[0]; s < cLenght; ++s)
		{

			u[r * hLenght + s] = t[s];

		}

	}

	return result;
}
*/

cv::Mat thresholdImage(cv::Mat img, float thresh)
{

	/** Objetivo: produzir uma imagem binária a partir de um threshold! */
	Mat result;
	int channels = img.channels();
	int nRows = img.rows;
	int nCols = img.cols * channels;
	int r, s;
	float *t;

	img.copyTo(result);

	for (r = 0; r < nRows; ++r)
	{

		t = result.ptr<float>(r);
		for (s = 0; s < nCols; ++s)
		{

			if (t[s] - thresh < 0)
			{
				t[s] = 0;
			}
			else t[s] = 255;

		}

	}
	
	return result;

}

/**
Normaliza a imagem. Faz os níveis de cinza se distribuírem entre 0 e 255.
\param [in] src imagem de entrada.
\return Imagem normalizada.
*/

Mat normalizaImgDebug(const Mat& src)
{
	Mat dst;
	dst = src.clone();
	
	/* Determina o máximo e o mínimo da imagem corrigida para normalizar a imagem. */
	double min, max;
	minMaxLoc(src, &min, &max);

	/* Subtrai o mínimo */
	for (int i = 0; i < src.rows; i++){
		for (int j = 0; j < src.cols; j++){
			dst.at<float>(i, j) = (float) (src.at<float>(i, j)-min);
		}
	}
	minMaxLoc(dst, &min, &max);
	/* Normaliza a imagem dividindo pelo máximo e multiplicando por 255, tornando-a visível */
	for (int i = 0; i < src.rows; i++){
		for (int j = 0; j < src.cols; j++){
			dst.at<float>(i, j) = (float) (255*(dst.at<float>(i, j)/max));
		}
	}

	// Converte a imageOriginalPow de float para CV_8U
	dst.convertTo(dst, CV_8U);

	return dst;

}

cv::Mat bytewiseUpdateImage(cv::Mat target, float value, int METHOD)
{

	/** Objetivo: produzir uma imagem binária a partir de um threshold! */
	Mat result;
	int channels = target.channels();
	int nRows = target.rows;
	int nCols = target.cols * channels;
	int r, s;
	float *t;

	target.copyTo(result);

	if( METHOD == METHOD_SUM )
	{

		for (r = 0; r < nRows; ++r)
		{

			t = result.ptr<float>(r);
			for (s = 0; s < nCols; ++s)
			{

				t[s] = t[s] + value;

			}

		}

	}
	else if ( METHOD == METHOD_SUBTRACT )
	{

		for (r = 0; r < nRows; ++r)
		{

			t = result.ptr<float>(r);
			for (s = 0; s < nCols; ++s)
			{

				t[s] = t[s] - value;

			}

		}

	}
	else if ( METHOD == METHOD_PRODUCT )
	{

		for (r = 0; r < nRows; ++r)
		{

			t = result.ptr<float>(r);
			for (s = 0; s < nCols; ++s)
			{

				t[s] = (float) t[s] * value;

			}

		}

	}
	else if ( METHOD == METHOD_DIVIDE )
	{

		for (r = 0; r < nRows; ++r)
		{

			t = result.ptr<float>(r);
			for (s = 0; s < nCols; ++s)
			{

				t[s] = (float) t[s] / value;

			}

		}

	}
	return result;

}

void matReadDouble(const char *file, std::vector<double>& v, const char *varName)
{
	
	// open MAT-file
	MATFile *pmat = matOpen(file, "r");
	if ( pmat == NULL ) return;

	// extract the specified variable
	mxArray *arr = matGetVariable(pmat, varName);
	if ( arr != NULL && mxIsDouble(arr) && !mxIsEmpty(arr) )
	{
		
		// copy data
		mwSize num = mxGetNumberOfElements(arr);
		double *pr = mxGetPr(arr);
		if (pr != NULL)
		{
			
			v.resize(num);
			v.assign(pr, pr + num);
			
		}
		
	}

	// cleanup
	mxDestroyArray(arr);
	matClose(pmat);
	
}

void matReadMat( const char *file, cv::Mat& v, const char *varName )
{
	
	int counter;
	int i;
	int j;
	MATFile *pmat = matOpen(file, "r");

	// open MAT-file
	pmat = matOpen(file, "r");
	if (pmat == NULL)
	{
		std::cout << "Could not open " << file << "." << std::endl;
		return;
	}
	
	// extract the specified variable
	mxArray *arr = matGetVariable(pmat, varName);
	if (arr == NULL)
	{
		std::cout << "Variable " << varName << "not found." << std::endl;
		return;
	}

	/* size[0] = rows, size[1] = cols */
	const size_t *size = mxGetDimensions(arr);
	v = Mat( (int) size[0], (int) size[1], CV_32F);

	counter = 0;
	if (arr != NULL && mxIsDouble(arr) && !mxIsEmpty(arr))
	{

		// copy data
		mwSize num = mxGetNumberOfElements(arr);
		double *pr = mxGetPr(arr);
		if (pr != NULL)
		{
			for (i = 0; i < size[0]; i++)
			{
				for (j = 0; j < size[1]; j++)
				{
					
					/* Assigning value by value. Notice OpenCV swaps columns with lines. */
					counter++;
					/** For unknown reasons, both iplementations below fails:
					 * v.at<float>(j, i) = (float) *(pr + counter);
					 * v.at<float>(counter) = (float) *(pr + counter);
					**/
					v.at<float>(j, i) = (float) pr[i * size[1] + j];
					
				}
				
			}
			
		}

	}
	
	// cleanup
	mxDestroyArray(arr);
	matClose(pmat);
}

void printVectorDouble(const char *vName, std::vector<double> v, int pauseMode, int indexes)
{

	std::cout << "--- Begin " << vName << " ---" << std::endl;
	for (size_t i = 0; i < v.size(); ++i)
	{

		if ( indexes == INDEXES_ON )
		{
			std::cout << i + 1 << ": " << v[i];
		}
		else std::cout << v[i];

		if ( (i + 1) % 6 == 0 )
		{
			std::cout << std::endl;
		}
		else std::cout << " ";

	}
	std::cout << std::endl << "--- End " << vName << " ---" << std::endl;
	if ( pauseMode == PAUSE_ON )
	{
		system("PAUSE");
	}
	else std::cout << "		No pause." << std::endl;

	return;

}

cv::Mat matProduct( cv::Mat input1, cv::Mat input2, const int METHOD )
{
	
	cv::Mat result;
	int k;
	if ( ( input1.size().height == input2.size().height ) && ( input1.size().width == input2.size().width ) )
	{

		result = cv::Mat( input1.size().height, input1.size().width, CV_32F );
		k = 0;
		while (k < result.size().area())
		{

			if ( METHOD == METHOD_PRODUCT )
			{
				result.at<float>(k) = input1.at<float>(k) * input2.at<float>(k);
			}
			else if ( METHOD == METHOD_DIVIDE )
			{
				result.at<float>(k) = input1.at<float>(k) / input2.at<float>(k);
			}
			else
			{
				cout << "Invalid method on matProduct();, returning" << std::endl;
				return result;
			}
			k++;

		}

	}
	else
	{

		std::cout << "WARNING: Dimentional sizes mismatch at product between matrices." << std::endl;

	}
	return result;

}

void matSinCos( cv::Mat Gxx, cv::Mat Gxy, cv::Mat Gyy, cv::Mat *sinTheta, cv::Mat *cosTheta )
{

	Mat denom;
	int k;
	if ( ( Gxx.size().height == Gxy.size().height ) && ( Gxy.size().height == Gyy.size().height ) &&
		 ( Gxx.size().width == Gxy.size().width ) && ( Gxy.size().width == Gyy.size().width ) )
	{

		k = 0;
		(*sinTheta) = cv::Mat( Gxx.rows, Gxx.cols, CV_32F );
		(*cosTheta) = cv::Mat( Gxx.rows, Gxx.cols, CV_32F );
		denom = cv::Mat( Gxx.rows, Gxx.cols, CV_32F );
		while (k < Gxx.size().area())
		{

			denom.at<float>(k) = sqrt( ( Gxy.at<float>(k) * Gxy.at<float>(k) ) +
						  ( Gxx.at<float>(k) - Gyy.at<float>(k) ) * ( Gxx.at<float>(k) - Gyy.at<float>(k) ) );
			(*sinTheta).at<float>(k) = Gxy.at<float>(k) / denom.at<float>(k);
			(*cosTheta).at<float>(k) = ( Gxx.at<float>(k) - Gyy.at<float>(k) ) / denom.at<float>(k);
			k++;

		}

	}
	else
	{

		std::cout << "WARNING: Dimentional sizes mismatch at product between matrices." << std::endl;

	}
	/*namedWindow("denom", CV_WINDOW_AUTOSIZE);
	imshow("denom", denom);*/
	return;
	
}

cv::Mat matAtan2( cv::Mat input1, cv::Mat input2 )
{
	
	Mat result;
	int k;
	
	/*
	const long double PI = 3.141592653589793238L;
	const double PI = 3.141592653589793;
	const float PI = 3.1415927;
	*/
	
	if ( ( input1.size().height == input2.size().height ) && ( input1.size().width == input2.size().width ) )
	{
		
		result = Mat( input1.size().height, input1.size().width, CV_32F );
		k = 0;
		while ( k < input1.size().area() )
		{
			
			result.at<float>( k ) = ( ( float ) 3.1415927 / 2 ) + ( atan2( input1.at<float>( k ), input2.at<float>( k ) ) / 2 );
			k++;
			
		}
		
	}
	return result;
	
}

Mat matFunction( Mat target, const char *function )
{
	
	Mat result;
	target.copyTo(result);
	int r;
	int s;
	float *t;
	
	if ( !strcmp(function, "sin") )
	{
		
		for (r = 0; r < result.rows; ++r)
		{
			
			t = result.ptr<float>(r);
			for (s = 0; s < result.cols; ++s)
			{
				
				t[s] = ( float ) sin(t[s]);
				
			}
			
		}
		
	}
	return result;
	
}

Mat gradientMat( Mat TargetImage )
{
	
	int k;
	int j;
	
	Mat Gx = Mat::zeros(TargetImage.rows, TargetImage.cols, CV_32F);;
	Mat Gy = Mat::zeros(TargetImage.rows, TargetImage.cols, CV_32F);;
	Mat result = Mat::zeros(TargetImage.rows, TargetImage.cols, CV_32F);

	Sobel(TargetImage, Gx, CV_32F, 1, 0, 3, 1, 0, BORDER_DEFAULT);
	Sobel(TargetImage, Gy, CV_32F, 0, 1, 3, 1, 0, BORDER_DEFAULT);
	//addWeighted(grad_x, 0.5, grad_y, 0.5, 0, grad);

	
	/* Gradient in A, bad performance - enhance with symmetry + memory access: */
	for (k = 1; k < TargetImage.rows - 1; k++)
	{
		
		for (j = 1; j < TargetImage.cols - 1; j++)
		{
			
			/*Gx = ( TargetImage.at<float>(k + 1, j - 1) + 2 * TargetImage.at<float>(k + 1, j) + TargetImage.at<float>(k + 1, j + 1) )
			   - ( TargetImage.at<float>(k - 1, j - 1) + 2 * TargetImage.at<float>(k - 1, j) + TargetImage.at<float>(k - 1, j + 1) );
			
			Gy = ( TargetImage.at<float>(k + 1, j + 1) + 2 * TargetImage.at<float>(k, j + 1) + TargetImage.at<float>(k - 1, j + 1) )
			   - ( TargetImage.at<float>(k + 1, j - 1) + 2 * TargetImage.at<float>(k, j - 1) + TargetImage.at<float>(k - 1, j - 1) );
			*/
			
			result.at<float>(k, j) = sqrt(Gx.at<float>(k, j) * Gx.at<float>(k, j) + Gy.at<float>(k, j) * Gy.at<float>(k, j));
			
		}
		
	}
	return result;
	
}

Mat tailMat( Mat TargetImage )
{
	
	Mat result;
	int i = 0;
	int j = 0;
	
	TargetImage.copyTo( result );
	while (j < TargetImage.cols)
	{
		result.at<float>(TargetImage.rows - 1, j) = 0.0;
		j++;
	}
	while (i < TargetImage.rows)
	{
		result.at<float>(i, TargetImage.cols - 1) = 0.0;
		i++;
	}
	return result;
	
}

Mat normalizaImg(const Mat& src)
{

	Mat dst;
	dst = src.clone();

	/* Determina o máximo e o mínimo da imagem corrigida para normalizar a imagem. */
	double min, max;
	minMaxLoc(dst, &min, &max);

	/* Normaliza a imagem dividindo pelo máximo e multiplicando por 255, tornando-a visível */
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			dst.at<float>(i, j) = (float) round( 255 * ( dst.at<float>(i, j)/max ) );
		}
	}

	return dst;

}

/* Retorna, separamente */
void findValueMat( cv::Mat& binary, std::vector<int> &idx, std::vector<int> &idy )
{
	
	const int M = binary.rows;
	const int N = binary.cols;
	float *t;
	
	for (int m = 0; m < M; ++m)
	{
		
		t = binary.ptr<float>(m);
		for (int n = 0; n < N; ++n)
		{
			
			if (t[n] > 0)
			{
				idx.push_back(n);
				idy.push_back(m);
			}
			
		}
		
	}
	
}

std::vector<int> findValueVectorInt(std::vector<int> &targetVector, int value)
{
	std::vector<int> result;
	size_t n;
	for (n = 0; n < targetVector.size(); ++n)
	{
		if (targetVector[n] == value)
		{
			result.push_back((int) n);
		}
	}
	return result;
	/*TODO*/
}

std::vector<int> subvectorInt(std::vector<int> targetVector, std::vector<int> chosenIndexes)
{
	std::vector<int> result;
	size_t i;
	for (i = 0; i < chosenIndexes.size(); ++i)
	{
		result.push_back( targetVector.at(chosenIndexes.at(i)) );
	}
	return result;
}

/* Debug conly, will be removed eventually: */
void printVectorInt(const char *vName, std::vector<int> v, int pauseMode, int indexes)
{

	std::cout << "--- Begin " << vName << " ---" << std::endl;
	for (size_t i = 0; i < v.size(); ++i)
	{

		if (indexes == INDEXES_ON)
		{
			std::cout << i + 1 << ": " << v[i];
		}
		else std::cout << v[i];

		if ((i + 1) % 6 == 0)
		{
			std::cout << std::endl;
		}
		else std::cout << " ";

	}
	std::cout << std::endl << "--- End " << vName << " ---" << std::endl;
	if (pauseMode == PAUSE_ON)
	{
		system("PAUSE");
	}
	else std::cout << "		No pause." << std::endl;

	return;

}

std::vector<Mat> createGaborFilters( double initialTheta, double increment, int size )
{
	
	/*
	* ksize Size of the filter returned.
	* sigma Standard deviation of the gaussian envelope.
	* theta Orientation of the normal to the parallel stripes of a Gabor function.
	* lambda Wavelength of the sinusoidal factor.
	* gamma Spatial aspect ratio.
	* psi Phase offset.
	* ktype Type of filter coefficients. It can be CV_32F or CV_64F.
	**/
	
	double sigma = 4.0;
	double lambd = 10.0;
	double gamma = 1; //0.5
	double psi = 0;// CV_PI*0.5;
	int ktype = CV_32F;
	Mat kernel;
	
	std::vector<Mat> result;
	
	Size ksize = Size(size, size);
	while (initialTheta - 180 < 0)
	{
		
		kernel = getGaborKernel(ksize, sigma, initialTheta, lambd, gamma, psi, ktype);
		result.push_back(kernel);
		initialTheta = initialTheta + increment;
		
	}
	return result;
	
}

void saveCoordinates( const std::string &filename, int coordX, int coordY )
{
	std::ofstream myfile;
	myfile.open(filename);
	myfile << coordX << std::endl;
	myfile << coordY;
	myfile.close();
}

vector<int> loadCoordinates( const std::string &filename )
{
	std::vector<int> result;
	std::ifstream myfile;
	std::string line;
	
	myfile.open(filename);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			/*
			** result.at(0) = X;
			** result.at(1) = Y;
			*/
			result.push_back((int) std::stol(line));
		}
		myfile.close();
	}
	
	return result;
}

void saveFeatures(const std::string &filename, vector<float> features)
{
	
	std::ofstream myfile;
	int k;
	
	myfile.open(filename);
	k = 0;
	while (k < features.size())
	{
		
		myfile << features.at(k) << std::endl;
		k++;
		
	}
	myfile.close();
	
}

vector<float> loadFeatures(const std::string &filename)
{
	
	std::vector<float> result;
	std::ifstream myfile;
	std::string line;
	
	myfile.open(filename);
	if (myfile.is_open())
	{
		
		while (getline(myfile, line))
		{
			
			result.push_back(std::stof(line));
			
		}
		myfile.close();
		
	}
	return result;
	
}

void salilNormalize(Mat ROI, float M0, float V0)
{
	int y;
	int x;

	Scalar ViVector;
	Scalar MiVector;

	float Mi;
	float Vi;
	
	/* Calcule mean and standard deviation of the block: */
	meanStdDev(ROI, MiVector, ViVector);
	Mi = (float) MiVector.val[0];
	Vi = (float) ViVector.val[0];
	
	for (y = 0; y < ROI.rows; y++)
	{
		for (x = 0; x < ROI.cols; x++)
		{
			/* I(x, y) > Mi */
			if (ROI.at<float>(y, x) - Mi > 0)
			{
				ROI.at<float>(y, x) = M0 + sqrt((V0 * pow(ROI.at<float>(y, x) - Mi, 2)) / Vi);
			}
			else ROI.at<float>(y, x) = M0 - sqrt((V0 * pow(ROI.at<float>(y, x) - Mi, 2)) / Vi);
		}
	}
}

float featureSector( Mat ROI, float M0, float V0, Mat filter )
{
	float result;
	float PiTheta;
	Mat dst;
	Scalar PiThetaScalar;
	int k;
	
	#ifdef DEBUG_FINGERCODES
		salilNormalize(ROI, M0, V0);
		ROI.copyTo(dst);
	#endif
	#ifndef DEBUG_FINGERCODES
		ROI.copyTo(dst);
		salilNormalize(dst, M0, V0);
	#endif
	filter2D(dst, dst, -1, filter, Point(-1, -1), 0, BORDER_DEFAULT);
	PiThetaScalar = cv::mean(dst);
	PiTheta = (float) PiThetaScalar.val[0];
	result = 0;
	
	/*
	** dst is Salil's F(i, theta):
	** PiTheta is Salil's P(i, theta)
	** dst.size().area() is Salil's n(i)
	** result is Salil's V(i, theta)
	*/
	for (k = 0; k < dst.size().area(); k++)
	{
		
		/* Definition of module function: */
		if ( dst.at<float>(k) - PiTheta > 0 ) { result += ( dst.at<float>(k) - PiTheta ); }
		else result += ( PiTheta - dst.at<float>(k) );
		
	}
	result = result / k; /* k = dst.size().area() */
	
	return result;
}

vector<float> featureVector( Mat ROI, int sectionWidth, int sectionHeight, float M0, float V0, vector<Mat>filterBank )
{
	vector<float> result;
	int horizontalLimit;
	int verticalLimit;
	int i;
	int j;
	int k;
	Rect subROIAux;
	Mat subROI;
	
	#ifdef DEBUG_FINGERCODES
		Mat visibleROI;
		ROI.convertTo(visibleROI, CV_8UC1);
	#endif

	horizontalLimit = ROI.cols / sectionWidth;
	verticalLimit = ROI.rows / sectionHeight;

	for (i = 0; i < verticalLimit; i++)
	{
		for (j = 0; j < horizontalLimit; j++)
		{
			subROIAux = Rect(j * sectionWidth, i * sectionHeight, sectionWidth, sectionHeight);
			subROI = ROI(subROIAux);
			for (k = 0; k < 8; k++)
			{
				result.push_back(featureSector(subROI, M0, V0, filterBank.at(k)));
				#ifdef DEBUG_FINGERCODES
					ROI.convertTo(visibleROI, CV_8UC1);
				#endif
			}
		}
	}

	/*if (true) { Mat outputSaver = normalizaImg(ROI);
	cv::imwrite("ResourceFiles/images/output/intermediate/Normalizada.jpg", outputSaver); }*/

	return result;
}

std::vector<Mat> displayFeatures(std::vector<float>vectorOfFeatures)
{
	int i;
	int theta;
	vector<Mat> result;
	Mat matIterator;
	
	for (theta = 0; theta < 8; theta++)
	{
		/*20 = (2 * SIZE_ROI) / Sector_Size = (2 * 160)/16*/
		matIterator = Mat::zeros(20, 20, CV_32F);
		for (i = 0; i < 400; i++) /* 400 = sectors_count */
		{
			matIterator.at<float>(i) = vectorOfFeatures.at(i * 8 + theta);
		}
		matIterator = normalizaImg(matIterator);
		matIterator.convertTo(matIterator, CV_8UC1);
		/*Sector_Size = 16*/
		resize(matIterator, matIterator, Size(matIterator.rows * 16, matIterator.cols * 16), 0, 0, INTER_NEAREST);
		result.push_back(matIterator);
	}
	
	return result;
}

std::vector<Mat> convertFeatures(std::vector<float>vectorOfFeatures)
{
	int i;
	int theta;
	vector<Mat> result;
	Mat matIterator;

	for (theta = 0; theta < 8; theta++)
	{
		matIterator = Mat::zeros(10, 8, CV_32F);
		for (i = 0; i < 80; i++)
		{
			matIterator.at<float>(i) = vectorOfFeatures.at(i * 8 + theta);
		}
		matIterator = normalizaImg(matIterator);
		result.push_back(matIterator);
	}

	return result;
}

float distanceFeatures(std::vector<Mat>input1, std::vector<Mat>input2)
{
	
	int iterator;
	float result;
	
	result = 0;

	if (input1.size() != input2.size())
	{
		std::cout << "Size mismatch." << std::endl;
		return -1;
	}
	
	for( iterator = 0; iterator < input1.size(); iterator++ )
	{
		result += ( float ) norm(input1.at(iterator), input2.at(iterator), NORM_L2);
	}
	
	return result;
	
}

void saveComparisons(const std::string &filename, std::vector<float> comparisonVector, float threshold)
{
	
	std::ofstream myfile;
	int k;
	float aux;
	int successCounter;
	
	myfile.open("ResourceFiles/comparisons/" + filename + std::to_string(threshold) + ".txt");
	successCounter = 0;
	for(k = 0; k < comparisonVector.size(); k++)
	{
		
		aux = comparisonVector.at(k);
		if ( aux - threshold > 0 )
		{
			myfile << aux << " False" << std::endl;
			continue;
		}
		myfile << comparisonVector.at(k) << " True" << std::endl;
		successCounter++;
		
	}
	myfile << successCounter;
	myfile.close();
	
}

void saveAllComparisons(const std::string &filename, std::vector<float> comparisonVector)
{
	int k;
	float threshold;

	threshold = 5000;
	for (k = 0; k < 13; k++)
	{
		threshold = threshold + 1000;
		saveComparisons(filename, comparisonVector, threshold);
	}
}

void saveFilteredRegions(cv::Mat imageRoi, std::vector<cv::Mat> filterBank)
{
	
	int counter;
	cv::Mat outputSaver;
	String str = "ResourceFiles/images/output/intermediate/Normalizada";
	Mat currentFilter;
	
	for (counter = 0; counter < filterBank.size(); counter++)
	{
		currentFilter = filterBank.at(counter);
		imageRoi.convertTo(imageRoi, CV_32F);
		filter2D(imageRoi, outputSaver, -1, filterBank.at(counter), Point(-1, -1), 0, BORDER_DEFAULT);
		outputSaver = normalizaImgDebug(outputSaver);
		//outputSaver.convertTo(outputSaver, CV_8UC1);
		cv::imwrite(str + std::to_string(counter) + ".jpg", outputSaver);
	}
	
}