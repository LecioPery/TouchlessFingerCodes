#ifndef __TRICKY_LIB
#define __TRICKY_LIB

	#include <opencv2/highgui/highgui.hpp>
	#include <opencv2/opencv.hpp>
	#include <opencv/cv.h>
	#include <iostream>
	#include <fstream>

	#define	METHOD_SUM 0
	#define	METHOD_SUBTRACT 1
	#define METHOD_PRODUCT 2
	#define METHOD_DIVIDE 3
	#define PAUSE_ON 1
	#define PAUSE_OFF 0
	#define INDEXES_ON 1
	#define INDEXES_OFF 0

	void matParse( cv::Mat &Target, const char *varname );
	cv::Mat setMatValuesFloat( cv::Mat img, cv::Range heightInterval, cv::Range widthInterval, float value );
	cv::string getImageType(int number);
	float mean(float *data, int n);
	float standard_deviation(float *data, int n);
	float std2(cv::Mat imagem);
	/*cv::Mat extractIndexes(cv::Mat target, int *lineIndex, int hLenght, int *columnIndex, int cLenght);*/
	cv::Mat thresholdImage(cv::Mat img, float thresh);
	cv::Mat normalizaImgDebug(const cv::Mat&);
	float standard_deviation_vector(cv::vector<float>);
	cv::Mat bytewiseUpdateImage(cv::Mat target, float value, int METHOD);
	void matReadDouble(const char *file, std::vector<double>& v, const char *varName);
	void matReadMat(const char *file, cv::Mat& v, const char *varName);
	void printVectorDouble(const char *vName, std::vector<double> v, int pauseMode, int indexes);
	cv::Mat matProduct(cv::Mat input1, cv::Mat input2, const int METHOD);
	void matSinCos(cv::Mat Gxx, cv::Mat Gxy, cv::Mat Gyy, cv::Mat *sinThetha, cv::Mat *cosThetha);
	cv::Mat matAtan2(cv::Mat input1, cv::Mat input2);
	cv::Mat matFunction(cv::Mat target, const char *function);
	cv::Mat gradientMat(cv::Mat TargetImage);
	cv::Mat tailMat(cv::Mat TargetImage);
	cv::Mat normalizaImg(const cv::Mat& src);
	void findValueMat( cv::Mat& binary, std::vector<int> &idx, std::vector<int> &idy );
	std::vector<int> findValueVectorInt(std::vector<int> &targetVector, int value);
	std::vector<int> subvectorInt(std::vector<int> targetVector, std::vector<int> chosenIndexes);
	void printVectorInt(const char *vName, std::vector<int> v, int pauseMode, int indexes);
	std::vector<cv::Mat> createGaborFilters(double initialTheta, double increment, int size);
	void saveCoordinates(const std::string &filename, int coordX, int coordY);
	std::vector<int> loadCoordinates(const std::string &filename);
	void saveFeatures(const std::string &filename, std::vector<float> features);
	std::vector<float> loadFeatures(const std::string &filename);
	void salilNormalize(cv::Mat ROI, float M0, float V0);
	float featureSector(cv::Mat ROI, float M0, float V0, cv::Mat filter);
	std::vector<float> featureVector(cv::Mat ROI, int sectionWidth, int sectionHeight, float M0, float V0, std::vector<cv::Mat>filterBank);
	std::vector<cv::Mat> displayFeatures(std::vector<float>vectorOfFeatures);
	std::vector<cv::Mat> convertFeatures(std::vector<float>vectorOfFeatures);
	float distanceFeatures(std::vector<cv::Mat>input1, std::vector<cv::Mat>input2);
	void saveAllComparisons(const std::string &filename, std::vector<float> comparisonVector);
	void saveFilteredRegions(cv::Mat imageRoi, std::vector<cv::Mat> filterBank);
#endif