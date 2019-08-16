#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>

#include "..\Header Files\num2strSpecial.hpp"
#include "..\Header Files\preparaim.hpp"
#include "..\Header Files\extraidirecional.hpp"
#include "..\Header Files\trickyLib.hpp"
#include "..\Header Files\findCenterRoi.hpp"
#include "..\Header Files\generateCenters.hpp"
#include "..\Header Files\commonConstants.hpp"
#include "..\Header Files\Falange\falange.hpp"

/* MATLAB constants: */
/*#define NIMG 25*/

#define PREFIX "SIRE-1_"
#define SUFIX "_0"
#define EXTENSION ".jpg"

/* MATLAB Boundaries: */
#define ORIENT_VEC_WIDTH 99

using namespace cv;
using namespace std;

std::vector<std::string> initializeAllowedList( void );
bool insideAllowedList(int person, int finger, std::vector<std::string> dictionary);

void finalPicture(Mat flow, Mat original);

/* Note: precision issues in regard of Matlab: */
int main(int argc, char **argv)
{
	
	//http://funvision.blogspot.com/2016/02/opencv-31-tutorial-optical-flow.html
	Mat predecessor;
	Mat sucessor;
	int k, i, j;
	String address = IMG_INPUT_PATH;
	//UMat uflow;
	Mat uflow;
	Mat carbonCopy;
	
	k = 1;
	i = 1;
	j = 1;
	predecessor = imread(address + "rotated/" + to_string(k) + "_" + to_string(i) + "_" + to_string(j) + SUFIX + EXTENSION, CV_8UC1);
	sucessor = imread(address + "rotated/" + to_string(k) + "_" + to_string(i) + "_" + to_string(j + 1) + SUFIX + EXTENSION, CV_8UC1);
	predecessor.convertTo(uflow, CV_32FC2, 1.0, 0);
	
	predecessor.copyTo(carbonCopy);
	carbonCopy.setTo(Scalar::all(0));
	sucessor.copyTo(carbonCopy(Rect(predecessor.cols - sucessor.cols, predecessor.rows - sucessor.rows, sucessor.cols, sucessor.rows)));
	
	calcOpticalFlowFarneback(predecessor, carbonCopy, uflow, 0.5, 3, 15, 3, 5, 1.2, 0);
	finalPicture(uflow, carbonCopy);
	cout << "Fim." << endl;
	
	return 0;
	
}

void finalPicture(Mat flow, Mat original)
{
	for (int y = 0; y < original.rows; y += 5) {
		for (int x = 0; x < original.cols; x += 5)
		{
			// get the flow from y, x position * 10 for better visibility
			const Point2f flowatxy = flow.at<Point2f>(y, x) * 10;
			// draw line at flow direction
			line(original, Point(x, y), Point(cvRound(x + flowatxy.x), cvRound(y + flowatxy.y)), Scalar(255, 0, 0));
			// draw initial point
			circle(original, Point(x, y), 1, Scalar(0, 0, 0), -1);


		}
	}
}

std::vector<std::string> initializeAllowedList( void )
{
	
	/**
	 * Before going, consider:
	 *	1) Problem on follow images for a region of 320 x 320:
	 *		subject-finger-sample
	 *		5-1-2
	 *		14-5-1
	 *		14-5-4
	 *		14-5-5
	*/

	std::vector<std::string> dictionary;
	
	dictionary.push_back("1_1");
	/*dictionary.push_back("1_2");
	dictionary.push_back("1_3");
	dictionary.push_back("1_4");
	dictionary.push_back("1_5");
	dictionary.push_back("1_7");
	dictionary.push_back("1_8");
	dictionary.push_back("1_9");
	dictionary.push_back("1_10");
	dictionary.push_back("2_1");
	dictionary.push_back("2_3");
	dictionary.push_back("2_4");
	dictionary.push_back("2_8");
	dictionary.push_back("3_1");
	dictionary.push_back("3_4");
	dictionary.push_back("3_5");
	dictionary.push_back("4_6");
	dictionary.push_back("5_1");
	dictionary.push_back("5_2");
	dictionary.push_back("5_4");
	dictionary.push_back("5_5");
	dictionary.push_back("5_7");
	dictionary.push_back("5_8");
	dictionary.push_back("6_1");
	dictionary.push_back("6_2");
	dictionary.push_back("6_6");
	dictionary.push_back("7_1");
	dictionary.push_back("7_3");
	dictionary.push_back("7_4");
	dictionary.push_back("7_7");
	dictionary.push_back("7_8");
	dictionary.push_back("7_10");
	dictionary.push_back("8_1");
	dictionary.push_back("8_2");
	dictionary.push_back("8_3");
	dictionary.push_back("8_4");
	dictionary.push_back("8_8");
	dictionary.push_back("8_10");
	dictionary.push_back("9_5");
	dictionary.push_back("9_7");
	dictionary.push_back("9_10");
	dictionary.push_back("10_1");
	dictionary.push_back("10_3");
	dictionary.push_back("10_4");
	dictionary.push_back("10_7");
	dictionary.push_back("10_10");
	dictionary.push_back("11_1");
	dictionary.push_back("11_3");
	dictionary.push_back("11_6");
	dictionary.push_back("11_8");
	dictionary.push_back("11_9");
	dictionary.push_back("12_5");
	dictionary.push_back("12_8");
	dictionary.push_back("13_3");
	dictionary.push_back("13_4");
	dictionary.push_back("13_5");
	dictionary.push_back("13_7");
	dictionary.push_back("13_8");
	dictionary.push_back("14_4");
	//dictionary.push_back("14_5");
	dictionary.push_back("14_7");
	dictionary.push_back("14_8");
	dictionary.push_back("15_5");
	dictionary.push_back("15_8");
	dictionary.push_back("16_7");
	dictionary.push_back("16_8");
	dictionary.push_back("16_10");
	dictionary.push_back("17_6");
	dictionary.push_back("18_2");
	dictionary.push_back("18_6");
	dictionary.push_back("18_7");
	dictionary.push_back("18_8");
	dictionary.push_back("19_3");
	dictionary.push_back("19_4");
	dictionary.push_back("19_5");
	dictionary.push_back("19_7");
	dictionary.push_back("19_8");
	dictionary.push_back("19_9");
	dictionary.push_back("20_1");
	dictionary.push_back("20_2");
	dictionary.push_back("20_5");
	dictionary.push_back("20_7");*/
	
	return dictionary;

}

bool insideAllowedList(int person, int finger, std::vector<std::string> dictionary)
{
	std::string entry(to_string(person) + "_" + to_string(finger));
	int i;
	for (i = 0; i < dictionary.size(); i++)
	{
		if (!dictionary.at(i).compare(entry)) return true;
	}
	return false;
}

void fingerCodes()
{
	/* Other Matlab variables */
	Mat imageTarget;
	Mat imageTarget2;
	String str = IMG_INPUT_PATH;
	String str2 = IMG_OUTPUT_PATH;
	String strW = IMG_INPUT_PATH_W;
	String str2W = IMG_OUTPUT_PATH_W;

	/* MIDAS(RIP) output: */
	//HT12 = create_h12("HT1.bmp", "HT2.bmp");

	/* Gabor filterbank: */
	vector<Mat> filterBank;

	/* Salil: */
	Mat normalizedROI;
	vector<float> vectorOfFeatures;

	vector<findCenterRoi> images;
	Rect regionOfInterest;
	Mat imageRoi;

	vector<int> x;
	vector<int> y;
	vector<int> coordinates;
	/* Old but gold:
	int highest;
	*/

	vector<Mat>features;
	int k;
	int i;
	int j;
	int l;

	/* Comparisons: */
	vector<vector<Mat>> fingerDatabase;
	vector<vector<Mat>> pickedFingerDatabase;
	vector<float> verticalComparisons;
	vector<float> horizontalComparisons;

	/* filterBank initialization: */
	filterBank = createGaborFilters(0, 22.5, 33);

	/* Center generation, uses MIDAS internally: */
	//centerGenerate( (char *) (str + "salum/").c_str(),  (char *) CENTERS_PATH);

	/************ Luan Caius test ***********
	//int h = falange( (str + "salum/1_1_1_0.jpg").c_str());
	int h = 33;
	imageTarget = imread(str + "rotated/1_1_1_0.jpg", CV_8UC1);
	regionOfInterest = Rect(0, h, imageTarget.cols, 20);
	imageTarget(regionOfInterest) = 0;*/

	std::vector<std::string> dictionary;

	dictionary = initializeAllowedList();

	for (k = 1; k <= NPEOPLE; k++)
	{

		//system( ("rmdir /s /q " + str2W + "centers\\" + std::to_string(k) ).c_str());
		system(("rmdir /s /q " + str2W + "fingerCodes\\" + std::to_string(k)).c_str());
		//system( ("mkdir " + str2W + "centers\\" + std::to_string(k) ).c_str() );
		system(("mkdir " + str2W + "fingerCodes\\" + std::to_string(k)).c_str());
		for (i = 1; i <= NFINGER; i++)
		{

			if (!insideAllowedList(k, i, dictionary)) { std::cout << "Skipped [k][i]: [" << k << "] [" << i << "]" << std::endl; continue; }
			//system( ("rmdir /s /q " + str2W + "centers\\" + std::to_string(k) + "\\" + std::to_string(i) ).c_str());
			system(("rmdir /s /q " + str2W + "fingerCodess\\" + std::to_string(k) + "\\" + std::to_string(i)).c_str());
			//system( ("mkdir " + str2W + "centers\\" + std::to_string(k) + "\\" + std::to_string(i) ).c_str());
			system(("mkdir " + str2W + "fingerCodes\\" + std::to_string(k) + "\\" + std::to_string(i)).c_str());
			for (j = 1; j <= NSAMPLE; j++)
			{

				system(("rmdir /s /q " + str2W + "fingerCodess\\" + std::to_string(k) + "\\" + std::to_string(i) + "\\" + std::to_string(j)).c_str());
				system(("mkdir " + str2W + "fingerCodes\\" + std::to_string(k) + "\\" + std::to_string(i) + "\\" + std::to_string(j)).c_str());

				imageTarget = imread(str + "rotated/" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + SUFIX + ".jpg", CV_8UC1);
				imageTarget2 = imread(str + "salum/" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + SUFIX + ".jpg", CV_8UC1);
				coordinates = loadCoordinates(CENTERS_PATH + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + ".txt");
				if (coordinates.at(0))
				{

					if ((coordinates.at(0) - SIZE_ROI >= 0) && (coordinates.at(0) + SIZE_ROI <= imageTarget.cols) && (coordinates.at(1) - SIZE_ROI >= 0) && (coordinates.at(1) + SIZE_ROI <= imageTarget.rows))
					{

						regionOfInterest = Rect(coordinates.at(0) - SIZE_ROI, coordinates.at(1) - SIZE_ROI, 2 * SIZE_ROI, 2 * SIZE_ROI);
						imageRoi = imageTarget(regionOfInterest);
						imageRoi.convertTo(imageRoi, CV_32F);
						vectorOfFeatures = featureVector(imageRoi, 16, 16, 100.0, 100.0, filterBank);

						features = displayFeatures(vectorOfFeatures);
						for (l = 0; l < features.size(); l++)
						{
							imwrite(str2 + "fingerCodes/" + std::to_string(k) + "/" + std::to_string(i) + "/" + std::to_string(j) + "/" + num2strSpecial(l + 1) + EXTENSION, features.at(l));
						}
						saveFeatures(str2 + "fingerCodes/" + std::to_string(k) + "/" + std::to_string(i) + "/" + std::to_string(j) + "/feature.txt", vectorOfFeatures);
						saveFeatures(str2 + "matcher/" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + ".txt", vectorOfFeatures);

						imageRoi.convertTo(imageRoi, CV_8UC1);
						features = convertFeatures(vectorOfFeatures);
						fingerDatabase.push_back(features);
						std::cout << "k: " << k << " i: " << i << " j: " << j << std::endl;
						imwrite(str2 + "centers/" + std::to_string(k) + "/" + std::to_string(i) + "/" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + SUFIX + "orig" + EXTENSION, imageTarget(regionOfInterest));
						imwrite(str2 + "centers/" + std::to_string(k) + "/" + std::to_string(i) + "/" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + SUFIX + "salum" + EXTENSION, imageTarget2(regionOfInterest));
						imageTarget(regionOfInterest) = 0;
						imageTarget2(regionOfInterest) = 0;
						imwrite(str2 + "centers/" + std::to_string(k) + "/" + std::to_string(i) + "/" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + SUFIX + EXTENSION, imageTarget);
						continue;

					}
					else
					{
						std::cout << "Boundaries exceded at image [k][i][j] = [" << k << "][" << i << "]" << "[" << j << "]" << std::endl;
						std::cout << "(coordinates.at(0) - SIZE_ROI > 0): " << coordinates.at(0) - SIZE_ROI << " (coordinates.at(0) + SIZE_ROI < imageTarget.cols): " << coordinates.at(0) + SIZE_ROI << std::endl;
						std::cout << "(coordinates.at(1) - SIZE_ROI > 0): " << coordinates.at(1) - SIZE_ROI << " (coordinates.at(1) + SIZE_ROI < imageTarget.rows): " << coordinates.at(1) + SIZE_ROI << std::endl;
					}
				}
				else
				{
					std::cout << "No center found." << std::endl;
				}

				/* Hack so that centers not found have at least a representation, even if fake. */
				fingerDatabase.push_back(features);
				std::cout << "No center to show at image [k][i][j] = [" << k << "][" << i << "]" << "[" << j << "]" << std::endl;

			}

		}

	}
	
	imageRoi.convertTo(imageRoi, CV_8UC1);
	saveFilteredRegions(imageRoi, filterBank);
	
	std::cout << "Fim." << std::endl;
	
}