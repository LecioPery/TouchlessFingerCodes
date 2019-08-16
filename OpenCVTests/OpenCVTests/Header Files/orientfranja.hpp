#ifndef __ORIENTFRANJA
#define __ORIENTFRANJA

#include <opencv2/highgui/highgui.hpp>

using namespace cv;

class orientfranja
{

	private:
	Mat orientim;

	static Mat f1Mat;
	static Mat f1xMat;
	static Mat f1yMat;
	static Mat f2Mat;
	static Mat f3Mat;
	
	public:
	/*function [orientim] = orientfranja(im)*/
	Mat getOrientim(void);
	orientfranja(Mat normim);

};

#endif