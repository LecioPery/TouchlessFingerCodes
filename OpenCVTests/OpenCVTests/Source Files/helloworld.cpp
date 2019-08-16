#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <cstdlib>

#include "..\Header Files\helloWorld.hpp"

using namespace cv;
using namespace std;

/*
*  author : v3nw
*  desc   : OpenCV "Hello World" Sample
*  date   : 30.05.2013
*  web    : v3nw.com
*  source : http://opencv.willowgarage.com/wiki/Getting_started
*/

void helloWorld()
{

	char* windowName = "HelloWorldWindow";
	cvNamedWindow(windowName, 0x01);
	IplImage* img = cvCreateImage(cvSize(0x96, 0x32), IPL_DEPTH_8U, 0x01);
	CvFont font;
	double hScale = 0x01;
	double vScale = 0x01;
	int lineWidth = 0x01;
	cvInitFont(&font, CV_FONT_HERSHEY_COMPLEX_SMALL | CV_AA, hScale, vScale, 0x00, lineWidth);
	cvPutText(img, "Hello World.", cvPoint(0x00, 0x1E), &font, cvScalar(0x00, 0x00, 0x00));
	cvShowImage(windowName, img);
	cvWaitKey();
	return;

}

void imshowTest( char *name )
{

	Mat image;
	image = imread( name, CV_LOAD_IMAGE_UNCHANGED );   // Read the file

	if (!image.data)                              // Check for invalid input
	{
		cout << "Could not open or find the image" << std::endl;
		return;
	}

	namedWindow( "Display window", WINDOW_AUTOSIZE );// Create a window for display.
	imshow( "Display window", image );                   // Show our image inside it.

	waitKey(0);
	return;

}