#ifndef __EXTRAIDIRECIONAL
#define __EXTRAIDIRECIONAL

#include <opencv2/highgui/highgui.hpp>

using namespace cv;

class extraidirecional
{

	Mat orientim;
	Mat mask;

	public:
	/*function [orientim, mask] = extraidirecional(fingerprint, blksze);*/
	extraidirecional(Mat fingerprint, int blksze);
	Mat getMask(void);
	Mat extraidirecional::getOrientim(void);

};

#endif