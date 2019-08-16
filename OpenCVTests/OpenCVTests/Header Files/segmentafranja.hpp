#ifndef __SEGMENTAFRANJA
#define __SEGMENTAFRANJA

#include <opencv2/highgui/highgui.hpp>

/* Isso ainda tem de ser repensado: */
//#define MAGIC_NUMBER 99 /*Unused*/

using namespace cv;

/*function[normim, mask, maskind] = segmentafranja(im, blksze)*/
class segmentafranja
{

	Mat normim;
	Mat mask;
	//int maskind[MAGIC_NUMBER];/*Unused*/
	
	public:
	segmentafranja(Mat im, int blksze);
	Mat getNormim();
	Mat getMask();
	//int *getMaskind(); /*Unused*/

};

#endif