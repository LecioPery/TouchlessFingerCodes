#ifndef __INSEREZEROS
#define __INSEREZEROS

#include <opencv2/highgui/highgui.hpp>

using namespace cv;

/*function [zeroIm, numBlocosL, numBlocosC] = inserezeros(imOrig, blksze)*/
class inserezeros
{
	private:
	Mat zeroIm;
	int numBlocosL;
	int numBlocosC;

	public:
	inserezeros(Mat imOrig, int blksze);
	Mat getzeroIM();
	int getnumBlocosL();
	int getnumBlocosC();
};

#endif