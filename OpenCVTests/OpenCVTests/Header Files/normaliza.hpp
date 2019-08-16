#ifndef __NORMALIZA
#define __NORMALIZA

#include <opencv2/highgui/highgui.hpp>
using namespace cv;

/*function im = normaliza(im)*/
class normaliza
{
	private:
	Mat im;

	public:
	normaliza::normaliza(Mat im);
	Mat getIm(void);
};

#endif