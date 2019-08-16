#ifndef __STDDVIM
#define __STDDVIM

#include <opencv2/highgui/highgui.hpp>
/*function imagemstdv = stddvim(imagem, blksze);*/
cv::Mat stddvim(cv::Mat imagem, int blksze);

#endif