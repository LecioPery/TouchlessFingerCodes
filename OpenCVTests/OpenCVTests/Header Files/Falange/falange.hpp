/*
 * falange.h
 *
 *  Created on: 19/09/2011
 *      Author: Caius
 */

#ifndef FALANGE_H_
#define FALANGE_H_

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <imgproc.hpp>
//#include "useful_macros.h"	  // opencv macros include file

#define CV_CLAHE_RANGE_FULL 0
#define CV_CLAHE_RANGE_INPUT 1

extern int falange(const char *str1);
extern int phalange_QT(IplImage *original);

#endif /* FALANGE_H_ */
