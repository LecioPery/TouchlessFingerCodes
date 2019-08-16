/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Bruno Macchiavello  bruno@cic.unb.br             %
%          Alexandre Zaghetto  alexandre@cic.unb.br         %
%          Mamede Lima-Marques limamarques@gmail.com        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 30/11/2010					                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UnB  - Universidade de Brasília                           %
% IE   - Institute of Exact Sciences                        %
% CIC  - Department of Computer Science                     %
% LISA - Laboratory of Image, Signal and Audio              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPAI - Centre for Research on Architecture of Information %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

//#include <opencv2/highgui/highgui.hpp>

//using namespace cv;

#define MaxFiltSize 7

cv::Mat create_h12(char *str1, char *str2);
int powerlaw(IplImage *, IplImage * );
int imdilate(IplImage *, IplImage *);
int cannythreshold(float *, float *, int *, IplImage *);
