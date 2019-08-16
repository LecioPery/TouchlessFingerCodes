#include "..\Header Files\orientfranja.hpp"
#include "..\Header Files\trickyLib.hpp"

#include <math.h>
#include <vector>

#include <opencv2/highgui/highgui.hpp>

//#define DEBUG_CODE1 1
//#define DEBUG_CODE2 1

#ifdef DEBUG_CODE1
#include <iostream>
#endif

#ifdef DEBUG_CODE2
#include <iostream>
#endif

using namespace cv;

orientfranja::orientfranja(Mat normim)
{
	
	/*std::vector<double> f1;
	std::vector<double> f1x;
	std::vector<double> f1y;
	std::vector<double> f2;
	std::vector<double> f3;*/

	Mat Gx;
	Mat Gy;
	Mat Gxx;
	Mat Gxy;
	Mat Gyy;
	Mat ImageF1;
	
	static Mat f1Mat;
	static Mat f1xMat;
	static Mat f1yMat;
	static Mat f2Mat;
	static Mat f3Mat;
	
	Mat sinTheta;
	Mat cosTheta;
	
	if (f1Mat.empty())
	{
		
		matReadMat("orientfiltros.mat", f1Mat, "f1");
		matReadMat("orientfiltros.mat", f1xMat, "f1x");
		matReadMat("orientfiltros.mat", f1yMat, "f1y");
		matReadMat("orientfiltros.mat", f2Mat, "f2");
		matReadMat("orientfiltros.mat", f3Mat, "f3");
		
	}
	
	/* Input = coder.load ('orientfiltros.mat', 'f1', 'f1x', 'f1y', 'f2', 'f3'); */
	

	/* Debug Purposes:*/
	#ifdef DEBUG_CODE1
	matReadDouble("orientfiltros.mat", f1, "f1");
	matReadDouble("orientfiltros.mat", f1x, "f1x");
	matReadDouble("orientfiltros.mat", f1y, "f1y");
	matReadDouble("orientfiltros.mat", f2, "f2");
	matReadDouble("orientfiltros.mat", f3, "f3");
	/*As funções acima podem ser melhoradas: trazendo matrizes completas ao invés de vetores*/

	printVectorDouble("f1", f1,  PAUSE_OFF, INDEXES_ON);
	printVectorDouble("f1x", f1x, PAUSE_OFF, INDEXES_ON);
	printVectorDouble("f1y", f1y, PAUSE_ON, INDEXES_ON);
	printVectorDouble("f2", f2, PAUSE_ON, INDEXES_OFF);
	printVectorDouble("f3", f3, PAUSE_ON, INDEXES_OFF);
	#endif
	
	cv::filter2D(normim, Gx, -1, f1xMat, cv::Point(-1, -1), 0, cv::BORDER_CONSTANT);
	cv::filter2D(normim, Gy, -1, f1yMat, cv::Point(-1, -1), 0, cv::BORDER_CONSTANT);

	Gxx = matProduct( Gx, Gx, METHOD_PRODUCT );
	Gxy = matProduct( Gx, Gy, METHOD_PRODUCT );
	Gyy = matProduct( Gy, Gy, METHOD_PRODUCT );

	#ifdef DEBUG_CODE2
	namedWindow("normim_before", CV_WINDOW_AUTOSIZE);
	namedWindow("Gx", CV_WINDOW_AUTOSIZE);
	namedWindow("Gy", CV_WINDOW_AUTOSIZE);
	namedWindow("Gxx", CV_WINDOW_AUTOSIZE);
	namedWindow("Gxy", CV_WINDOW_AUTOSIZE);
	namedWindow("Gyy", CV_WINDOW_AUTOSIZE);
	namedWindow("sinTheta", CV_WINDOW_AUTOSIZE);
	namedWindow("cosTheta", CV_WINDOW_AUTOSIZE);
	namedWindow("orientim", CV_WINDOW_AUTOSIZE);
	#endif
	
	cv::filter2D(Gxx, Gxx, -1, f2Mat, cv::Point(-1, -1), 0, cv::BORDER_CONSTANT);
	cv::filter2D(Gxy, Gxy, -1, f2Mat, cv::Point(-1, -1), 0, cv::BORDER_CONSTANT);
	cv::filter2D(Gyy, Gyy, -1, f2Mat, cv::Point(-1, -1), 0, cv::BORDER_CONSTANT);
	/* Gxy = 2 * Gxy: */
	Gxy = bytewiseUpdateImage( Gxy, 2, METHOD_PRODUCT );

	matSinCos(Gxx, Gxy, Gyy, &sinTheta, &cosTheta);
	cv::filter2D(sinTheta, sinTheta, -1, f3Mat, cv::Point(-1, -1), 0, cv::BORDER_CONSTANT);
	cv::filter2D(cosTheta, cosTheta, -1, f3Mat, cv::Point(-1, -1), 0, cv::BORDER_CONSTANT);

	/* Return Theta/2 + Pi/2: */
	this->orientim = matAtan2(sinTheta, cosTheta);

	#ifdef DEBUG_CODE2
	imshow("normim_before", normim);
	imshow("Gx", Gx);
	imshow("Gy", Gy);
	imshow("Gxx", Gxx);
	imshow("Gxy", Gxy);
	imshow("Gyy", Gyy);
	imshow("sinTheta", sinTheta);
	imshow("cosTheta", cosTheta);
	imshow("orientim", orientim);
	#endif

}

Mat orientfranja::getOrientim()
{
	return orientim;
}



/**
Matlab equivalent code:
function [orientim] = orientfranja(im)

% Carrega os filtros gaussianos Input.f1 (e gradientes Input.f1x e Input.f1y), Input.f2 e Input.f3;
Input = coder.load ('orientfiltros.mat', 'f1', 'f1x', 'f1y', 'f2', 'f3');

% Filtra a imagem com o gradiente do filtro Input.f1 em x, ou seja, Input.f1x
Gx = filter2(Input.f1x, im, 'same'); % Gradient da imagem em x

% Filtra a imagem com o gradiente do filtro Input.f1 em y, ou seja, Input.f1y
Gy = filter2(Input.f1y, im, 'same'); % Gradient da imagem em y

% Estimate the local orientation of each block by finding the axis
% orientation that minimises the area moment.

Gxx = Gx.^2;       % Momentos
Gxy = Gx.*Gy;
Gyy = Gy.^2;

% Filtra Gxx, Gxy e Gyy com o filtro gaussiano 2
Gxx = filter2(Input.f2, Gxx);  % Momentos suavizados
Gxy = 2*filter2(Input.f2, Gxy);
Gyy = filter2(Input.f2, Gyy);

% Determina sin2theta e cos2theta
denom = sqrt(Gxy.^2 + (Gxx - Gyy).^2);
sin2theta = Gxy./denom;
cos2theta = (Gxx-Gyy)./denom;

% Filtra sin2theta e cos2theta com o filtro gaussiano 3
cos2theta = filter2(Input.f3, cos2theta); % Senos e cossenos do dobro do angulo suavizados
sin2theta = filter2(Input.f3, sin2theta);

% Determina a orientaçao a partir de sin2thera e cos2theta
orientim = pi/2 + atan2(sin2theta,cos2theta)/2;

return
*/