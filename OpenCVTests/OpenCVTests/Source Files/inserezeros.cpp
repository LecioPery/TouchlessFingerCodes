#include <cmath>
#include <iostream>
#include <opencv2/highgui/highgui.hpp>

#include "..\Header Files\inserezeros.hpp"

using namespace cv;

/* Incompleto */
inserezeros::inserezeros(Mat imOrig, int blksze)
{
	int C;
	int L;
	int numL;
	int numC;

	double min;
	double max;

	C = imOrig.cols;
	L = imOrig.rows;
	
	numBlocosL = ceil(L / blksze);
	numBlocosC = ceil(C / blksze);

	numL = blksze*numBlocosL;
	numC = blksze*numBlocosC;

	cv::minMaxLoc(imOrig, &min, &max);

	/* zeroIm = max(max(imOrig))*ones(numL, numC); */
	zeroIm = ( (float) max ) * Mat::ones(numL, numC, CV_32FC3);

	/* zeroIm(1:L, 1:C) = imOrig; */
	imOrig.copyTo(zeroIm);

	/*std::cout << "--------i---------" << std::endl;
	std::cout << "C: " << C << " L: " << L << std::endl;
	std::cout << "numC: " << numC << " numL: " << numL << std::endl;
	std::cout << "numBlocosC: " << numBlocosC << " numBlocosL: " << numBlocosL << std::endl;
	std::cout << "-----------------" << std::endl;*/
}

Mat inserezeros::getzeroIM()
{
	return zeroIm;
}

int inserezeros::getnumBlocosL()
{
	return numBlocosL;
}

int inserezeros::getnumBlocosC()
{
	return numBlocosC;
}

/**
Matlab equivalent code:
function [zeroIm, numBlocosL, numBlocosC] = inserezeros(imOrig, blksze)

% Funcao que, dadas a imagem original imOrig e a dimensao do sub-bloco
% blksze, insere zeros para que imagem contenha um numero inteiro de
% sub-blocos.

% Determina as dimensoes de imOrig
[L, C] = size(imOrig);

% Determina num inteiro de sub-blocos para a imagem
numBlocosL = ceil(L/blksze);
numBlocosC = ceil(C/blksze);

% Determina a nova dimensao para a imagem, baseando-se no numero de
% sub-bloco
numL = blksze*numBlocosL;
numC = blksze*numBlocosC;

% Cria uma imagem com a dimensao e a cor de fundo necessarios
zeroIm = max(max(imOrig))*ones(numL, numC);

% Insere a imOrig dentro de zeroIm
zeroIm(1:L, 1:C) = imOrig;

return


*/