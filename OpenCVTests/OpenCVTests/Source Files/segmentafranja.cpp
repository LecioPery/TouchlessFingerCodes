#include <iostream>
#include <fstream>

#include "..\Header Files\segmentafranja.hpp"
#include "..\Header Files\normaliza.hpp"
#include "..\Header Files\stddvim.hpp"
#include "..\Header Files\trickyLib.hpp"

//#define DEBUG_CODE3

float evaluateMean(cv::Mat mask, cv::Mat im);
float evaluateDeviation(cv::Mat mask, cv::Mat im);

segmentafranja::segmentafranja(Mat im, int blksze)
{

	float thresh;
	float mean;
	float deviation;
	
	Mat desvpadrim;

	thresh = (float) 0.1;
	im = normaliza(im).getIm(); /* Pronto. */
	desvpadrim = stddvim(im, blksze); /* Pronto. */
	
	/** Exibe na tela: */
	#ifdef DEBUG_CODE3
	namedWindow("desvpadrim", CV_WINDOW_AUTOSIZE);
	imshow("desvpadrim", desvpadrim);
	#endif
	
	/* Máscara */
	mask = thresholdImage(desvpadrim, thresh);
	/* Talvez isso nunca mais seja necessário. Muito provavelmente, pois chequei "no braço":
	namedWindow("mask2", CV_WINDOW_AUTOSIZE);
	imshow("mask2", mask);*/

	mean = evaluateMean(mask, im);
	im = bytewiseUpdateImage(im, mean, METHOD_SUBTRACT);
	deviation = evaluateDeviation(mask, im);
	normim = bytewiseUpdateImage(im, deviation, METHOD_DIVIDE);

	/*std::cout << "Mean: " << mean << std::endl;
	std::cout << "deviation: " << deviation << std::endl;*/

}

/*
 * Public:
*/
Mat segmentafranja::getNormim()
{
	return normim;
}

Mat segmentafranja::getMask()
{
	return mask;
}

/* Unused
int *segmentafranja::getMaskind()
{
	return maskind;
}*/

/*
* Private:
*/
float evaluateMean( cv::Mat mask, cv::Mat im )
{

	float *t;
	float *t2;

	double bufferAccumulator;
	
	int ammountOfValues;
	int nRows;
	int nCols;
	int r;
	int s;

	nRows = mask.rows;
	nCols = mask.cols;

	bufferAccumulator = 0;
	ammountOfValues = 0;
	for (r = 0; r < nRows; ++r)
	{

		t = mask.ptr<float>(r);
		t2 = im.ptr<float>(r);
		for (s = 0; s < nCols; ++s)
		{

			if (t[s] == 255)
			{

				bufferAccumulator = bufferAccumulator + t2[s];
				ammountOfValues++;
				
			}

		}

	}
	return (float) bufferAccumulator / ammountOfValues;
}

float evaluateDeviation(cv::Mat mask, cv::Mat im)
{

	float mean;
	float *t;
	float *t2;

	double bufferAccumulator;

	int ammountOfValues;
	int nRows;
	int nCols;
	int r;
	int s;

	nRows = mask.rows;
	nCols = mask.cols;

	bufferAccumulator = 0;
	ammountOfValues = 0;
	mean = evaluateMean( mask, im ); /* TEM que ser recalculado. */
	for (r = 0; r < nRows; ++r)
	{

		t = mask.ptr<float>(r);
		t2 = im.ptr<float>(r);
		for (s = 0; s < nCols; ++s)
		{

			/*if (t[s] == 255) std::cout << "r: " << r << " s: " << s << std::endl;*/
			if (t[s] == 255)
			{

				bufferAccumulator = bufferAccumulator + (t2[s] - mean) * (t2[s] - mean);
				ammountOfValues++;
				
			}

		}

	}
	return (float) sqrt( bufferAccumulator / (ammountOfValues - 1) );
}

/* As duas funções acima são similares ao "mean" de TrickyLib.cpp. Mas nem todos os pixels entram no cálculo. */

/**
Matlab equivalent code:
function [normim, mask, maskind] = segmentafranja(im, blksze)

% Determina o limiar de corte para geraçao da mascara
thresh = 0.1;

% Normaliza a imagem para que tenha media 0 e desvio padrao 1
im = normaliza(im);

% Gera uma imagem onde cada elemento de um sub-bloco de dimensoes blksze eh substituido
% pelo desvio padrao do sub-bloco
desvpadrim = stddvim(im, blksze);

% Gera uma mascara para identificar as regioes onde o desvio padrao
% local e inferior a um determinado limiar
mask = desvpadrim  > thresh;

% Identifica os indices onde a mascara é 1
maskind = find(mask);

% Com o auxilio da mascara, renormaliza apenas as regioes onde ha
% franjas para media 0 e desvio padrao 1
im = im - mean(im(maskind));
normim = im/std(im(maskind));

return
*/