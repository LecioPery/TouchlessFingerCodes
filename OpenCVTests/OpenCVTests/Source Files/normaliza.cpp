#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include "..\Header Files\normaliza.hpp"
#include "..\Header Files\trickyLib.hpp"

using namespace cv;
using namespace std;

/* Completo. */
normaliza::normaliza(Mat im)
{
	int channels = im.channels();
	int nRows = im.rows;
	int nCols = im.cols;
	int i, j, k;
	float average;
	float deviation = 0;
	float *p;

	if (getImageType(im.type()) != "CV_32FC3")
	{
		im.convertTo(im, CV_32FC3, 1.0, 0.0);
	}
	
	/* std(im(:)); é o desvio padrão de todos os valores da imagem. */
	deviation = standard_deviation(im.ptr<float>(0), nRows * nCols * channels);
	average = mean(im.ptr<float>(0), nRows * nCols * channels);

	for (k = 0; k < channels; k++)
	{
		/** "channel: k"*/
		for (i = 0; i < nRows; ++i)
		{
			p = im.ptr<float>(i);
			for (j = k * nCols; j < (k + 1) * nCols; ++j)
			{
				/* n = im - mean(im(:)); */
				p[j] = p[j] - average;
				/* n = n/std(im(:)); */
				p[j] = p[j] / deviation;
				
			}
		}
	}
	im.copyTo(this->im);
}

Mat normaliza::getIm()
{
	return im;
}

/**
Matlab equivalent code:
% NORMALIZA - Normaliza os valores da imagem para media 0 e desvio padrao 1

function n = normaliza(im)

% Verifica se a imagem foi convertida de uint8 para double.
% Em caso negativo realiza a conversao.
if ~isa(im,'double')
im = double(im);
end

% Normaliza para media zero e desvio padrao unitario
n = im - mean(im(:));
n = n/std(im(:));

return

*/