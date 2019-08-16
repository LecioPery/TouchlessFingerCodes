#include "..\Header Files\extraidirecional.hpp"
#include "..\Header Files\segmentafranja.hpp"
#include "..\Header Files\orientfranja.hpp"
#include "..\Header Files\trickyLib.hpp"

#include <iostream>

using namespace cv;
using namespace std;

extraidirecional::extraidirecional(Mat fingerprint, int blksze)
{

	segmentafranja *franjaSegmentada = new segmentafranja(fingerprint, blksze); /* Feito */
	
	this->mask = franjaSegmentada->getMask();
	this->orientim = orientfranja( franjaSegmentada->getNormim() ).getOrientim();

	delete franjaSegmentada;

}

Mat extraidirecional::getMask(void)
{
	return mask;
}

Mat extraidirecional::getOrientim(void)
{
	return orientim;
}

/**
Matlab equivalent code:
function [orientim, mask] = extraidirecional(fingerprint, blksze);

% Normaliza a imagem
[normim, mask] = segmentafranja(fingerprint, blksze);

% Determinar a orientacao
orientim = orientfranja(normim);

return
*/