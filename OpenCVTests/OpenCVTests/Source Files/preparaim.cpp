#include <cmath>
#include <opencv2/highgui/highgui.hpp>
#include "..\Header Files\preparaim.hpp"
#include "..\Header Files\trickyLib.hpp"

#include <iostream>
#include <cstdlib>

using namespace cv;
using namespace std;

/** Prepares image for mathematical floating-point operations: */
cv::Mat preparaim( cv::Mat img )
{
	Mat immagem;

	/* immagem=double(img): */
	img.convertTo(immagem, CV_32FC3, 1.0, 0.0);
	/* Went from 8UC3 to 32FC3 */

	Mat fingerprint = immagem;

	int graylevmax = 0; /*Unused*/

	if (getImageType(img.type()) == "CV_8UC3")
	{
		graylevmax = ( int ) pow(2, 8) - 1;
	}
	else if (getImageType(img.type()) == "CV_16UC3")
	{
		graylevmax = ( int ) pow(2, 16) - 1;
	}
	else if (getImageType(img.type()) == "CV_32UC3")
	{
		graylevmax = ( int ) pow(2, 32) - 1;
	}

	return fingerprint;

	/* Eu realmente não entendi o que essa função supostamente faz.
	* Tudo o que compreendi é que ela transforma cada valor em um double e então retorna perdendo o valor. */
}

/**
Matlab equivalent code:
function fingerprint = preparaim(img)

immagem=double(img);

%ajustar a intensidade da imagem
if isa(img,'uint8')
graylevmax=2^8-1;
end

if isa(img,'uint16')
graylevmax=2^16-1;
end

if isa(img,'uint32')
graylevmax=2^32-1;
end

fingerprint = immagem;

return
*/