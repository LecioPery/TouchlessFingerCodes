#include "..\Header Files\num2strSpecial.hpp"

#include <iostream>
#include <cstdlib>

#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

String num2strSpecial(int n)
{
	String s = "";
	int d;
	int ncopy = n;

	while (n > 0)
	{
		d = n % 10;
		s = to_string(d) + s;
		n = (n - d) / 10;
	}
	if (ncopy < 10) s = "0" + s; /* "Bonus". */
	return s;
}

/**
Matlab equivalent code:
%Note: this might be eventualy deleted!
%Main reason is that, while num2str does that, this function was created due to "coder" plugin from Matlab.

function s = num2strSpecial(n)
s = [];
while n > 0
d = mod(n, 10);
s = [char(48 + d), s];
n = (n - d) / 10;
end
end
*/