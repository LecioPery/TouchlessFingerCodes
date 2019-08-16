#ifndef __FINDCENTERROI
#define __FINDCENTERROI

#include <opencv2/highgui/highgui.hpp>

using namespace cv;

class findCenterRoi
{

private:
	std::vector<vector<int>> CandRefiFragments;
	std::vector<vector<int>> CandRefjFragments;
	vector<int> x;
	vector<int> y;
	int indexOfHighest;

	vector<vector<int>> fragmentVector(vector<int> input);

public:
	void findCenterRoi::getResult(std::vector<int> *xcenter, std::vector<int> *ycenter);
	findCenterRoi(Mat imageTarget, std::string filename, bool useCrop);
	std::vector<std::vector<int>> getCandRefi();
	std::vector<std::vector<int>> getCandRefj();
	int getIndexOfHighest();

};

#endif