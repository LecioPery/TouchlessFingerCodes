#include "..\Header Files\generateCenters.hpp"
#include "..\Header Files\findCenterRoi.hpp"
#include "..\Header Files\trickyLib.hpp"
#include "..\Header Files\commonConstants.hpp"

#define PREFIX "SIRE-1_"
#define SUFIX "_0"
#define EXTENSION ".jpg"

/*
RIP MIDAS:
#include "..\Header Files\Deprecated\h12functions.hpp"
Use:
Mat HT12;
HT12 = create_h12((char *)(address + "SIRE-1_" + std::to_string(i) + "_" + std::to_string(j) + "_HT1.jpg").c_str(), (char *)(address + "SIRE-1_" + std::to_string(i) + "_" + std::to_string(j) + "_HT2.bmp").c_str());
*/

/*My own constants: */
#define USE_CROP false

void centerGenerate( char *imgPath, char *centersPath )
{
	
	Mat imageTarget;
	int k;
	int i;
	int j;
	
	vector<int> x;
	vector<int> y;
	
	vector<findCenterRoi> images;
	
	std::string address(imgPath);
	
	for (k = 1; k <= NPEOPLE; k++)
	{
		
		for (i = 1; i <= NFINGER; i++)
		{
			
			for (j = 1; j <= NSAMPLE; j++)
			{
				
				try
				{
					
					imageTarget = imread(address + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + SUFIX + EXTENSION, CV_8UC1);
					findCenterRoi *image = new findCenterRoi(imageTarget, std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j), USE_CROP);
					image->getResult(&x, &y);
					if (x.empty())
					{
						saveCoordinates(centersPath + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + ".txt", 0, 0);
					}
					else saveCoordinates(centersPath + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) + ".txt", x.at(0), y.at(0));
					delete image;
					
				}
				catch (Exception e)
				{
					
					std::cout << e.msg << std::endl;
					std::cout << i;
					system("PAUSE");
					exit(-1);
					
				}
				std::cout << address + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j) << std::endl;
				
			}
			
		}
		
	}
	return;
	
}
