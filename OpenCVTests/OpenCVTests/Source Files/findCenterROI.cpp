#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <cstdlib>

#include "..\Header Files\num2strSpecial.hpp"
#include "..\Header Files\preparaim.hpp"
#include "..\Header Files\extraidirecional.hpp"
#include "..\Header Files\trickyLib.hpp"
#include "..\Header Files\commonConstants.hpp"

#include "..\Header Files\findCenterRoi.hpp"

/* MATLAB constants: */
#define VERBOSE 1
//#define DEBUG_CODE
//#define DEBUG_CODE2

/* MATLAB Boundaries: */
#define ORIENT_VEC_WIDTH 99

using namespace cv;
using namespace std;

/* Note: precision issues in regard of Matlab: */
//int main(int argc, char **argv)
findCenterRoi::findCenterRoi( Mat imageTarget, std::string filename, bool useCrop )
{
	
	/* Matlab variables */
	Mat orientVec;
	
	/* Other Matlab variables */
	Mat img;
	Mat imgCrop;
	Mat orientimCrop;
	Mat maskCrop;
	Mat fingerprint;
	Mat orientim;
	Mat mask;
	Mat NHOOD;
	Mat E;
	Mat M;
	Mat Mfilter;
	Mat Mbin;
	Mat outputSaver;
	String test;
	int imgheight;
	int imgwidth;
	int lmask;
	int cmask;
	int Bordax;
	int Borday;
	int Limiar;
	/*int Liml;
	int Limi;
	int Lims;*/
	const float PercentBorda = (float) 0.3;//(float) 0.07;
	std::vector<int> CandRefi;
	std::vector<int> CandRefj;
	std::vector<int> iCand;
	std::vector<int> Refi;
	std::vector<int> Refj;
	
	String str = IMG_OUTPUT_PATH;

	/* My own variables (helpers): */
	/* Old but gold:
	int lastY;
	*/
	int currentX;
	int currentY;
	/* Old but gold:
	int counter;
	*/
	#ifdef DEBUG_CODE
	Mat imageRoi;
	Rect regionOfInterest;
	#endif
	
	img = imageTarget;
	fingerprint = preparaim(img); /* Pronto */
	extraidirecional *direcional = new extraidirecional(fingerprint, 16);
	mask = direcional->getMask();
	orientim = direcional->getOrientim();
	imgheight = orientim.size().height;
	imgwidth = orientim.size().width;
	Bordax = (int)floor(PercentBorda * imgwidth);
	Borday = (int)floor(PercentBorda * imgheight);
	
	if(useCrop) Borday += (int) floor(0.1 * imgheight); /* Cut more 10% */

	/* Subtract 1 in both coordinates due to Matlab's indexing system that differs from C++'s. */
	(img(Range(Borday - 1, imgheight - Borday), Range(Bordax - 1, imgwidth - Bordax))).copyTo(imgCrop);
	#ifdef DEBUG_CODE
	namedWindow("imgCrop", CV_WINDOW_AUTOSIZE);
	imshow("imgCrop", imgCrop);
	#endif
	
	/* Subtract 1 in both coordinates due to Matlab's indexing system that differs from C++'s. */
	(orientim(Range(Borday - 1, imgheight - Borday), Range(Bordax - 1, imgwidth - Bordax))).copyTo(orientimCrop);
	if (OUTPUT_INTERMEDIATE) { outputSaver = normalizaImg(orientim); cv::imwrite(str + "intermediate/orientim" + filename + ".jpg", outputSaver); }
	if (OUTPUT_INTERMEDIATE) { outputSaver = normalizaImg(orientimCrop); cv::imwrite(str + "intermediate/orientimCrop" + filename + ".jpg", outputSaver); }
	
	NHOOD = Mat::ones(55, 55, CV_32F);
	lmask = mask.size().height;
	cmask = mask.size().width;
	
	// imshow(mask);
	mask = setMatValuesFloat(mask, Range(1, 15), Range(1, mask.size().width), 0);
	mask = setMatValuesFloat(mask, Range(1, mask.size().height), Range(1, 15), 0);
	mask = setMatValuesFloat(mask, Range(lmask - 15, lmask), Range(1, mask.size().width), 0);
	mask = setMatValuesFloat(mask, Range(1, mask.size().height), Range(cmask - 15, cmask), 0);
	if (OUTPUT_INTERMEDIATE) cv::imwrite(str + "intermediate/mask" + filename + ".jpg", mask);
	
	/*#ifdef DEBUG_CODE
	namedWindow("mask", CV_WINDOW_AUTOSIZE);
	imshow("mask", direcional.getMask());
	#endif*/
	
	erode(mask, mask, NHOOD);
	
	/* Subtract 1 in both coordinates due to Matlab's indexing system that differs from C++'s. */
	(mask(Range(Borday - 1, imgheight - Borday), Range(Bordax - 1, imgwidth - Bordax))).copyTo(maskCrop);
	
	#ifdef DEBUG_CODE
	namedWindow("maskCrop", CV_WINDOW_AUTOSIZE);
	imshow("maskCrop", maskCrop);
	#endif
	
	E = matFunction(orientimCrop, "sin");
	if (OUTPUT_INTERMEDIATE) { outputSaver = normalizaImg(E); cv::imwrite(str + "intermediate/E = sin(orientim)" + filename + ".jpg", outputSaver); }
	
	#ifdef DEBUG_CODE2
	namedWindow("E", CV_WINDOW_AUTOSIZE);
	imshow("E", E);
	#endif
	#ifdef DEBUG_CODE
	namedWindow("E * mask", CV_WINDOW_AUTOSIZE);
	imshow("E * mask", matProduct(E, bytewiseUpdateImage(maskCrop, 255, METHOD_DIVIDE), METHOD_PRODUCT));
	#endif
	
	M = gradientMat(E);
	M = tailMat(M);
	if (OUTPUT_INTERMEDIATE) { outputSaver = normalizaImg(M); cv::imwrite(str + "intermediate/M = Sobel(E)" + filename + ".jpg", outputSaver); }
	
	#ifdef DEBUG_CODE
	namedWindow("M * mask", CV_WINDOW_AUTOSIZE);
	imshow("M * mask", matProduct(M, bytewiseUpdateImage(maskCrop, 255, METHOD_DIVIDE), METHOD_PRODUCT));
	#endif
	
	M.copyTo(Mfilter);
	Mfilter = normalizaImg(Mfilter);
	if (OUTPUT_INTERMEDIATE) { outputSaver = normalizaImg(Mfilter); cv::imwrite(str + "intermediate/Mfilter" + filename + ".jpg", outputSaver); }
	
	#ifdef DEBUG_CODE
	Mat MfilterAux = matProduct(Mfilter, bytewiseUpdateImage(maskCrop, 255, METHOD_DIVIDE), METHOD_PRODUCT);
	MfilterAux.convertTo(MfilterAux, CV_8U);
	namedWindow("uint8(maskCrop.*Mfilter)", CV_WINDOW_AUTOSIZE);
	imshow("uint8(maskCrop.*Mfilter)", MfilterAux);
	#endif
	
	Limiar = 140;
	
	Mbin = matProduct(Mfilter, bytewiseUpdateImage(maskCrop, 255, METHOD_DIVIDE), METHOD_PRODUCT);
	cv::threshold(Mbin, Mbin, Limiar, 1, THRESH_BINARY);
	if (OUTPUT_INTERMEDIATE) { outputSaver = normalizaImg(Mbin); cv::imwrite(str + "intermediate/Mbin" + filename + ".jpg", outputSaver); }
	
	#ifdef DEBUG_CODE
	namedWindow("Mbin", CV_WINDOW_AUTOSIZE);
	imshow("Mbin", Mbin);
	#endif

	findValueMat(Mbin, CandRefj, CandRefi);

	/* Old but gold:
	CandRefiFragments = fragmentVector(CandRefi);
	CandRefjFragments = fragmentVector(CandRefj);*/
	
	if (CandRefi.empty())
	{
		orientVec = Mat::zeros(1, ORIENT_VEC_WIDTH, CV_32F);
		indexOfHighest = 0;
		x.push_back(0);
		y.push_back(0);
		std::cout << "No center has been found" << std::endl;
	}
	else
	{
		
		indexOfHighest = 0;
		/* Old but gold:
		counter = 0;
		lastY = 2000; /* Any big number
		while( counter < CandRefiFragments.size() )
		{*/
			
			//A = A + 1;
			/* jCand = [1, 1, ..., 1] => don't care */
			iCand = findValueVectorInt(CandRefi, *std::min_element(std::begin(CandRefi), std::end(CandRefi)));
			Refi = subvectorInt(CandRefi, iCand);
			Refj = subvectorInt(CandRefj, iCand);

			/* Old but gold:
			iCand = findValueVectorInt(CandRefiFragments.at(counter), *std::min_element(std::begin(CandRefiFragments.at(counter)), std::end(CandRefiFragments.at(counter))));
			Refi = subvectorInt(CandRefiFragments.at(counter), iCand);
			Refj = subvectorInt(CandRefjFragments.at(counter), iCand);
			*/
			
			/**
			* Refi - 5:Refi + 5, Refj - 5:Refj + 5:
			* Refi - 5 => Rect(x, ..., 10, ...);
			* Refj - 5 => Rect(..., y, ..., 10);
			* (Matlab) -5 = -6 (OpenCV)
			**/
			currentY = (*std::min_element(std::begin(Refi), std::end(Refi))) - 6;
			currentX = (*std::min_element(std::begin(Refj), std::end(Refj))) - 6;

			//Passo 1: Calcular o centro da imagem (mediana);
			/*int centerX;
			int centerY;
			int i;
			int lastDistance;
			int currentDistance;
			centerX = Mbin.size().width / 2;
			centerY = Mbin.size().height / 2;
			currentX = 0;
			currentY = 0;
			lastDistance = 0x7FFFFFFF; /* Infinite: */
			/*i = 0;
			//Passo 2: para todo elemento candidato:
			do
			{
				//Passo 3: pelo quadrado da distância euclidiana, tomar a que tem o menor:
				currentDistance = (CandRefi[i] - centerX) * (CandRefi[i] - centerX) + (CandRefj[i] - centerY) * (CandRefj[i] - centerY);
				if ( currentDistance < lastDistance )
				{
					currentX = CandRefi[i];
					currentY = CandRefj[i];
					lastDistance = currentDistance;
				}
				i++;
			} while ( i < CandRefi.size() );*/
			
			#ifdef DEBUG_CODE
			regionOfInterest = Rect(y, x, 11, 11);
			imageRoi = imgCrop(regionOfInterest);
			imageRoi = 0; /* Matlab style */
			#endif
			
			/* Correct to un-cropped image: */
			currentY += Borday - 1;
			currentX += Bordax - 1;
			
			#ifdef DEBUG_CODE
			namedWindow("imgCrop", CV_WINDOW_AUTOSIZE);
			imshow("imgCrop", imgCrop);
			#endif
			
			/* The following is incomplete: */
			/*
			Liml = 50;
			Limi = 20;
			Lims = 100;
			
			/*regionOfInterest = Rect(x, y, 11, 11);
			RegInteresse = orientimCrop(Refi - Limi:Refi + Lims, Refj - Liml : Refj + Liml);
			
			/* Extrai a região de interesse da imagem orifinal */
			/*RegInteresseImg = imgCrop(Refi - Limi:Refi + Lims, Refj - Liml : Refj + Liml);*/
			x.push_back(currentX);
			y.push_back(currentY);
			/* Old but gold:
			if (currentY < lastY) indexOfHighest = counter;
			lastY = currentY;
			counter++;
			*/
			
		/* Old but gold:
		}*/
		
	}
	#ifdef DEBUG_CODE
	waitKey(0);
	#endif // DEBUG_CODE
	
	delete direcional;
	
}

vector<vector<int>> findCenterRoi::fragmentVector(vector<int> input)
{
	
	vector<vector<int>> result;
	int aux1;
	int aux2;
	int counter;
	
	counter = 0;
	aux1 = 0;
	aux2 = 0;
	while (counter < input.size())
	{
		
		vector<int> newVector;
		/* Enquanto a diferença não for muito grande: */
		do
		{
			
			aux1 = input[counter];
			newVector.push_back(aux1);
			counter++;
			if (counter == input.size()) break;
			aux2 = input[counter];
			
		} while ((aux1 - aux2 > -10) && (aux1 - aux2 < 10));
		result.push_back(newVector);
		
	}
	return result;
	
}

std::vector<vector<int>> findCenterRoi::getCandRefi()
{
	return CandRefiFragments;
}
std::vector<vector<int>> findCenterRoi::getCandRefj()
{
	return CandRefjFragments;
}

int findCenterRoi::getIndexOfHighest()
{
	return indexOfHighest;
}

/*void findCandidate()
{
	/* jCand = [1, 1, ..., 1] => don't care *//*
	iCand = findValueVectorInt(CandRefi, *std::min_element(std::begin(CandRefi), std::end(CandRefi)));
	Refi = subvectorInt(CandRefi, iCand);
	Refj = subvectorInt(CandRefj, iCand);
}*/

void findCenterRoi::getResult(vector<int> *xcenter, vector<int> *ycenter)
{
	*xcenter = findCenterRoi::x;
	*ycenter = findCenterRoi::y;
	return;
}

void save()
{
	/* Para salvar o resultado final:
	imwrite("Desviopadrao2.bmp", imagemstdv);*/
	return;
}

/*
Matlab equivalent code:
% Imagens direcionais
%-------------------------
% Programa que, dado um grupo de fingerprints, extrai suas imagens
% direcionais

clc;
clear all
close all

% Mostra resultados
verbose = 1;

% Debug
debug_code = 1;

A = 0;

% Numero de imagens
Nimg = 2; % 25

% Extrai as imagens direcionais
% -----------------------------
orientMAT = [];
conta = 1;

% Pré-alocação:
imagem = zeros([240 336], 'uint8');
map = zeros([256 3], 'uint8');
M = 3.7;

% Abre a imagem
for i=1:Nimg %i=1:Nimg

if i<10
ind = ['0' num2strSpecial(i)];
else
ind = num2strSpecial(i);
end

[imagem,map]=imread(['images/' ind '.bmp']);

img = imagem;

% Prepara Imagem
fingerprint = preparaim(img);


% Utiliza o algoritmo de extração de minucias para gerar a imagem
% direcional orietim
[orientim, mask] = extraidirecional(fingerprint,16);

% Determina a largura da borda
[imgheight, imgwidth] = size(orientim);

% Determina limietes para realizar Crop na imagem
PercentBorda = 0.07;
Bordax = floor(PercentBorda*imgwidth);
Borday = floor(PercentBorda*imgheight);

% Faz crop na imagem original
imgCrop = img(Borday:(imgheight-Borday), Bordax:(imgwidth-Bordax));

if debug_code
figure(1);
imshow(imgCrop);
pause(1.0)
end

% Determina as direções das franjas com base na imagem binarizada já
% com o Crop
orientimCrop = orientim(Borday:(imgheight-Borday), Bordax:(imgwidth-Bordax));
%         hold on
%         s_orient = plotridgeorient(orientimCrop, 10);
%         pause(1.0)
%         hold off
%
% Erode a máscara

NHOOD = ones(55,55);
[lmask, cmask] = size(mask);

% imshow(mask);
mask(1:15,:) = 0;
mask(:,1:15) = 0;
mask(lmask-15:lmask,:) = 0;
mask(:,cmask-15:cmask) = 0;

%         figure
%         imshow(mask)

mask = imerode(mask,NHOOD);

% Faz o Crop da máscara
maskCrop = mask(Borday:(imgheight-Borday), Bordax:(imgwidth-Bordax));

if debug_code
figure(2);
imshow(maskCrop);
pause(1.0)
end


% Aplica o seno na imagem de direções
E = sin(orientimCrop);

if debug_code
figure(3);
imshow(maskCrop.*E);
pause(1.0)
end

% Gradiente em A
[l, c] = size(E);
for k=2:l-1
for j=2:c-1
Gx = (E(k+1, j-1)+2*E(k+1, j)+E(k+1,j+1))-(E(k-1, j-1)+2*E(k-1, j)+E(k-1,j+1));
Gy = (E(k+1, j+1)+2*E(k, j+1)+E(k-1,j+1))-(E(k+1, j-1)+2*E(k, j-1)+E(k-1,j-1));
M(k,j) = sqrt(Gx^2 + Gy^2);
end
end

M(l,:)=0;
M(:,c)=0;

if debug_code
figure(4)
imshow(maskCrop.*M)
end

%         M = round(255*M/max(max(M)));
%         pause(1.0)
%
%         % Gradiente em A
%         [l, c] = size(E);
%         for k=2:l-1
%             for j=2:c-1
%                 Gx = (M(k+1, j-1)+2*M(k+1, j)+M(k+1,j+1))-(M(k-1, j-1)+2*M(k-1, j)+M(k-1,j+1));
%                 Gy = (M(k+1, j+1)+2*M(k, j+1)+M(k-1,j+1))-(M(k+1, j-1)+2*M(k, j-1)+M(k-1,j-1));
%                 Mfilter(k,j) = sqrt(Gx^2 + Gy^2);
%             end
%         end
%
%         Mfilter(l,:)=0;
%         Mfilter(:,c)=0;
%
Mfilter = M;

Mfilter = round(255*Mfilter/max(max(Mfilter)));

% H = fspecial('sobel');
% Mfilter = round(imfilter(M,H));

if debug_code
figure(5)
imshow(uint8(maskCrop.*Mfilter));
end

Limiar = 140;

% Limiarização para detectar a referência
Mbin = Mfilter.*maskCrop > Limiar;

if debug_code
figure(6);
imshow(Mbin);
end

% Integra as direções
% A = integrapix(E);

% Determina as coordenadas da referência
[CandRefi,CandRefj] = find(Mbin == 1);

% Se não existem cadidatos a core, o fingerprint é do tipo arch

if isempty(CandRefi)
orientVec = zeros(1,99);
else

% disp('Opaaaa!');
A = A+1;
% pause

% Se existem candidatos define qual candidato será escolhido
[iCand, jCand] = find(CandRefi == min(CandRefi));

Refi = CandRefi(iCand);
Refj = CandRefj(iCand);

% Mostra a imagem com o core
imgCrop(Refi-5:Refi+5, Refj-5:Refj+5) = 0;

if debug_code
figure(7);
imshow(imgCrop)
pause
end

% Extrai a região de interesse na imagem direcional
Liml = 50;
Limi = 20;
Lims = 100;
RegInteresse = orientimCrop(Refi-Limi:Refi+Lims, Refj-Liml:Refj+Liml);

% Extrai a região de interesse da imagem orifinal
RegInteresseImg = imgCrop(Refi-Limi:Refi+Lims, Refj-Liml:Refj+Liml);

%             figure(7);
%             imshow(RegInteresseImg);

% Orientaçao sobre a Região de Interesse
% RegInteresse_SS é a matriz RegInteresse subamostrada em intervalos
% de 8 pixels na horizontal e na vertical

RegInteresse_SS = ssamporientimg(RegInteresse, 10);

%             hold on
%             s_orient = plotridgeorient(RegInteresse, 10);
%             pause(1.0)
%             hold off

[h,w] = size(RegInteresse_SS);
orientVec = reshape(RegInteresse_SS,1,h*w);

end % Fim isempty(CandRefi)

% Define a matriz onde as direções para cada imagem serão armazenadas
orientMAT(conta,:) = orientVec;

conta = conta+1;

end % Fim Nimg

close all;

% save orient.mat orientMAT;
*/