#include <iostream>

#include "..\Header Files\inserezeros.hpp"
#include "..\Header Files\stddvim.hpp"
#include "..\Header Files\trickylib.hpp"

using namespace cv;

Mat stddvim(cv::Mat imagem, int blksze)
{

	int C;
	int L;
	int numBlocosL;
	int numBlocosC;
	int i;
	int j;

	float desviopadrao;
	
	Mat imagemstdv;
	Mat ImagemN;
	Mat helper;

	/* [L, C] = size(imagem); */
	C = imagem.cols;
	L = imagem.rows;

	/* [ImagemN, numBlocosL, numBlocosC]  = inserezeros(imagem, blksze); */
	inserezeros inserezerosObj(imagem, blksze);
	ImagemN = inserezerosObj.getzeroIM();
	numBlocosL = inserezerosObj.getnumBlocosL();
	numBlocosC = inserezerosObj.getnumBlocosC();

	for (i = 0; i < numBlocosL; i++)
	{

		for (j = 0; j < numBlocosC; j++)
		{

			/**
			Ocorre que, implementar linha por linha do código em Matlab é impossível.
			Ao invés disso, interpretei o código e tentei reproduzir as mesmas variáveis de saída.
			Da documentação do matlab: std2 computes the standard deviation of the array A using std(A(:)).
			Ou seja, desvio padrão considerando todos os valores da imagem recebida como argumento.
			No caso, a imagem é o subbloco que vai de 0 até LimL e de 0 até LimC:
			**/
			desviopadrao = std2(ImagemN(Rect(j * blksze, i * blksze, blksze, blksze)));
			ImagemN(Rect(j * blksze, i * blksze, blksze, blksze)) = desviopadrao;

		}

	}
	imagemstdv = ImagemN( Rect(0, 0, C, L) );

	/* Se conveniente...
	imwrite("Desviopadrao2.bmp", imagemstdv);*/

	return imagemstdv;
}

/**
Matlab equivalent code:

function imagemstdv = stdvim(imagem, blksze);
coder.extrinsic('std2');
% Determina dimensoes da imagem
[L, C] = size(imagem);

% Ajusta a dimensao da imagem, inserindo zeros, de maneira a gerar um numero
% inteiro de sub-blocos de dimensao blksze
[ImagemN, numBlocosL, numBlocosC]  = inserezeros(imagem, blksze);

% Processa a imagem sub-bloco a sub-bloco, de blksze pixels, trocando os
% elementos do sub-bloco pelo desvio padrao do subbloco
for i=0:(numBlocosL-1)
for j=0:(numBlocosC-1)
LimL = (i*blksze+1):((i+1)*blksze);
LimC = (j*blksze+1):((j+1)*blksze);
desviopadrao = std2(ImagemN(LimL, LimC));
ImagemN(LimL, LimC) = desviopadrao;
end
end

% Retorna uma imagem com a mesma dimensao de imagem
imagemstdv = ImagemN(1:L, 1:C);


return
*/