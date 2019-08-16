/*
 * falange.c
 *
 *  Created on: 19/09/2011
 *      Author: Caius
 */

#include "..\..\Header Files\Falange\falange.hpp"
#include "..\..\Header Files\Falange\clahe.h"


using namespace cv;

void showProjH(IplImage *src,int *proj)
{
    cvZero(src);
    int j;
    //uchar *uproj=(uchar*)src->imageData;
    for(j=0;j<src->height;j++)
    {
        cvLine(src,cvPoint(src->width-1,j),cvPoint(src->width-1-proj[j],j),cvScalarAll(255),1,8,0);
    }
}


int falange(const char *str1)
{

    IplImage *img,*dst,*dst2,*dst3;
    int adaptX=3,adaptY=1,threshAdapt=80,histAdapt=180;

    int i,j,sum,maxSum,index,mean;
    CvMemStorage *storage = cvCreateMemStorage(0);
    char key;
    int flags = 8 + (255 << 8) + CV_FLOODFILL_FIXED_RANGE;

    CvConnectedComp comp;

    printf("Abrindo %s\n",str1);
    //obriga a imagem a ser nivel de cinza
    img=cvLoadImage(str1,0);
    if(!img)
    {
        printf("faltou imagem de entrada\n");
        return 0;
    }
    dst=cvCreateImage(cvGetSize(img),img->depth,1);
    dst2=cvCreateImage(cvGetSize(img),img->depth,1);
    dst3=cvCreateImage(cvGetSize(img),img->depth,1);

    int vetor[]={2,4,8,16};

    uchar* pdst2=(uchar*)dst2->imageData;
	//int projH[dst2->height];
	int *projH;
	projH = ( int * ) malloc ( dst2->height * sizeof(int) );

    /*
        cvCreateTrackbar("Block X","Adapt",&adaptX,3,callback);
        cvCreateTrackbar("Block Y","Adapt",&adaptY,3,callback);
        cvCreateTrackbar("bins hist","Adapt",&histAdapt,256,callback);
        cvCreateTrackbar("thresh","thresh_adapt",&threshAdapt,255,callback);
         */

    cvAdaptEqualize(img,dst,vetor[adaptX],vetor[adaptY],histAdapt,CV_CLAHE_RANGE_FULL);

    cvThreshold(dst,dst2,threshAdapt,255,CV_THRESH_BINARY_INV);

    //funcaoX(dst2,dst3,0.7,100);
    //cvShowImage("teste1",dst3);
    cvMorphologyEx(dst2,dst2,NULL,NULL,CV_MOP_OPEN,2);
    cvMorphologyEx(dst2,dst2,NULL,NULL,CV_MOP_CLOSE,2);





    //cvDilate(dst2,dst2,NULL,2);
    //cvErode(dst2,dst2,NULL,2);


    for(i=1;i<10;i++)
    {
        cvFloodFill(dst2,cvPoint((i*dst2->width)/11,dst2->height*0.95),cvScalarAll(0),cvScalarAll(0),cvScalarAll(200),&comp,flags,NULL);
        cvFloodFill(dst2,cvPoint((i*dst2->width)/11,dst2->height*0.05),cvScalarAll(0),cvScalarAll(0),cvScalarAll(200),&comp,flags,NULL);
    }

    //funcaoX(dst2,dst3,0.8,100);
    //cvShowImage("teste2",dst3);

    IplImage *IprojH = cvCreateImage(cvGetSize(dst2),dst2->depth,1);
    uchar *uProjH = (uchar*)IprojH->imageData;

    for(j=0;j<dst2->height;j++)
    {
        projH[j]=0;
        for(i=0;i<dst2->width;i++)
        {
            if(pdst2[i+j*dst2->widthStep]==255)
            {
                projH[j]++;
                uProjH[IprojH->width-projH[j]-1+j*IprojH->widthStep]=255;

            }
        }
    }
    //cvShowImage("ProjH",IprojH);

    for(index=0,maxSum=0,sum=0,i=0;i<dst2->height-3;i+=3)
    {
        sum=projH[i]+projH[i+1]+projH[i+2];
        if(sum>maxSum&&(i>0.1*dst2->height&&i<0.9*dst2->height))
        {
            maxSum=sum;
            index=i+1;

        }
    }

    for(mean=0,i=dst2->height*0.1;i<dst2->height*0.9;i++)
    {
        mean+=projH[i];

    }

    maxSum=maxSum/3;
    mean=mean/(0.8*dst2->height);

    //printf("(%d,%d) mean = %d\t\tsum = %d\n",dst2->width,dst2->height,mean,maxSum);

    for(i=dst2->height*0.1;i<dst2->height*0.9;i++)
    {
        if(projH[i]<mean)
        {
            projH[i]=0;
        }
    }

    if(maxSum>4*mean&&(maxSum>dst3->width*0.1||maxSum-mean>0.05*dst2->width))
    {
        cvLine(img,cvPoint(0,index),cvPoint(dst3->width-1,index),cvScalarAll(255),10,8,0);

    }

    showProjH(IprojH,projH);

    //cvShowImage("ProjH2",IprojH);
    cvShowImage("Result",dst2);

    cvShowImage("Source",img);

    key=cvWaitKey(0);
    cvReleaseImage(&dst);
    cvReleaseImage(&dst2);
    cvReleaseImage(&dst3);
    cvReleaseImage(&img);
    cvReleaseMemStorage(&storage);

    return index;
}

int phalange_QT(IplImage *original)
{
    int i,j,sum,maxSum,index=0,mean;
    CvMemStorage *storage = cvCreateMemStorage(0);

    int flags = 8 + (255 << 8) + CV_FLOODFILL_FIXED_RANGE;

    IplImage *img,*dst,*dst2,*dst3;
    int adaptX=3,adaptY=1,threshAdapt=80,histAdapt=180;

    CvConnectedComp comp;

    img=cvCreateImage(cvGetSize(original),original->depth,1);
    cvCvtColor(original,img,CV_BGR2GRAY);

    dst=cvCreateImage(cvGetSize(img),img->depth,1);
    dst2=cvCreateImage(cvGetSize(img),img->depth,1);
    dst3=cvCreateImage(cvGetSize(img),img->depth,1);

    int vetor[]={2,4,8,16};

    uchar* pdst2=(uchar*)dst2->imageData;
    //int projH[dst2->height];
	int *projH;
	projH = (int *)malloc(dst2->height * sizeof(int));

    //cvAdaptEqualize(img,dst,vetor[adaptX],vetor[adaptY],histAdapt,CV_CLAHE_RANGE_FULL);

    cvThreshold(dst,dst2,threshAdapt,255,CV_THRESH_BINARY_INV);

    cvMorphologyEx(dst2,dst2,NULL,NULL,CV_MOP_OPEN,2);
    cvMorphologyEx(dst2,dst2,NULL,NULL,CV_MOP_CLOSE,2);

    for(i=1;i<10;i++)
    {
        cvFloodFill(dst2,cvPoint((i*dst2->width)/11,dst2->height*0.95),cvScalarAll(0),cvScalarAll(0),cvScalarAll(200),&comp,flags,NULL);
        cvFloodFill(dst2,cvPoint((i*dst2->width)/11,dst2->height*0.05),cvScalarAll(0),cvScalarAll(0),cvScalarAll(200),&comp,flags,NULL);
    }


    IplImage *IprojH = cvCreateImage(cvGetSize(dst2),dst2->depth,1);
    uchar *uProjH = (uchar*)IprojH->imageData;

    for(j=0;j<dst2->height;j++)
    {
        projH[j]=0;
        for(i=0;i<dst2->width;i++)
        {
            if(pdst2[i+j*dst2->widthStep]==255)
            {
                projH[j]++;
                uProjH[IprojH->width-projH[j]-1+j*IprojH->widthStep]=255;

            }
        }
    }

    // Main horizontal project count
    for(index=0,maxSum=0,sum=0,i=0;i<dst2->height-4;i+=3)
    {
        sum=projH[i]+projH[i+1]+projH[i+2];
        if(sum>maxSum&&(i>0.1*dst2->height&&i<0.9*dst2->height))
        {
            maxSum=sum;
            index=i+1;

        }
    }

    for(mean=0,i=dst2->height*0.1;i<dst2->height*0.9;i++)
    {
        mean+=projH[i];

    }

    maxSum=maxSum/3;
    mean=mean/(0.8*dst2->height);


    for(i=dst2->height*0.1;i<dst2->height*0.9;i++)
    {
        if(projH[i]<mean)
        {
            projH[i]=0;
        }
    }

    if(maxSum>4*mean&&(maxSum>dst3->width*0.1||maxSum-mean>0.05*dst2->width))
    {
        //cvLine(img,cvPoint(0,index),cvPoint(dst3->width-1,index),cvScalarAll(255),10,8,0);
    }
    else
    {
        index=0;
    }


    cvReleaseImage(&dst);
    cvReleaseImage(&dst2);
    cvReleaseImage(&dst3);
    cvReleaseImage(&img);
    cvReleaseMemStorage(&storage);

    return index;
}
