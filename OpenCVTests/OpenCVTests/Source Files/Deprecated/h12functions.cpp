/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Bruno Macchiavello  bruno@cic.unb.br             %
%          Alexandre Zaghetto  alexandre@cic.unb.br         %
%          Mamede Lima-Marques limamarques@gmail.com        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 30/11/2010					                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UnB  - Universidade de Brasília                           %
% IE   - Institute of Exact Sciences                        %
% CIC  - Department of Computer Science                     %
% LISA - Laboratory of Image, Signal and Audio              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPAI - Centre for Research on Architecture of Information %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cv.h"
#include "highgui.h"
#include "..\..\Header Files\Deprecated\h12functions.hpp"


cv::Mat create_h12(char *str1,char *str2)
{

	cv::Mat result;

	/**
		The following is the Legacy Code of MIDAS.
	*/

    char *JpegFile,*strmindtct,*strcjpegb,*strtam,*BmpFile,*strdir, *strdel, *ptr;
    float  higthreshold, lowthreshold;
    int filterSze;
    
    IplImage *Imgh1; //  Imagem origial HT1
    IplImage *Imgh1Proc;    //  Processed HT1
    IplImage *Imgh2; //  Imagem origial HT2
    IplImage *Binh2; //  Binary version of HT2
    IplImage *NBinh2; //  Negative version of Binary image of HT2
    IplImage *Mask; //  Binary Mask
    IplImage *Mask2; //  Binary Mask
    IplImage *Edges; //  Edges Mask
    IplImage *NEdges; //  Inverse Edges Mask
    IplImage *FMask2; //  smooth version of Binary Mask
    IplImage *Imgh12; //  Final HT12 image
    
    
    JpegFile = (char*)malloc( sizeof(char)* (strlen(str2)+10));
    strmindtct = (char*)malloc( sizeof(char)* (2*strlen(str2)+20));
    strcjpegb = (char*)malloc( sizeof(char)* (strlen(str2)+50));
    strtam = (char*)malloc( sizeof(char)* 5);
    strdir = (char*)malloc( sizeof(char)* (strlen(str2)+80));
    strdel = (char*)malloc( sizeof(char)* (strlen(str2)+100));
    BmpFile = (char*)malloc( sizeof(char)* (strlen(str2)+10));
    
    //load image h1
     Imgh1 = cvLoadImage(str1,0);     
     
     // Verify if load was successful
    if(!Imgh1) { 
                 printf("Could not load input H1 image file: '%s'. \n",str1);   
                 system("PAUSE");
                 exit(-1);
    }
    
     //load image h2
     Imgh2 = cvLoadImage(str2,0);     
     
     // Verify if load was successful
    if(!Imgh2) { 
                 printf("Could not load input H2 image file: '%s'. \n",str2);   
                 system("PAUSE");
                 exit(-1);
    }
    
   
   //pre-processing of HT1 image (create binary image)
   //-------------------------------------------------
   
   //create a JPEG version
   strcpy(JpegFile,str2);
   strcpy(JpegFile+(strlen(str2)-4),".jpg");
   
   if(!cvSaveImage(JpegFile, Imgh2, NULL)) printf("Could not save: %s\n",JpegFile);
   
   //create command string for mindtct
   // system(['mindtct -b ' filename_B '.jpg ' filename_B '_output'] );
   strcpy(strmindtct,"mindtct -b ");
   strcat(strmindtct,JpegFile);
   strcat(strmindtct," ");
   strcat(strmindtct,JpegFile);
   strcpy(strmindtct+(strlen(strmindtct)-4),"_output");
   
   system(strmindtct);
   free(strmindtct);
   
    //create command string for cjpegb
   // system(['cjpegb 95 jpg ' filename_B '_output.brw -r ' num2str(w) ',' num2str(h) ',8']);
   strcpy(strcjpegb,"cjpegb 95 jpg ");
   strcat(strcjpegb,JpegFile);
   strcpy(strcjpegb+(strlen(strcjpegb)-4),"_output.brw -r ");
   _itoa(Imgh2->width,strtam,10);
   strcat(strcjpegb,strtam);
   strcat(strcjpegb,",");
   _itoa(Imgh2->height,strtam,10);  
   strcat(strcjpegb,strtam);
   strcat(strcjpegb,",8");
   
   system(strcjpegb);
   free(strcjpegb);
   free(strtam);
   
   //creation of HT12 image
   //----------------------
   
   
   //load binay image;
   strcpy(JpegFile+(strlen(str2)-4),"_output.jpg");
   Binh2 = cvLoadImage(JpegFile,0);     
     
   // Verify if load was successful
   if(!Binh2) { 
                 printf("Could not load input H1 image file: '%s'. \n",str1);   
                 system("PAUSE");
                 exit(-1);
   }
   
   free(JpegFile);
   
   cvThreshold(Binh2,Binh2, 127, 255, CV_THRESH_BINARY);
   //create inverse and mask
   NBinh2=cvCreateImage(cvSize(Binh2->width,Binh2->height),IPL_DEPTH_8U,1);
   Mask=cvCreateImage(cvSize(Binh2->width,Binh2->height),IPL_DEPTH_8U,1);
   cvNot(Binh2,NBinh2); 
   cvThreshold(NBinh2,Mask, 127, 1, CV_THRESH_BINARY);

    
   //Power-law transformation   
   Imgh1Proc=cvCreateImage(cvSize(Binh2->width,Binh2->height),IPL_DEPTH_8U,1);
   powerlaw(Imgh1, Imgh1Proc);

   
   //create Mask2 
   Mask2=cvCreateImage(cvSize(Binh2->width,Binh2->height),IPL_DEPTH_8U,1);
   cvMul(Imgh1Proc, Mask, Mask2, 1.0);
   cvAdd(Mask2,Binh2, Mask2, NULL);

   
   //create edges image
   Edges=cvCreateImage(cvSize(Binh2->width,Binh2->height),IPL_DEPTH_8U,1);
   cannythreshold(&higthreshold, &lowthreshold, &filterSze, NBinh2);
   cvCanny(NBinh2, Edges, lowthreshold, higthreshold, filterSze); 
   
   // Dilate Edge Image
   imdilate(Edges, Edges);
     
   //smooth mask2
   FMask2=cvCreateImage(cvSize(Binh2->width,Binh2->height),IPL_DEPTH_8U,1);
   cvSmooth(Mask2, FMask2, CV_GAUSSIAN, 5, 5, 1.0, 1.0);
      
   
   //create final HT12 image
   Imgh12=cvCreateImage(cvSize(Binh2->width,Binh2->height),IPL_DEPTH_8U,1);
   NEdges=cvCreateImage(cvSize(Binh2->width,Binh2->height),IPL_DEPTH_8U,1);
   cvNot(Edges,NEdges);
   cvThreshold(Edges,Edges, 127, 1, CV_THRESH_BINARY);
   cvThreshold(NEdges,NEdges, 127, 1, CV_THRESH_BINARY);
   cvMul(FMask2,Edges,FMask2, 1.0);  
   cvMul(Mask2,NEdges,Mask2, 1.0);   
   cvAdd(Mask2,FMask2, Imgh12, NULL);
   
      
   //clear other images and write a BMP version of final image
   strcpy(strdir,str2);
   ptr=strrchr(strdir,'\\');
   std::string delcmd;
   if (ptr == NULL) ptr=strrchr(strdir,'/');
   if (ptr != NULL) 
   {        //create name of final image
            strcpy(BmpFile,ptr+1);
            strcpy(BmpFile+(strlen(ptr)-8),"HT12.bmp");
            
            //delete unsuded images
            //strcpy(ptr+1,"*.brw"); /*Original*/
			strcpy(ptr + (strlen(ptr) - 4), "_output.brw"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
            system(delcmd.c_str()); /*3rd modification*/
			//system(strdel);
            
            //strcpy(ptr+1,"*.dm");
			strcpy(ptr + (strlen(ptr) - 4), ".dm"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
			system(delcmd.c_str()); /*3rd modification*/
			//system(strdel);
            
            //strcpy(ptr+1,"*.hcm");
			strcpy(ptr + (strlen(ptr) - 3), ".hcm"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
			system(delcmd.c_str()); /*3rd modification*/
			//system(strdel);
            
            //strcpy(ptr+1,"*.lcm");
			strcpy(ptr + (strlen(ptr) - 4), ".lcm"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
			system(delcmd.c_str()); /*3rd modification*/
			//system(strdel);
            
            //strcpy(ptr+1,"*.lfm");
			strcpy(ptr + (strlen(ptr) - 4), ".lfm"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
			system(delcmd.c_str()); /*3rd modification*/
            //system(strdel);
            
            //strcpy(ptr+1,"*.min");
			strcpy(ptr + (strlen(ptr) - 4), ".min"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
			system(delcmd.c_str()); /*3rd modification*/
            //system(strdel);
            
            //strcpy(ptr+1,"*.qm");
			strcpy(ptr + (strlen(ptr) - 4), ".qm"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
			system(delcmd.c_str()); /*3rd modification*/
            //system(strdel);
            
            //strcpy(ptr+1,"*.xyt");
			strcpy(ptr + (strlen(ptr) - 3), ".xyt"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
			system(delcmd.c_str()); /*3rd modification*/
            //system(strdel);
            
            //strcpy(ptr+1,"*.jpg");
			strcpy(ptr + (strlen(ptr) - 4), ".jpg"); /*1st modification*/
            strcpy(strdel,"del ");
            strcat(strdel, strdir);
			delcmd = std::string(strdel); /*2nd modification*/
			std::replace(delcmd.begin(), delcmd.end(), '/', '\\'); /*2nd modification*/
			system(delcmd.c_str()); /*3rd modification*/
            //system(strdel);
   }
   else
   {   //create name of final image
       strcpy(BmpFile,str2);
       strcpy(BmpFile+(strlen(str2)-7),"HT12.bmp");
       
       //delete unsuded images
       strcpy(strdel,"del *.brw *.dm *.hcm *.lcm *.lfm *.min *.qm *.xyt *.jpg"); 
       system(strdel);       
   }  
   //write final HT12 image
   
   if(!cvSaveImage(BmpFile, Imgh12, NULL)) printf("Could not save: %s\n",BmpFile);
   free(strdel);
   free(strdir);
   free(BmpFile);
   
   result = cv::cvarrToMat(Imgh12);
   
   cvReleaseImage( &FMask2);
   cvReleaseImage( &Imgh1Proc );
   cvReleaseImage( &Imgh1 );
   cvReleaseImage( &Imgh2 );
   cvReleaseImage( &Binh2 );
   cvReleaseImage( &NBinh2 );
   cvReleaseImage( &Mask );
   cvReleaseImage( &Mask2 );
   cvReleaseImage( &Edges );
   cvReleaseImage( &NEdges );     
   
   return result;
    
}

int powerlaw(IplImage *imgIN, IplImage *imgOUT){

   IplImage *imgIN32 = cvCreateImage(cvSize(imgIN->width, imgIN->height), IPL_DEPTH_32F, 1); 

   float gamma = 1.2;
   
   int i, j;

   cvConvertScale(imgIN, imgIN32, 1.0/255, 0.0 ); 
   
   // Power-law transformation
    for(i=0; i<imgIN32->height; i++ ){
             for(j=0; j<imgIN32->width; j++ ){                  
               ((float *)(imgIN32 ->imageData + i*imgIN32 ->widthStep) )[j] = 
                pow( ((float *)(imgIN32->imageData + i*imgIN32->widthStep) )[j], gamma);                              
             }                                            
    }     

   /*cvNormalize(imgIN32,imgIN32,0,255,CV_MINMAX);*/
   cvNormalize(imgIN32, imgIN32, 0, 255, CV_MINMAX, CV_8UC1);
   cvConvertScale(imgIN32, imgOUT, 1.0, 0.0);
   
   cvReleaseImage( &imgIN32 ); 
    
    return 0;
}

int imdilate(IplImage *ImgSegmentMask, IplImage *ImgSegmentMorph){

      int sizeSE = 2;
      int iterations = 1;
   	  IplConvKernel* structuringElement = cvCreateStructuringElementEx(sizeSE, sizeSE, cvFloor(sizeSE/2), cvFloor(sizeSE/2), CV_SHAPE_CROSS, NULL );
   	  cvDilate(ImgSegmentMask, ImgSegmentMorph, structuringElement, iterations);
   	  cvReleaseStructuringElement(&structuringElement); 
   	  
   	  return 0;    
}


int cannythreshold(float *higthreshold, float *lowthreshold, int *filterSzeOut, IplImage *NBinh2){

    // Determine filter size, based on sigma
   
   // cvCanny allowe filter of sizes  3, 5 or 7
   // if sigma == 0.5 => filterSze = 3
   // if sigma == 1.0 => filterSze = 5
   // if sigma == 1.5 => filterSze = 7
   // sigma == 1.5 results in pour edge detection, best is sigma == 1.0
    
    float pw[MaxFiltSize], sigma = 1.0, ssq, GaussianDieOff =  0.0001, ThresholdRatio = 0.4, aux[MaxFiltSize], t;
    float **x, **y, PercentOfPixelsNotEdges = 0.7;
    CvMat *gau, *gaut, *dgau2D, *dgau2Dt, *mag;
    IplImage *ax, *ay, *NBinh2Smooth;

    int i, j, filterSze;
       
    for(i = 0; i<MaxFiltSize;i++) pw[i] = (float)i+1;
    ssq = pow(sigma,2);
    
    for(i = 0; i<MaxFiltSize;i++) 
          aux[i] = exp(-(pw[i]*pw[i])/(2*ssq));
    
    for(i = 0; i<MaxFiltSize;i++) 
          if (aux[i] > GaussianDieOff)
               filterSze = i+1;
               
    if( filterSze%2 == 0 && filterSze != 0) filterSze++;

    if (filterSze==0) filterSze = 3;
    
    *filterSzeOut = filterSze;
                     
    //  1D gaussian filter
    gau = cvCreateMat(1, 2*filterSze+1, CV_32F );
    float *gaudata = gau->data.fl;
    for(j = -filterSze; j<=filterSze; j++) {              
       t = (double)j;
       gaudata[j+filterSze] =  exp( -(float)(t*t)/(2*ssq))/(2*3.1415*ssq); 
   }
   
    // Transposed 1D gaussian filter
    gaut = cvCreateMat(2*filterSze+1, 1, CV_32F );
    cvTranspose(gau, gaut);    

    // Meshgrid    
    x      = (float **)malloc( (2*filterSze+1)*sizeof(float *) );
    y      = (float **)malloc( (2*filterSze+1)*sizeof(float *) );
    for(i = 0; i<2*filterSze+1; i++){
          *(x+i)      = (float*)malloc((2*filterSze+1)*sizeof(float));
          *(y+i)      = (float*)malloc((2*filterSze+1)*sizeof(float));
    }
    for(i = -filterSze; i<=filterSze; i++)    
      for(j = -filterSze; j<=filterSze; j++){
        x[i+filterSze][j+filterSze] = j;
        y[i+filterSze][j+filterSze] = i;
      }
  
   /* Find the directional derivative of 2D Gaussian (along X-axis)
      Since the result is symmetric along X, we can get the derivative along
      Y-axis simply by transposing the result for X direction.
   */   
    dgau2D = cvCreateMat(2*filterSze+1, 2*filterSze+1, CV_32F );
    float *dgau2Ddata = dgau2D->data.fl;

    for(i = 0; i < 2*filterSze+1; i++) 
      for(j = 0; j < 2*filterSze+1; j++)
              dgau2Ddata[i*dgau2D->cols+j] = -x[i][j]*exp( -(pow(x[i][j],2)+pow(y[i][j],2)) / (2* ssq) ) / (3.1415*ssq);                            

    dgau2Dt = cvCreateMat(2*filterSze+1, 2*filterSze+1, CV_32F );
    cvTranspose(dgau2D, dgau2Dt);
                           
    // Smooth the image out
    IplImage *NBinh232 = cvCreateImage(cvSize(NBinh2->width, NBinh2->height), IPL_DEPTH_32F, 1); 
    cvConvertScale(NBinh2, NBinh232, 1.0/255, 0.0 ); 
        
    NBinh2Smooth=cvCreateImage(cvSize(NBinh2->width,NBinh2->height),IPL_DEPTH_32F, 1);
    cvFilter2D(NBinh232, NBinh2Smooth, gau, cvPoint(-1, -1));
    cvFilter2D(NBinh2Smooth, NBinh2Smooth, gaut, cvPoint(-1, -1));
    
    cvReleaseMat( &gau );
    cvReleaseMat( &gaut );
    cvReleaseImage( &NBinh232 );
              
    // Apply directional derivatives   
    ax = cvCreateImage(cvSize(NBinh2->width,NBinh2->height),IPL_DEPTH_32F,1);
    cvFilter2D(NBinh2Smooth, ax, dgau2D, cvPoint(-1,-1));           
    ay=cvCreateImage(cvSize(NBinh2->width,NBinh2->height),IPL_DEPTH_32F,1);   
    cvFilter2D(NBinh2Smooth, ay, dgau2Dt, cvPoint(-1,-1));
    
    cvReleaseMat( &dgau2D );
    cvReleaseMat( &dgau2Dt );
    cvReleaseImage( &NBinh2Smooth );
                       
    // Calculate magnitude                                
    mag = cvCreateMat(NBinh2->height, NBinh2->width, CV_32F );
    float *data_mag = mag->data.fl;
     
    for (i = 0; i < ax->height; i++)
         for (j = 0; j < ax->width; j++){
                    data_mag[i*mag->cols+j] = sqrt(                    
                    pow( ((float *)(ax ->imageData + i*ax ->widthStep) )[j],2) +
                    pow( ((float *)(ay ->imageData + i*ay ->widthStep) )[j],2) );
      }
      
      
      
    cvReleaseImage( &ax );
    cvReleaseImage( &ay );

    // Find maximum magnitude
     float maxmag = 0;
     for (i = 0; i < NBinh2->height; i++)
          for (j = 0; j < NBinh2->width; j++)
                    if(data_mag[i*mag->cols+j] > maxmag)   
                               maxmag =   data_mag[i*mag->cols+j];
                               
     // Create image structure and scale magnitudes to fit into the interval [0, 255];
     IplImage *imgmag = cvCreateImage(cvSize(NBinh2->width,NBinh2->height),IPL_DEPTH_8U,1);
     
     
     if(maxmag>0){
     for (i = 0; i < imgmag->height; i++)
          for (j = 0; j < imgmag->width; j++)
                              //data_mag[i*mag->cols+j] = 255*data_mag[i*mag->cols+j]/maxmag;
                              ((uchar *)(imgmag ->imageData + i*imgmag ->widthStep) )[j] = (uchar)(255*data_mag[i*mag->cols+j]/maxmag);
     }

     //Compute histogram                   
     int histsize = 64;
     float histrange[] = {0, 255};
     float *ranges[] = {histrange};
     float count, cumsum[65];
     
    cumsum[0] = 0;

    CvHistogram* histmag = cvCreateHist(1, &histsize, CV_HIST_ARRAY, ranges, 1);          
    cvCalcHist(&imgmag, histmag, 0, NULL);
     
    for(i = 0; i < histsize; i++ ){
          count = cvQueryHistValue_1D(histmag,i);
          cumsum[i+1] = cumsum[i]+count;
    }

     
    for(i = 1; i < histsize+1; i++ ){
        if(cumsum[i] > PercentOfPixelsNotEdges*imgmag->height*imgmag->width) break;
    }
       
    *higthreshold = 255.0*i/64;
    *lowthreshold = 0.4*(*higthreshold);

    cvReleaseMat  ( &mag );
    cvReleaseImage( &imgmag );
    cvReleaseHist ( &histmag );

    return 0;
    
}

