#include <iostream>
#include "maths/iavector.h"
#include "maths/matrix.h"
#include <opencv2/core.hpp>
#include "opencv2/opencv.hpp"
#include <chrono>
#include <vector>
#include "./utils/mat2matrix.h"
#include "pyramid/pyramid.h"
#include "cmath"

using namespace cv;


int main(int args, char ** argv){


    VideoCapture cap(argv[1]); 
    setUseOptimized(false);
    for (int i = 0; i < args ; i++){
        string arg = argv[i];
        if (arg.compare("-o")==0){
            std::cout << "optimized!" << std::endl;
            setUseOptimized(true);
        }
    }

    int framecount = cap.get(cv::CAP_PROP_FRAME_COUNT);
    int height = cap.get(cv::CAP_PROP_FRAME_HEIGHT);
    int width = cap.get(cv::CAP_PROP_FRAME_WIDTH);
    int number_of_channels = 3;
    int frame_rate = cap.get(cv::CAP_PROP_FPS);

    int ex = static_cast<int>(cap.get(CAP_PROP_FOURCC));

    Size size = Size(width,height);

    std::cout << "Frame count: " << framecount << std::endl;
    std::cout << "H: " << height << std::endl;
    std::cout << "W: " << width << std::endl;
    std::cout << "frame rate: " << frame_rate << std::endl;
    // Check if camera opened successfully
    if(!cap.isOpened()){
        cout << "Error opening video stream or file" << endl;
        return -1;
    }

    VideoWriter outputVideo;
    string posfix = "out.avi";
    outputVideo.open(argv[1] + posfix, ex, frame_rate,size,true);
    if (!outputVideo.isOpened())
    {
        cout  << "Could not open the output video for write: " << endl;
        return -1;
    }

    Mat src;
    Pyramid pyrI;
    Pyramid pyrY;
    Pyramid pyrQ;

    Pyramid pI;
    Pyramid pY;
    Pyramid pQ;
    Matrix y_layer;
    Matrix i_layer;
    Matrix q_layer;

    Pyramid lowpass1I;
    Pyramid lowpass2I;
    Pyramid lowpass1Y;
    Pyramid lowpass2Y;
    Pyramid lowpass1Q;
    Pyramid lowpass2Q;

    //cout << "Begin" << endl;

    for (int i = 0; i < framecount-10; i++){
        cap >> src;

        cout << "Total: " << (i*100)/(framecount-10) << endl;
        
        //RGB TO NTCS
        // let's define the matrix for RGB -> YIQ conversion
        // Matx33f matYIQ( 0.299f,  0.587f,  0.114f,
        //                 0.596f, -0.274f, -0.322f,
        //                 0.211f, -0.523f,  0.312f);

        Matx33f matYIQ( 0.299f,  0.587f,  0.114f,
                        0.596f, -0.274f, -0.322f,
                        0.211f, -0.523f,  0.312f 
                        );
        // let's define the matrix for YIQ -> RGB conversion
        ///https://la.mathworks.com/help/images/ref/ntsc2rgb.html
        Matx33f matRGB( 1.000f,  0.956f,  0.621,
                        1.000f, -0.272f, -0.647f,
                        1.000f, -1.106f,  1.703f);

        // I assume you have a source image of type CV_8UC3 
        // CV_8UC3: 3 channels, each on unsigned char, so [0,255]
        // here is a dummy one, black by default, 256x256
        //Mat ImgRGB_8UC3(256, 256, CV_8UC3);
        // We need to convert this to a new image of type CV_32FC3
        // CV_32FC3: 3 channels each on 32bit float [-inf, +inf]
        // we need to do this because YIQ result will be in [-1.0, 1.0] (I & Q)
        // so this obviously cannot be stored in CV_8UC3
        // At the same time, we will also divide by 255.0 to put values in [0.0, 1.0]
        Mat ImgYIQ_32FC3(src);
        src.convertTo(ImgYIQ_32FC3, CV_32FC3, 1.0/255.0);
        // at this point ImgYIQ_32FC3 contains pixels made of 3 RGB float components in [0-1]
        // so let's convert to YIQ
        // (cv::transform will apply the matrix to each 3 component pixel of ImgYIQ_32FC3)
        cv::transform(ImgYIQ_32FC3, ImgYIQ_32FC3, matYIQ);

        Mat frameNtsc = ImgYIQ_32FC3;

        std::vector<Mat> yiq;

        cv::split(frameNtsc,yiq);

        y_layer = Mat2Matrix::mat2Matrix(yiq[0]);
        i_layer = Mat2Matrix::mat2Matrix(yiq[1]);
        q_layer = Mat2Matrix::mat2Matrix(yiq[2]);

        

        float chromAttenuation=0.1;
        float alpha = 10;
        float lambda_c = 16;
        float delta = lambda_c/8.0/(1.0+alpha);
        float exaggeration_factor = 5;
        float lambda = pow((std::pow(height,2)+ pow(width,2)),0.5)/3.0;
        //cout << "Build L Pyr" << endl;
        pI.BuildLpyr(i_layer);
        pY.BuildLpyr(y_layer);
        pQ.BuildLpyr(q_layer);

        if (i == 0){
            lowpass1I = pI;
            lowpass2I = pI;
            lowpass1Y = pY;
            lowpass2Y = pY;
            lowpass1Q = pQ;
            lowpass2Q = pQ;
        }
        

        if (i > 0){
            //cout << "Filtering" << endl;
            float r1= 0.5;
            float r2= 0.01;
            float rr1 = (1.0-r1);
            float rr2 = (1.0-r2);

            //cout << "Filtering 2" << endl;

            lowpass1I = lowpass1I*(1.0-r1) + pI*r1;
            lowpass2I = lowpass2I*(1.0-r2) + pI*r2;
            lowpass1Y = lowpass1Y*(1.0-r1) + pY*r1;
            lowpass2Y = lowpass2Y*(1.0-r2) + pY*r2;
            lowpass1Q = lowpass1Q*(1.0-r1) + pQ*r1;
            lowpass2Q = lowpass2Q*(1.0-r2) + pQ*r2;
            Pyramid filteredI = lowpass1I - lowpass2I;
            Pyramid filteredY = lowpass1Y - lowpass2Y;
            Pyramid filteredQ = lowpass1Q - lowpass2Q;


            // filteredI.IgnoreHighestAndLowestFrequencies();
            // filteredY.IgnoreHighestAndLowestFrequencies();
            // filteredQ.IgnoreHighestAndLowestFrequencies();
            for(int i=pI._length-1; i >= 0;i--){
                float currAlpha = lambda/delta/8.0 - 1.0;
                currAlpha = currAlpha*exaggeration_factor;
                if (i==pI._length-1 || i==0 ){
                    filteredI._pyramid[i].fill(0);
                    filteredY._pyramid[i].fill(0);
                    filteredQ._pyramid[i].fill(0);
                }
                if(currAlpha > alpha){
                    //cout << "ON alpha" << endl;
                    filteredI._pyramid[i] = filteredI._pyramid[i] * alpha; 
                    filteredY._pyramid[i] = filteredY._pyramid[i] * alpha; 
                    filteredQ._pyramid[i] = filteredQ._pyramid[i] * alpha; 
                }
                else{
                    //cout << "ON Curr alpha" << endl;
                    filteredI._pyramid[i] = filteredI._pyramid[i] * currAlpha; 
                    filteredY._pyramid[i] = filteredY._pyramid[i] * currAlpha; 
                    filteredQ._pyramid[i] = filteredQ._pyramid[i] * currAlpha; 
                }
                lambda = lambda/2.0;             
            }
            //cout << "Reconstruct" << endl;
            Matrix outputY  = filteredY.reconLpyr();
            Matrix outputI  = filteredI.reconLpyr()*chromAttenuation;
            Matrix outputQ  = filteredQ.reconLpyr()*chromAttenuation;
            //cout << y_layer.rows() << "x" << y_layer.cols() << " <---> " << outputY.rows() << "x" << outputY.cols() << endl;
            //cout << i_layer.rows() << "x" << i_layer.cols() << " <---> " << outputI.rows() << "x" << outputI.cols() << endl;
            //cout << q_layer.rows() << "x" << q_layer.cols() << " <---> " << outputQ.rows() << "x" << outputQ.cols() << endl;
            //cout << outputY << endl;
            //exit(0);
            y_layer = y_layer + outputY;
            i_layer = i_layer + outputI;
            q_layer = q_layer + outputQ;


            //cout << "Override" << endl;

            //cout << endl << rr1 << endl;
            //cout << rr2 << endl;

            // cout << pyr*r1  << endl;
            // cout << pyr*r2  << endl;
            //cout << p._pyramid[6]  << endl;
            //cout << filtered._pyramid[6]  << endl;
            //exit(0);
        }

        //end

        pyrI = pI;
        pyrY = pY;
        pyrQ = pQ;
        
        
        //cout << filtered  << endl;

        
        
        //build lpyr

        Mat y_layer_ret = Mat2Matrix::matrix2Mat(y_layer);
        Mat i_layer_ret = Mat2Matrix::matrix2Mat(i_layer);
        Mat q_layer_ret = Mat2Matrix::matrix2Mat(q_layer);

        
        
        std::vector<Mat> yiq_ret;

        yiq_ret.push_back(y_layer_ret);
        yiq_ret.push_back(i_layer_ret);
        yiq_ret.push_back(q_layer_ret);
        

        Mat ImgYIQ_32FC3_2;
        cv::merge(yiq_ret, ImgYIQ_32FC3_2);
        cv::transform(ImgYIQ_32FC3_2, ImgYIQ_32FC3_2, matRGB);

        Mat mask = ImgYIQ_32FC3_2 < 0;
        ImgYIQ_32FC3_2.setTo(0,mask);
        mask = ImgYIQ_32FC3_2 > 1;
        ImgYIQ_32FC3_2.setTo(1,mask);

        ImgYIQ_32FC3_2.convertTo(ImgYIQ_32FC3_2, CV_8UC3,255);

        outputVideo << ImgYIQ_32FC3_2;

        

    }

    return 0;
}



// int main(){
//     Pyramid p;
//     Matrix m(500,900);
//     for (int i = 0;i<100000;i++){
//         m = Matrix(500,900);
//         Matrix m2(200,200);
//         p.BuildLpyr(m);
//     }        

// }


// void myf();
// int main()
// {
//     Matrix a(5, 6);
//     Matrix b(5, 5);
//     for (int i = 0; i < 30; i++)
//     {
//         a.set(i / 6, i % 6, i + 1);
//     }
//     for (int i = 0; i < 25; i++)
//     {
//         b.set(i / 5, i % 5, 0);
//     }
//     b.set(2, 2, 1);
//     cout << a << endl;
//     cout << b << endl;
//     Matrix c = a.conv(b);
//     Matrix d = a.conv(b, false);
//     cout << c << endl;
//     cout << d << endl;

//     Pyramid p;
//     Matrix *m = p.rconv2(a, b);
//     cout << *m << endl;
//     myf();
// }

void myf(){
    float a[9][8] = {24.113,   27.352,   27.712,   26.047,   24.803,   24.051,   24.623,   24.939
    ,28.037,   34.556,   37.659,   34.699,   30.418,   26.335,   25.610,   25.242
    ,33.025,   39.671,   41.303,   41.859,   42.040,   38.367,   36.576,   34.452
    ,36.542,   38.518,   33.377,   41.719,   49.070,   47.973,   46.127,   44.155
    ,39.744,   33.747,   22.719,   42.069,   52.855,   52.923,   50.489,   49.400
    ,39.809,   31.849,   19.011,   37.002,   50.942,   52.567,   50.067,   51.406
    ,44.523,   41.495,   31.673,   34.011,   42.170,   49.433,   46.161,   49.521
    ,56.054,   55.105,   47.349,   42.915,   39.292,   42.687,   45.104,   50.093
    ,61.637,   62.583,   60.294,   57.213,   48.799,   43.976,   47.547,   52.640};

    float d[9][15] = {
   14.9794, 17.7576,   20.4335,   20.3034,   19.6268,   18.9349,   18.3787,   17.8871,   17.4965,   17.2371,   16.7603,   16.8638,   17.6457,   17.7549,   17.4699,
   17.6406, 19.9103,   26.0388,   27.6412,   26.8718,   25.8979,   24.6391,   23.1764,   21.5740,   19.8451,   17.9702,   17.6638,   18.5126,   18.2303,   17.1185,
   20.2588, 24.0065,   30.0178,   30.9747,   28.5375,   28.0868,   29.8015,   30.9240,   30.1853,   28.6110,   26.5748,   25.8595,   26.5559,   25.6313,   21.9356,
   23.1352, 26.6985,   30.5160,   25.7892,   19.5990,   24.0095,   30.3146,   34.7600,   35.4294,   34.8680,   33.7632,   33.0995,   32.8711,   32.3175,   29.2120,
   26.8665, 28.7749,   29.1284,   14.1438,    8.4955,   21.9408,   32.5999,   36.5798,   37.7802,   38.6235,   37.8948,   35.8721,   35.6291,   35.5466,   33.8781,
   27.3237, 28.9374,   27.4726,   11.5392,    6.2696,   18.8580,   28.4022,   32.3803,   36.9888,   39.5555,   38.2637,   33.7198,   35.0506,   36.6577,   36.3707,
   30.0582, 32.2407,   32.7202,   23.7628,   19.0650,   22.7840,   25.0425,   23.8098,   29.0925,   36.0888,   37.9065,   32.1262,   29.8773,   34.9334,   36.8411,
   38.6526, 40.0477,   40.9420,   36.7280,   32.0308,   30.9868,   31.7050,   28.4294,   25.6028,   28.4087,   31.8618,   30.7553,   29.5237,   35.1210,   37.7869,
   42.7221, 43.8497,   45.1041,   44.2815,   42.1774,   41.2894,   41.7009,   39.0639,   33.4971,   30.6030,   30.7448,   31.1594,   32.5120,   37.1468,   38.8923
    };
   
    
    float b[5] = {0.088388,0.353553,0.530330,0.353553,0.088388};

    Matrix im(9,15);
    Matrix m(9,8);
    Matrix m2(Vector(5,b).toMatrix().transpose());

    for (int i = 0; i < 9;i++){
        for (int  j = 0; j < 8;j++){
            float f = a[i][j];
            m.set(i,j,f);
        }
    }

    for (int i = 0; i < 9;i++){
        for (int  j = 0; j < 15;j++){
            float f = d[i][j];
            im.set(i,j,f);
        }
    }

    int index21[2] = {2,1};
    int index12[2] = {1,2};
    int index11[2] = {0,0};
    Pyramid p;

    Matrix lo = p.corrDn(im,m2.transpose(),"",index12, index11);
    cout << lo << endl;
    Matrix lo2 = p.corrDn(m,m2,"",index21, index11);
    cout << lo2 << endl;
}