#include "mat2matrix.h"

Matrix Mat2Matrix::mat2Matrix(cv::Mat mat){

    int rows = mat.rows;
    int cols = mat.cols;

    Matrix m(rows,cols);


    for (int i = 0 ; i < rows; ++i){
        for (int j = 0 ; j < cols; ++j){
            m.set(i,j,mat.at<float>(i,j));
        }
    }

    return m;

}

cv::Mat Mat2Matrix::matrix2Mat(Matrix mat){

    cv::Mat m(mat.rows(), mat.cols(), CV_32FC1);
    for (int i = 0 ; i < mat.rows(); ++i){
        for (int j = 0 ; j < mat.cols(); ++j){
            m.at<float>(i,j)= mat.get(i,j);
        }
    }

    return m;

}