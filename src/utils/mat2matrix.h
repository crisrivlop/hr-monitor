
#include <opencv2/core.hpp>
#include "opencv2/opencv.hpp"
#include "../maths/matrix.h"

#ifndef MAT_2_MATRIX_H_
#define MAT_2_MATRIX_H_

class Mat2Matrix {
    public:
        static Matrix mat2Matrix(cv::Mat mat);
        static cv::Mat matrix2Mat(Matrix mat);
};



#endif //MAT_2_MATRIX_H_