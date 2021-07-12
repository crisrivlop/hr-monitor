#include "size.h"

MatrixSize::MatrixSize(int width, int height){
    this->_width = width;
    this->_height = height;
}
int MatrixSize::width(){
    return this->_width;
}
int MatrixSize::height(){
    return this->_height;
}