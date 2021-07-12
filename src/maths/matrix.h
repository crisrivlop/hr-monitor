
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include "iavector.h"
using namespace std;


class Matrix {
        friend Vector;
        float ** _data = 0;
        unsigned int _rows;
        unsigned int _cols;
        Matrix add(const Matrix& b);
        Matrix sub(const Matrix& b);
        Matrix mul(float b);
        static ostream& print(ostream& os, Matrix const & v);
        void clone(const Matrix& b);

    public:
        Matrix(unsigned int rows, unsigned int cols);
        Matrix(const Matrix& m);
        Matrix();
        
        float& get(unsigned int i,unsigned int j) {return this->_data[i][j];};
        void set(unsigned int i,unsigned int j,float data) {this->_data[i][j] = data;};
        Matrix operator+(const Matrix& b){return this->add(b);};
        Matrix operator-(const Matrix& b){return this->sub(b);};
        Matrix operator*(float b){return this->mul(b);};
        void   operator=(const Matrix& b){this->clone(b);};
        friend ostream& operator<<(ostream& os, Matrix const & m){return m.print(os,m);};

        Vector dot(const Vector& b);
        Matrix transpose();
        Matrix conv(const Matrix &b, bool isValid = true);

        void fill(float number);


        unsigned int cols() const;
        unsigned int rows() const;

        float ** data();

        ~Matrix();
};

#endif //_MATRIX_H_