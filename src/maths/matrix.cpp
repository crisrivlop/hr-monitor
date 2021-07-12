#include "matrix.h"
#include "iavector.h"
#include <iostream>


#if defined(HAVE_NEON)

#include <arm_neon.h>
#include <omp.h>

#endif


Matrix::Matrix(unsigned int rows, unsigned int cols){

    this->_rows = rows;
    this->_cols = cols;
    this->_data = new float*[rows]();
    for (int i = 0; i < rows; ++i)
        this->_data[i] = new float[cols];
}



Matrix::Matrix() : _rows(0),_cols(0),_data(0){
    this->_data = 0;
}

Matrix::Matrix(const Matrix& m){
    this->_data = 0;
    this->_rows = 0;
    this->clone(m);

}

void Matrix::clone(const Matrix& m){
    for (int i = 0; i < this->_rows; ++i)
        delete [] this->_data[i];
    delete [] this->_data;

    this->_rows = m._rows;
    this->_cols = m._cols;
    this->_data = new float*[m._rows]();
    for (int i = 0; i < m._rows; i++){
        this->_data[i] = new float[m._cols];
        for(int j = 0; j < m._cols; j++)
            this->_data[i][j] = m._data[i][j];
    }
}



Matrix::~Matrix(){
    for (int i = 0; i < this->_rows; ++i){
        delete [] this->_data[i];
        this->_data[i]=0;
    }        
    delete [] this->_data;
};


float ** Matrix::data(){
    return this->_data;
};


Matrix Matrix::add(const Matrix& b){
    Matrix m{this->_rows,this->_cols};
#if defined(HAVE_NEON)
    unsigned int cols4 = (this->_cols<4)? 0 : this->_cols - this->_cols%4;
    for (unsigned int i = 0; i < this->_rows; i++){
        int j = 0;
#pragma omp parallel private(j) num_threads(4) 
        {
            #pragma omp for
            for (j = 0; j < cols4; j+=4){
                vst1q_f32(m._data[i] + j,  vaddq_f32( vld1q_f32(this->_data[i] + j), vld1q_f32(b._data[i] + j)));
            }

            for (j = cols4; j < this->_cols; j++){
                m._data[i][j] = this->_data[i][j] + b._data[i][j];
            }
        }
    }


#else
    for (int i = 0; i < this->_rows; i++)
        for (int j = 0; j < this->_cols; j++)
            m._data[i][j] = this->_data[i][j] + b._data[i][j];
#endif

    return m;
}


Matrix Matrix::sub(const Matrix& b){

    Matrix m{this->_rows,this->_cols};

#if defined(HAVE_NEON)
    unsigned int cols4 = (this->_cols<4)? 0 : this->_cols - this->_cols%4;
    float * rowThis, * rowB;
    unsigned int i;
#pragma omp parallel private(i,rowThis,rowB) num_threads(4) 
{
#pragma omp for
    for (i = 0; i < this->_rows; i++){
        int j = 0;
        rowThis = this->_data[i];
        rowB =  b._data[i];
        for (j = 0; j < cols4; j+=4){
            vst1q_f32(m._data[i] + j,  vsubq_f32( vld1q_f32(rowThis + j), vld1q_f32(rowB + j)));
        }
        for (j = cols4; j < this->_cols; j++){
            m._data[i][j] = this->_data[i][j] - b._data[i][j];
        }
    }
}


#else
    for (int i = 0; i < this->_rows; i++)
        for (int j = 0; j < this->_cols; j++)
            m._data[i][j] = this->_data[i][j] - b._data[i][j];
#endif


    
    return m;
}


Matrix Matrix::mul(const float b){
    Matrix m{this->_rows,this->_cols};

#if defined(HAVE_NEON)
    unsigned int cols4 = (this->_cols<4)? 0 : this->_cols - this->_cols%4;
    float32x4_t constant = vdupq_n_f32(b);
    unsigned int i;
#pragma omp parallel private(i) num_threads(4) 
{
    #pragma omp for
    for (i = 0; i < this->_rows; i++){
        int j = 0;
            for (j = 0; j < cols4; j+=4){
                vst1q_f32(m._data[i] + j,  vmulq_f32( vld1q_f32(this->_data[i] + j), constant ) );
            }

            for (j = cols4; j < this->_cols; j++){
                m._data[i][j] = this->_data[i][j] * b;
            }
    }
}

#else
    for (int i = 0; i < this->_rows; i++)
        for (int j = 0; j < this->_cols; j++)
            m._data[i][j] = this->_data[i][j] * b;
#endif

    return m;
}

Vector Matrix::dot(const Vector& b){
    unsigned int len = this->rows();
    Vector v{len};
    Vector tmp{b.len()};

#if defined(HAVE_NEON)
    float *acum = new float[4];
    float32x4_t acumulator; 
    int j = 0;
    for (int i = 0; i < this->_rows; i++){
        v[i] =dotHelper( this->_data[i], b._data, b._size);
    }

#else
    for (int i = 0; i < this->_rows; i++){
        for (int j = 0; j < this->_cols; j++)
            tmp[j] = this->_data[i][j];
        v[i] = (tmp.dot(b));
    }
#endif
    return v;
}



ostream& Matrix::print(ostream& os, Matrix const & v){
    os << std::endl << v.rows() << "x" << v.cols()<< "\n[";
    for (int i = 0; i < v._rows; i++){
        if(i != 0)os << " ";
        for (int j = 0; j < v._cols; j++){
            os << v._data[i][j];
            if(j+1 != v._cols)os << ",\t";
        }

        if(i+1 == v._rows)os << "]";
        
        os << std::endl;
    }
    return os;
}

unsigned int Matrix::cols() const{
    return this->_cols;
}
unsigned int Matrix::rows() const{
    return this->_rows;
}


Matrix Matrix::transpose(){
    Matrix m{_cols,_rows};


#if defined(HAVE_NEON)
    unsigned int cols4 = (this->_cols<4)? 0 : this->_cols - this->_cols%4;
    for (unsigned int i = 0; i < this->_rows; i++){
        int j = 0;
#pragma omp parallel private(j) num_threads(4) 
        {
            #pragma omp for
            for (j = 0; j < cols4; j+=4){
                vst1q_f32(m._data[j] + i ,vld1q_f32(this->_data[i] + j));
            }

            for (j = cols4; j < this->_cols; j++){
                m._data[j][i] = this->_data[i][j];
            }
        }
    }


#else
    for (int i = 0; i < this->_rows; i++)
        for (int j = 0; j < this->_cols; j++){
            m._data[j][i] = this->_data[i][j];
        }
#endif
    
    return m;
}

Matrix Matrix::conv(const Matrix &b, bool isValid){
    unsigned int resultRows = _rows + ((isValid)?(-(int)(b._rows)+1):b._rows  - 1 );
    unsigned int resultCols = _cols + ((isValid)?-(int)(b._cols)+1:b._cols - 1);
    Matrix m{resultRows, resultCols};
    int rowOffset = ((isValid)?0:1-b._rows);
    int colOffset = ((isValid)?0:1-b._cols);
    
    for(int i = 0; i < resultRows; i++){
       for(int j = 0; j < resultCols; j++){
           float res = 0;
           for(int k = 0; k < b._rows; k++){
               for(int l = 0; l < b._cols; l++){
                   int x = i+rowOffset+k;
                   int y = j+colOffset+l;
                   if(0 <= x && x < _rows && 0 <= y && y < _cols){
                        res += _data[x][y]*b._data[k][l];
                   }
               }
           }
           m._data[i][j] = res;
       } 
    }
    return m;
}

void Matrix::fill(float number){
    for (int i = 0; i < this->_rows; i++)
        for (int j = 0; j < this->_cols; j++){
            this->_data[i][j] = number;
        }
}