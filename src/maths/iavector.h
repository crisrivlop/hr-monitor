
#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <iostream>
using namespace std;

class Matrix;

float dotHelper( float *a, float *b, unsigned int len);

class Vector {
        friend Matrix;
        float * _data = 0;
        unsigned int _size;
        Vector add(const Vector& b);
        Vector add(const float& b);
        void  addi(const Vector& b);
        Vector sub(const Vector& b);
        Vector sub(const float& b);
        Vector mul(const Vector& b);
        Vector mul(const float& b);
        Vector div(const Vector& b);
        Vector div(const float& b);
        static ostream& print(ostream& os, Vector const & v);
        void clone(const Vector& b);

    public:
        Vector(unsigned int size);
        Vector(unsigned int size,const float * input);
        Vector(const Vector &v);
        Vector();
        
        float& operator[](unsigned int i) {return this->_data[i];};
        Vector subset(unsigned int i, unsigned int f);
        Vector operator+(const Vector& b){return this->add(b);};
        Vector operator+(const float& b){return this->add(b);};
        void   operator=(const Vector& b){this->clone(b);};
        void   operator+=(const Vector& b){this->addi(b);};
        Vector operator-(const Vector& b){return this->sub(b);};
        Vector operator-(const float& b){return this->sub(b);};
        Vector operator*(const Vector& b){return this->mul(b);};
        Vector operator*(const float& b){return this->mul(b);};
        Vector operator/(const Vector& b){return this->div(b);};
        Vector operator/(const float& b){return this->div(b);};
        bool operator==(const Vector& b){
            bool r = true; 
            for (int i = 0; i < this->_size; i++ ) 
                r &= b._data[i] == this->_data[i];
            return r;
        };
        friend ostream& operator<<(ostream& os, Vector const & v){return v.print(os,v);};

        float dot(const Vector& b);

        Vector inv();
        Vector exponential();

        Vector absolute();

        float maximum();
        
        Matrix extern_product(const Vector& b);

        Matrix toMatrix();

        Vector append(const float& e);

        void appendInPlace(const float& e);

        unsigned int len() const;


        Vector relu();
        Vector relu_derivated();

        ~Vector();
};

#endif //_VECTOR_H_