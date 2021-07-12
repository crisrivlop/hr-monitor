

#ifndef _PYRAMID_H_
#define _PYRAMID_H_

#include <iostream>
using namespace std;
#include "../maths/matrix.h"
#include "../maths/iavector.h"
#include "../utils/size.h"

class Pyramid {

    Vector namedFilter(string filtername);
    int maxPyramidHeight(int rows, int cols,int filterSize);
    
    

    Pyramid add(const Pyramid& b);
    Pyramid sub(const Pyramid& b);
    Pyramid mul(float b);
    void clone(const Pyramid& b);

    void BuildLpyrAux(Matrix im, int height, Matrix f1, Matrix f2, string edges);
    void reconLpyrAux(int level, int last, Matrix filter, string edges);

    Pyramid(int l);

    public:
        Matrix * _pyramid;
        int _length;
        Matrix rconv2(Matrix a, Matrix b, int center = 0);
        Matrix corrDn(Matrix im, Matrix filter, string edges, int* step = NULL, int* start = NULL, int* stop = NULL);
        Matrix upConv(Matrix im, Matrix filter, string edges, MatrixSize s, int* step = NULL, int* start = NULL, int* stop = NULL, Matrix* res = NULL);
    
        Pyramid();
        Pyramid(const Pyramid& m);
        void BuildLpyr(Matrix im, string filter1 = "", string filter2="", string edges="");
        Matrix reconLpyr();
        Pyramid operator+(const Pyramid& b){return this->add(b);};
        Pyramid operator-(const Pyramid& b){return this->sub(b);};
        Pyramid operator*(float b){return this->mul(b);};
        void   operator=(const Pyramid& b){this->clone(b);};
        void IgnoreHighestAndLowestFrequencies();


        ~Pyramid();
};

#endif //_PYRAMID_H_