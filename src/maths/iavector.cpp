#include "iavector.h"
#include "matrix.h"

#include <exception>
#include <cmath>

 #ifdef __ARM_NEON

#include <arm_neon.h>

#endif

using namespace std;

Vector::Vector(unsigned int size){
    this->_size = size;
    this->_data = new float[size];
};

Vector::Vector(unsigned int size,const float * input){
    this->_size = size;
    this->_data = new float[size];

 #ifdef __ARM_NEON
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(this->_data + i, vld1q_f32(input + i) );
    }
    for (; i < this->_size; i++){
        this->_data[i] = input[i];
    }
#else
    for (int i = 0; i < this->_size; i++){
        this->_data[i] = input[i];
    }
#endif

}


Vector::Vector() : _size(0),_data(0){
};


Vector::~Vector(){
    delete [] this->_data;
};

unsigned int Vector::len() const{
    return this->_size;
};

Vector Vector::subset(unsigned int i, unsigned int f){
    if (f < i || f > this->_size){
        throw new exception();
    }
    Vector v{f-i+1};
 #ifdef __ARM_NEON
    int j = 0;
    if (this->_size>=4)
    for (; j  + 4 <= f-i; j+=4){
        vst1q_f32(v._data + j, vld1q_f32(this->_data + i + j) );
    }
    for (; j <= f-i; j++){
        v._data[j] = this->_data[i+j];
    }
#else
    for (int j = 0; j <= f-i; j++){
        v._data[j] = this->_data[i+j];
    }
#endif
    return v;
}


Vector Vector::add(const Vector& b){
    Vector v{this->_size};
#ifdef __ARM_NEON
    int i = 0;
    int j = 0;
    if (this->_size>=4)
    for (; i + 4 <= this->_size; i += 4){
        vst1q_f32(v._data + i,  vaddq_f32( vld1q_f32(this->_data + i), vld1q_f32(b._data + i)));
        j++;
    }

    for (; i < this->_size; i++){
        v._data[i] = this->_data[i] + b._data[i];
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = this->_data[i] + b._data[i];
    }
#endif

    return v;
}

Vector Vector::add(const float& b){
    Vector v{this->_size};

 #ifdef __ARM_NEON
    float32x4_t constant;
    const float init[4] = {b, b, b, b};
    constant = vld1q_f32(init);
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4<= this->_size; i += 4){
        vst1q_f32(v._data + i,  vaddq_f32( constant, vld1q_f32(this->_data + i)));
    }
    for (; i < this->_size; i++){
        v._data[i] = this->_data[i] + b;
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = this->_data[i] + b;
    }
#endif


    return v;
}

void Vector::addi(const Vector& b){

 #ifdef __ARM_NEON
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4<= this->_size; i += 4){
        vst1q_f32(this->_data + i,  vaddq_f32( vld1q_f32(this->_data + i), vld1q_f32(b._data + i)));
    }
    for (; i < this->_size; i++){
        this->_data[i] = this->_data[i] + b._data[i];
    }
#else
    for (int i = 0; i < this->_size; i++){
        this->_data[i] = this->_data[i] + b._data[i];
    }
#endif
}


Vector Vector::sub(const Vector& b){
    Vector v{this->_size};


 #ifdef __ARM_NEON
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(v._data + i,  vsubq_f32( vld1q_f32(this->_data + i), vld1q_f32(b._data + i)));
    }
    for (; i < this->_size; i++){
        v._data[i] = this->_data[i] - b._data[i];
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = this->_data[i] - b._data[i];
    }
#endif

    return v;
}

Vector Vector::sub(const float& b){
    Vector v{this->_size};

 #ifdef __ARM_NEON
    float32x4_t constant;
    const float init[4] = {b, b, b, b};
    constant = vld1q_f32(init);
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(v._data + i,  vsubq_f32( vld1q_f32(this->_data + i), constant));
    }
    for (; i < this->_size; i++){
        v._data[i] = this->_data[i] - b;
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = this->_data[i] - b;
    }
#endif

    return v;
}


Vector Vector::mul(const Vector& b){
    Vector v{this->_size};

 #ifdef __ARM_NEON
    int i = 0;
    if (this->_size>=4)
    for (; i + 4 <=  this->_size; i += 4){
        vst1q_f32(v._data + i,  vmulq_f32( vld1q_f32(this->_data + i), vld1q_f32(b._data + i)));
    }
    for (; i  < this->_size; i++){
        v._data[i] = this->_data[i] * b._data[i];
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = this->_data[i] * b._data[i];
    }
#endif
    return v;
}


Vector Vector::mul(const float& b){
    Vector v{this->_size};

 #ifdef __ARM_NEON
    
    float32x4_t constant;
    const float init[4] = {b, b, b, b};
    constant = vld1q_f32(init);
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= v._size; i += 4){
        vst1q_f32(v._data + i,  vmulq_f32( vld1q_f32(this->_data + i), constant));
    }
    for (; i < v._size; i++){
        v._data[i] = this->_data[i] * b;
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = this->_data[i] * b;
    }
#endif
    return v;
}


Vector Vector::div(const Vector& b){
    Vector v{this->_size};
 #ifdef __ARM_NEON
    float32x4_t rep;
    float32x4_t reciprocal;
    int i = 0;
    if (this->_size>=4)
    for (; i +4  <= this->_size; i += 4){
        rep = vld1q_f32(b._data+i);
        reciprocal = vrecpeq_f32(rep);
        reciprocal = vmulq_f32(vrecpsq_f32(rep, reciprocal), reciprocal);
        reciprocal = vmulq_f32(vrecpsq_f32(rep, reciprocal), reciprocal);
        vst1q_f32(v._data + i,  vmulq_f32( vld1q_f32(this->_data + i), reciprocal));
    }
    for (; i < this->_size; i++){
        v._data[i] = this->_data[i] * b._data[i];
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = this->_data[i] * b._data[i];
    }
#endif

    return v;
}


Vector Vector::div(const float& b){
    Vector v{this->_size};
 #ifdef __ARM_NEON
    
    float32x4_t constant;
    const float init[4] = {b, b, b, b};
    constant = vld1q_f32(init);
    int i = 0;
    float32x4_t reciprocal = vrecpeq_f32(constant);
    reciprocal = vmulq_f32(vrecpsq_f32(constant, reciprocal), reciprocal);
    reciprocal = vmulq_f32(vrecpsq_f32(constant, reciprocal), reciprocal);
    if (this->_size>=4)
    for (; i  + 4 <= v._size; i += 4){
        vst1q_f32(v._data + i,  vmulq_f32( vld1q_f32(this->_data + i), reciprocal));
    }
    for (; i < v._size; i++){
        v._data[i] = this->_data[i] / b;
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = this->_data[i] / b;
    }
#endif
    return v;
}


Vector Vector::inv(){
    Vector v{this->_size};
 #ifdef __ARM_NEON
    float32x4_t rep;
    float32x4_t reciprocal;
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 < this->_size; i += 4){
        rep = vld1q_f32(this->_data + i);
        reciprocal = vrecpeq_f32(rep);
        reciprocal = vmulq_f32(vrecpsq_f32(rep, reciprocal), reciprocal);
        reciprocal = vmulq_f32(vrecpsq_f32(rep, reciprocal), reciprocal);
        vst1q_f32(v._data + i,  reciprocal);
    }
    for (; i < this->_size; i++){
        v._data[i] =1/this->_data[i] ;
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = 1/this->_data[i];
    }
#endif

    return v;
}

ostream& Vector::print(ostream& os, Vector const & v){
    int i;
    os << "[";
    for (i = 0; i < v._size-1; i++){
        os << v._data[i] << ", ";
    } 
    os << v._data[i] << "]";
    return os;
}


float dotHelper( float *a, float *b, unsigned int len){
    unsigned int i;
    float sum = 0.0;
    unsigned int len4 = len - len%4;
#ifdef __ARM_NEON
    float32x4_t acc = vdupq_n_f32(0.0f);
#pragma omp parallel private(i,acc) num_threads(4) 
    {
        #pragma omp for reduction(+:sum)
          for (i=0; i<len4; i+=4) 
            { 
                acc = vld1q_f32(a + i) * vld1q_f32(b + i);
                sum +=  acc[0] + acc[1] + acc[2] + acc[3];
            }
    }
    for (i = len4; i<len; i++) sum +=  a[i] * b[i];

    #endif
    return sum;

}



float Vector::dot(const Vector& b){
float r = 0.0f;

 #ifdef __ARM_NEON

    r = dotHelper(this->_data,b._data,this->_size);
    
#else
    for (int i = 0; i < this->_size; i++){
        r += this->_data[i] * b._data[i];
    }
#endif

    
    return r;
}

Vector::Vector(const Vector &v){
    this->_size = v._size;
    this->_data = new float[v._size];
 #ifdef __ARM_NEON
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(this->_data + i, vld1q_f32(v._data + i) );
    }
    for (; i < this->_size; i++){
        this->_data[i] = v._data[i];
    }

#else
    for (int i = 0; i < this->_size; i++){
        this->_data[i] = v._data[i];
    }
#endif
}

void Vector::clone(const Vector& v){
    delete [] this->_data;
    this->_size = v._size;
    this->_data = new float[v._size];
 #ifdef __ARM_NEON
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(this->_data + i, vld1q_f32(v._data + i) );
    }
    for (; i < this->_size; i++){
        this->_data[i] = v._data[i];
    }

#else
    for (int i = 0; i < this->_size; i++){
        this->_data[i] = v._data[i];
    }
#endif
}




Vector Vector::exponential(){
    Vector v{this->_size};
    for (int i = 0; i < this->_size; i++){
        v._data[i] = exp(this->_data[i]);
    }
    return v;
}



Vector Vector::absolute(){
    Vector v{this->_size};

 #ifdef __ARM_NEON
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(v._data + i, vabsq_f32(vld1q_f32(this->_data + i)));
    }
    for (; i < this->_size; i++){
        v._data[i] = abs(this->_data[i]);
    }

#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = abs(this->_data[i]);
    }
#endif
    return v;
}



Vector Vector::append(const float& e){
    Vector v{this->_size+1};
    int i = 0;
 #ifdef __ARM_NEON
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(v._data + i, vld1q_f32(v._data + i));
    }
    for (; i < this->_size; i++){
        v._data[i] = this->_data[i];
    }
#else
    for (; i < this->_size; i++){
        v._data[i] = this->_data[i];
    }
#endif
    v._data[i] = e;
    return v;
}

void Vector::appendInPlace(const float& e){
    unsigned int size = this->_size + 1;
    float *data = new float[size];
    int i = 0;

 #ifdef __ARM_NEON
    if (this->_size>=4)
    for (; i + 4 <= this->_size; i += 4){
        vst1q_f32(data + i, vld1q_f32(this->_data + i));
    }
    for (; i < this->_size; i++){
        data[i] = this->_data[i];
    }
#else
    for (; i < this->_size; i++){
        data[i] = this->_data[i];
    }
#endif
    data[i] = e;
    delete [] this->_data;
    this->_data = data;
    this->_size = size;
}


Vector Vector::relu(){
    Vector v{this->_size};

 #ifdef __ARM_NEON
    float32x4_t constant;
    const float init[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    constant = vld1q_f32(init);

    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(v._data + i, vmaxq_f32(constant,vld1q_f32(v._data + i)) );
    }
    for (; i < this->_size; i++){
        v._data[i] = max(this->_data[i],0.0f);
    }
#else
    for (int i = 0; i < this->_size; i++){
        v._data[i] = max(this->_data[i],0.0f);
    }
#endif
    return v;
}


Vector Vector::relu_derivated(){
    Vector v{this->_size};
    float m;
    for (int i = 0; i < this->_size; i++){
        m = max(this->_data[i],0.0f);
        v._data[i] = m/(m==0.0f?1.0f:m);
    }
    return v;
}


float Vector::maximum(){
    float m = this->_data[0];


 #ifdef __ARM_NEON
    int i = 0;
    float32x4_t maxi;
    float *init = new float[4];
    maxi = vld1q_f32(init);
    i = 0;
    if (this->_size>=4){
        for (i = 0; i < 4; i++)init[i] = this->_data[i];
        for (i = 0; i  + 4 <= this->_size; i += 4) maxi = vmaxq_f32(maxi,vld1q_f32(this->_data + i));
    }
    for (; i < this->_size; i++) m = max(this->_data[i], m);
    if (this->_size>=4){
        vst1q_f32(init, maxi);
        for (i = 0; i < 4; i++)m = max(m, init[i]);
    }
    
    delete [] init;

#else
    for (int i = 0; i < this->_size; i++){
        m = max(this->_data[i],m);   
    }
#endif

    return m;
}

Matrix Vector::extern_product(const Vector& b){
    Matrix m{b._size,this->_size};


 #ifdef __ARM_NEON
    int i = 0;
    float32x4_t vectorA;
    int j = 0;
    unsigned int cols4 = (m._cols<4)?  0 : m._cols - m._cols%4;
    if (b._size>=4)
        for(int i = 0; i + 4 <= b._size;i+=4){
            vectorA = vld1q_f32(b._data + i);
            j = 0;
#pragma omp parallel private(j) num_threads(4) 
            {
                #pragma omp for
                for(j = 0; j  < cols4;j+=4){
                    vst1q_f32(m._data[i] + j, vmulq_f32(vectorA,vld1q_f32(this->_data + j)) );
                }
            }
            for(; j < this->_size;j++){
                m.set(i,j,this->_data[j]*b._data[i]);
            }
        }
    for(; i < b._size;i++){
        for(j = 0; j < this->_size;j++)
            m.set(i,j,this->_data[j]*b._data[i]);
    }
    

#else
    for(int i = 0; i < b._size;i++){
        for(int j = 0; j < this->_size;j++){
            m.set(i,j,this->_data[j]*b._data[i]);
        }
    }
#endif

    return m;
}


Matrix Vector::toMatrix(){
    Matrix m{1, this->_size};


 #ifdef __ARM_NEON
    int i = 0;
    if (this->_size>=4)
    for (; i  + 4 <= this->_size; i += 4){
        vst1q_f32(m._data+i, vld1q_f32(this->_data + i) );
    }
    for (; i < this->_size; i++){
        m._data[0][i] = this->_data[i];
    }
#else
    for (int i = 0; i < this->_size; i++){
        m._data[0][i] = this->_data[i];
    }
#endif

    return m;
}
