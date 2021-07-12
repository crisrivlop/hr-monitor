#include "pyramid.h"
#include <exception>
#include <cmath>
#include "../utils/size.h"



Pyramid::Pyramid(){
    this->_length = 0;
    this->_pyramid = 0;

}


Pyramid::Pyramid(int l){
    this->_length = l;
    this->_pyramid = new Matrix[_length]();

}


Pyramid::Pyramid(const Pyramid& p){
    _pyramid = 0;
    _length = 0;
    this->clone(p);
}


void Pyramid::clone(const Pyramid& p){
    delete [] _pyramid;
    this->_length = p._length;
    this->_pyramid = new Matrix[_length]();
    for (int i = 0; i < this->_length;i++)
        this->_pyramid[i] = p._pyramid[i];
}


Pyramid Pyramid::add(const Pyramid& p){
    Pyramid p2(p._length);
    for (int i = 0; i < this->_length;i++)
        p2._pyramid[i] = this->_pyramid[i] + p._pyramid[i];
    return p2;
}

Pyramid Pyramid::sub(const Pyramid& p){
    Pyramid p2(p._length);
    for (int i = 0; i < this->_length;i++)
        p2._pyramid[i] = this->_pyramid[i] - p._pyramid[i];
    return p2;
}

Pyramid Pyramid::mul(float b){
    Pyramid p2(this->_length);
    for (int i = 0; i < this->_length;i++)
        p2._pyramid[i] = this->_pyramid[i]*b;
    return p2;
}


void Pyramid::IgnoreHighestAndLowestFrequencies(){
    if (_length > 0){
        _pyramid[0].fill(0);
        _pyramid[_length-1].fill(0);
    }
}



void Pyramid::BuildLpyr(Matrix im, string filter1, string filter2, string edges){

    if (!im.rows() || !im.cols()){
        //matriz invalida
        string errorMessage = "Invalid matrix with size (0,0)";
        cout << errorMessage << endl;
        throw exception();
    }

    MatrixSize im_sz(im.cols(),im.rows());

    filter1 = filter1.empty()?"binom5":filter1;
    Matrix f1 = namedFilter(filter1).toMatrix().transpose();
    
    if(f1.rows()>1 && f1.cols()>1){
        string errorMessage = "Filter1 should be a 1D filter (i.e., a vector)";
        cout << errorMessage << endl;
        throw exception();
    }

    filter2 = filter2.empty()?"binom5":filter2;
    Matrix f2 = namedFilter(filter2).toMatrix().transpose();

    if(f2.rows()>1 && f2.cols()>1){
        string errorMessage = "Filter2 should be a 1D filter (i.e., a vector)";
        cout << errorMessage << endl;
        throw exception();
    }



    int max_height =  1 + maxPyramidHeight(im_sz.width(),im_sz.height(), max(f1.rows(),f2.rows()));
    this->_length = max_height;
    delete [] _pyramid;
    this->_pyramid = new Matrix[_length]();
    BuildLpyrAux(im,max_height,f1,f2,edges);
    
}



void Pyramid::BuildLpyrAux(Matrix im, int height, Matrix f1, Matrix f2, string edges=""){

    MatrixSize im_sz(im.cols(),im.rows());
    int max_height =  1 + maxPyramidHeight(im_sz.width(),im_sz.height(), max(f1.cols(),f2.cols()));
    int ht = height;
    
    if (height <= 0){
        ht = max_height;
    }
    else if( height > max_height) {
        string errorMessage = "Cannot build pyramid higher than ";
        cout << errorMessage << max_height <<" levels." << endl;
        throw exception();
    }
    edges = edges.empty()?"reflect1":edges;

    Matrix lo,lo2;
    Matrix hi,hi2;
    int index21[2] = {2,1};
    int index12[2] = {1,2};
    int index11[2] = {0,0};

    Matrix pyr;
    MatrixSize pind(0,0);
    MatrixSize int_sz(0,0);
    
    if (ht <=1){
        this->_pyramid[this->_length -1] = im;
    }
    else{
        
        if(im.cols() == 1){
            lo2 = corrDn(im,f1,edges,index21,index11);
        }
        else if(im.rows() == 1){
            lo2 = corrDn(im, f1.transpose(), edges, index12, index11);
        }
        else{
            lo = corrDn(im, f1.transpose(), edges, index12, index11);
            int_sz = MatrixSize(lo.cols(),lo.rows());
            lo2 = corrDn(lo, f1, edges, index21, index11);
            
        }
        
        this->BuildLpyrAux(lo2, ht-1, f1, f2, edges);
        if (im_sz.height() == 1)
             hi2 = upConv(lo2, f2.transpose(), edges, im_sz, index12, index11);
        else if (im_sz.width() == 1)
             hi2 = upConv(lo2, f2, edges, im_sz, index21, index11);
        else{
             hi = upConv(lo2, f2, edges, int_sz, index21, index11);
             hi2 = upConv(hi, f2.transpose(), edges, im_sz, index12, index11);

        }
        hi2 = im - hi2;
        
        this->_pyramid[this->_length - ht] = hi2;
    }

    
}


Vector Pyramid::namedFilter(string filterName){
    float sqrt_of_2 = sqrt(2);
    
    if (filterName.compare("binom5") == 0){
        //kernel = sqrt(2) * binomialFilter(str2num(name(6:size(name,2))));
        float kernel_values[] = {0.088388,0.353553,0.530330,0.353553,0.088388};
        return Vector(5, kernel_values);
    }
    else if(filterName.compare("qmf5") == 0){
        float kernel_values[] = {-0.076103, 0.3535534, 0.8593118, 0.3535534, -0.076103};
        return Vector(5,kernel_values);
    }
    else if (filterName.compare("qmf9") == 0) {
        float kernel_values[] = {0.02807382, -0.060944743, -0.073386624, 0.41472545, 0.7973934,
            0.41472545, -0.073386624, -0.060944743, 0.0280738};
        return Vector(9,kernel_values);
    }
    else if (filterName.compare("qmf13") == 0){
        float kernel_values[] =  {-0.014556438, 0.021651438, 0.039045125, -0.09800052,
	      -0.057827797, 0.42995453, 0.7737113, 0.42995453, -0.057827797,
	      -0.09800052, 0.039045125, 0.021651438, -0.014556438};
        return Vector(13,kernel_values);
    }

    else if (filterName.compare("qmf8") == 0){
        float kernel_values[] = {0.00938715, -0.07065183, 0.06942827, 0.4899808,
        0.4899808, 0.06942827, -0.07065183, 0.00938715};
        
        for(int i = 0; i < 8;i++) kernel_values[i] *=  sqrt_of_2;
        return Vector(8,kernel_values);
    }
    else if (filterName.compare("qmf12") == 0){
        float kernel_values[] = {-0.003809699, 0.01885659, -0.002710326, -0.08469594,
	    0.08846992, 0.4843894, 0.4843894, 0.08846992, -0.08469594, -0.002710326,
	    0.01885659, -0.003809699 };
        for(int i = 0; i < 12;i++) kernel_values[i] *=  sqrt_of_2;
        return Vector(12,kernel_values);
    }
    else if (filterName.compare("qmf16") == 0){
        float kernel_values[] = {0.001050167, -0.005054526, -0.002589756, 0.0276414, -0.009666376,
	        -0.09039223, 0.09779817, 0.4810284, 0.4810284, 0.09779817, -0.09039223, -0.009666376,
	        0.0276414, -0.002589756, -0.005054526, 0.001050167 };
        for(int i = 0; i < 16;i++) kernel_values[i] *=  sqrt_of_2;
        return Vector(16,kernel_values);
    }
    else if (filterName.compare("haar") == 0){
        float kernel_values[] = {1.0, 1.0};
        for(int i = 0; i < 2;i++) kernel_values[i] /=  sqrt_of_2;
        return Vector(2,kernel_values);
    }
    else if (filterName.compare("daub2") == 0){
        float kernel_values[] = {0.482962913145, 0.836516303738, 0.224143868042, -0.129409522551};
        return Vector(4,kernel_values);
    }
    else if (filterName.compare("daub3") == 0){
        
        float kernel_values[] = {0.332670552950, 0.806891509311, 0.459877502118, -0.135011020010,
	                -0.085441273882, 0.035226291882};
        return Vector(6,kernel_values);
    }
    else if (filterName.compare("daub4") == 0){
        float kernel_values[] = {0.230377813309, 0.714846570553, 0.630880767930, -0.027983769417,
	            -0.187034811719, 0.030841381836, 0.032883011667, -0.010597401785};
        return Vector(8,kernel_values);
    }
    else if (filterName.compare("gauss5") == 0){
        float kernel_values[] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
        for(int i = 0; i < 5;i++) kernel_values[i] *=  sqrt_of_2;
        return Vector(5,kernel_values);
    }
    else if (filterName.compare("gauss3") == 0){
        float kernel_values[] = {0.25, 0.5, 0.25};
        for(int i = 0; i < 3;i++) kernel_values[i] *=  sqrt_of_2;
        return Vector(5,kernel_values);
    }
    else {
        throw new exception();
    }
}


int Pyramid::maxPyramidHeight(int rows, int cols, int filterSize){

    if (rows < filterSize || cols < filterSize){
        return 0;
    }
    else{
        int rows_floor = (int)(floor(float(rows)/2.0f));
        int cols_floor = (int)(floor(float(cols)/2.0f));
        return 1 + maxPyramidHeight(rows_floor,  cols_floor, filterSize);
    }
}


Matrix Pyramid::rconv2(Matrix a, Matrix b, int center){

    MatrixSize as(a.cols(),a.rows());
    MatrixSize bs(b.cols(),b.rows());
    Matrix small,large;

    if (as.height() <= bs.height() && as.width() <= bs.width()){
        small = a;
        large = b;
    }
    else if (as.height() >= bs.height() && as.width() >= bs.width()){
        small = b;
        large = a;
    }
    else {
        cout << as.height() << "x" << as.width() << " and " << bs.height() << "x" << bs.width() << endl;
        cout << "one arg must be larger than the other in both dimensions!" << endl;
        throw exception();
    }

    int ly = large.rows();
    int lx = large.cols();
    int sy = small.rows();
    int sx = small.cols();

    int sy2 = (sy+center-1)/2;
    int sx2 = (sx+center-1)/2;

    Matrix clarge(large.rows()+sy2*2,large.cols()+sx2*2);
    for (int i = 1;i<=sy2;i++){
        for(int j= 1;j<=sx2;j++){
            //esquina superior izquierda
            float sle=large.get(i,j);
            clarge.set(sy2-i,sx2-j,sle);
            //esquina superior derecha
            float sre=large.get(i,large.cols()-j-1);
            clarge.set(sy2-i,large.cols()+1+j,sre);
            // //esquina inferior izquierda
            float ile=large.get(large.rows()-i-1,j);
            clarge.set(large.rows()+i+1,sx2-j,ile);
            // //esquina inferior derecha
            float ire=large.get(large.rows()-i-1,large.cols()-j-1);
            clarge.set(large.rows()+i+1,large.cols()+1+j,ire);

        }

    }

    for (int i = 0;i < large.rows();i++){
        for(int j= 0;j<sx2;j++){
            //lado izquierdo
            clarge.set(i+sy2,sx2-j-1,large.get(i,j+1));
            //derecha
            clarge.set(i+sy2,sx2+large.cols()+j,large.get(i,large.cols()-j-2));
            
        }

    }

    for (int i = 0;i<sy2;i++){
        for(int j= 0;j<large.cols();j++){
            //arriba
            float up=large.get(i+1,j);
            clarge.set(sy2-i-1,sx2+j,up);
            //abajo
            float bottom=large.get(large.rows()-sy2-i,j);
            clarge.set(large.rows()+(sy2+i),sx2+j,bottom);
        }
    }

    for (int i = 0;i<large.rows();i++){
        for(int j= 0;j<large.cols();j++){
            //center
            float cell=large.get(i,j);
            clarge.set(sy2+i,sx2+j,cell);
        }
    }

    Matrix m(clarge.conv(small));
    return m;

}


Matrix Pyramid::corrDn(Matrix im, Matrix filter, string edges, int* step, int* start, int* stop){

    edges = edges.empty()?"reflect1":edges;

    int * defaultConfig = new int[2];
    defaultConfig[0] = 0;
    defaultConfig[1] = 0;

    int * defaultStep = new int[2];
    defaultStep[0] = 1;
    defaultStep[1] = 1;

    int * defaultStop = new int[2];
    defaultStop[0] = im.rows();
    defaultStop[1] = im.cols();
    

    step = step?step:defaultStep;
    start = start?start:defaultConfig;
    stop = stop?stop:defaultStop;

    Matrix filterInverted(filter.rows(),filter.cols());

    
    //se invierte el vector para generar la correlación
    for (int i = 0; i < filter.rows(); i++){
        for (int j = 0; j < filter.cols(); j++){
            filterInverted.set(i,j,filter.get(filter.rows()-1-i, filter.cols()-1-j));
        }
    }


    Matrix m(rconv2(im,filterInverted));

    int finalRows = ceil((float)(stop[0]-start[0])/(float)step[0]);
    int finalCols = ceil((float)(stop[1]-start[1])/(float)step[1]); 
    Matrix finalCorrelationalMatrix(finalRows, finalCols);

    for (int i = 0; i < finalCorrelationalMatrix.rows(); i++){
        for (int j = 0; j < finalCorrelationalMatrix.cols(); j++){
            finalCorrelationalMatrix.set(i,j,m.get(i*step[0], j*step[1]));
        }
    }
    

    delete defaultConfig;
    delete defaultStep;
    delete defaultStop;
    return finalCorrelationalMatrix;

}

Matrix Pyramid::upConv(Matrix im, Matrix filter, string edges, MatrixSize s, int* step, int* start, int* stop, Matrix* res){
    edges = edges.empty()?"reflect1":edges;

    int * defaultConfig = new int[2];
    defaultConfig[0] = 0;
    defaultConfig[1] = 0;
    int * defaultStep = new int[2];
    defaultStep[0] = 1;
    defaultStep[1] = 1;

    int * defaultStop = new int[2];
    defaultStop[0] = im.rows();
    defaultStop[1] = im.cols();
    
    step = step?step:defaultStep;
    start = start?start:defaultConfig;
    
    defaultStop[0] = step[0] * ((start[0]-1)/step[0]+(int)im.rows());
    defaultStop[1] = step[1] * ((start[1]-1)/step[1]+(int)im.cols());
    stop = stop?stop:defaultStop;

    if ( (int)std::ceil((float)(s.height()-start[0])/(float)step[0]) != im.rows()){
        std::cout << "Bad y dimensions" << std::endl;
        throw exception();
    }
    if ( (int)std::ceil((float)(s.width() -start[1] )/(float)step[1]) != im.cols()){
        std::cout << "Bad x dimensions" << std::endl;
        throw exception();
    }
    Matrix * resDefault = new Matrix((int)s.height()-(int)start[0], (int)s.width()-(int)start[1]);
    resDefault->fill(0);
    res = res?res:resDefault;

    Matrix tmp(res->rows(), res->cols());
    tmp.fill(0);
    
    for(int i = 0;i < im.rows();i++){
        for(int j = 0;j<im.cols();j++){
            tmp.set(start[0]+step[0]*i,start[1]+j*step[1],im.get(i,j));
        }
    }

    Matrix rcon(rconv2(tmp, filter));
    Matrix result =  rcon + (*res);

    delete defaultConfig;
    delete defaultStep;
    delete defaultStop;
    delete resDefault;


    return result;


}



Matrix Pyramid::reconLpyr(){
    int last=0; 
    string edges; 
    string filter;
    int level = _length;
    filter = filter.empty()?"binom5":filter;
    Matrix f = namedFilter(filter).toMatrix();
    edges = edges.empty()?"reflect1":edges;

    Matrix res;
    int index21[2] = {2,1};
    int index12[2] = {1,2};
    int index11[2] = {0,0};
    MatrixSize s(1,1);
    MatrixSize s2(1,1);

    for(int i=this->_length-1; i >= 0;i--){
        s = MatrixSize(_pyramid[i].cols(),_pyramid[i].rows());
        if(i==this->_length-1){
            res=Matrix(_pyramid[i].rows(),_pyramid[i].cols());
            res.fill(0);
        }
        else if (_pyramid[i].rows() == 1){
            res = upConv(res, f.transpose(), edges,s, index12, index11);
        }
        else if (_pyramid[i].cols() == 1){
            res = upConv(res, f, edges,s, index21, index11);
        }
        else{
            s2 = MatrixSize(_pyramid[i+1].cols(),_pyramid[i].rows());
            Matrix hi = upConv(res, f, edges,s2, index21, index11);
            res = upConv(hi, f.transpose(), edges,s, index12, index11);
        }
        if (res.cols() == _pyramid[i].cols() && res.rows() == _pyramid[i].rows())
            res = res + _pyramid[i];
        
    }
    return res;

}


void Pyramid::reconLpyrAux(int level, int last, Matrix filter, string edges){

}

Pyramid::~Pyramid(){
    delete [] _pyramid;
    this->_pyramid = 0;
    this->_length = 0;
}