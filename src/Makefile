compiler := g++
opencvlibs := `pkg-config --cflags --libs opencv4`
openmp := -fopenmp -lpthread -D WITH_OPENMP=ON -mavx2 -DENABLE_AVX2=ON

clean:
	rm -rf build
	rm -rf *.o

all: clean main.o


maths.o:
	$(compiler) -O3  -c maths/iavector.cpp
	$(compiler) -O3  -c maths/matrix.cpp

utils.o: maths.o
	$(compiler) -O3 -c utils/mat2matrix.cpp utils/size.cpp $(opencvlibs) 

pyramid.o: utils.o
	$(compiler) -O3 -c pyramid/pyramid.cpp $(opencvlibs)

main.o: maths.o utils.o pyramid.o
	mkdir build
	$(compiler) -O3 main.cpp -o hr-monitor iavector.o matrix.o mat2matrix.o size.o pyramid.o $(opencvlibs) 
	mv *.o build/
	mv hr-monitor build/

