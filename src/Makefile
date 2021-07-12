compiler := g++
neon := -DHAVE_NEON=1
opencvlibs := `pkg-config --cflags --libs opencv`
openmp := -fopenmp -lpthread -D WITH_OPENMP=ON

clean:
	rm -rf build
	rm -rf *.o

all: clean main.o


maths.o:
	$(compiler) -O3  -c maths/iavector.cpp ${neon}
	$(compiler) -O3  -c maths/matrix.cpp ${neon}

utils.o: maths.o
	$(compiler) -O3 -c utils/mat2matrix.cpp utils/size.cpp $(opencvlibs)  ${neon}

pyramid.o: utils.o
	$(compiler) -O3 -c pyramid/pyramid.cpp $(opencvlibs) ${neon}

main.o: maths.o utils.o pyramid.o
	mkdir build
	$(compiler) -O3 main.cpp -o hr-monitor iavector.o matrix.o mat2matrix.o size.o pyramid.o $(opencvlibs)  ${neon}
	mv *.o build/
	mv hr-monitor build/

