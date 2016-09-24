SRC = test.cpp \
      writePPM.cpp \
      parseScene.cpp \
	  tinyXML/tinystr.cpp \
	  tinyXML/tinyxml.cpp \
	  tinyXML/tinyxmlerror.cpp \
	  tinyXML/tinyxmlparser.cpp 

OBJ = $(SRC:.cpp=.o)
CFLAGS = -O2
LIBS = 

.cpp.o:
	g++ -c $< $(CFLAGS)

all: test

test: $(OBJ) 
	g++ $(OBJ) -o raytracer $(CFLAGS) $(LIBS)

clean:
	rm -f *.o raytracer test.ppm
