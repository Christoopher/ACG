HEADERS += \
    ../../../src/RigidBody.h \
    ../../../src/Quaternion.h \
    ../../../src/Physics.h \
    ../../../src/OpenGLViewer.h \
    ../../../src/ObjLoader.h

SOURCES += \
    ../../../src/main.cpp

LIBS += -lGL -lGLU -lglut -lglfw -lGLEW -llapack -lblas -larmadillo
	
