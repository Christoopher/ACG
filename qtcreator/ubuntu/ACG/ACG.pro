HEADERS += \
    ../../../src/RigidBody.h \
    ../../../src/Quaternion.h \
    ../../../src/Physics.h \
    ../../../src/OpenGLViewer.h \
    ../../../src/ObjLoader.h \
    ../../../src/Contact.h \
    ../../../src/CollisionResponse.h \
    ../../../src/CollisionDetection.h \
    ../../../src/UTLog.h

SOURCES += \
    ../../../src/main.cpp

LIBS += -lGL -lGLU -lglut -lglfw -lGLEW -llapack -lblas -larmadillo -fopenmp
	
