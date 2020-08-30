TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        custommodel.cpp \
        generator.cpp \
        integrator.cpp \
        linearalgebra.cpp \
        main.cpp \
        mathmodel.cpp

HEADERS += \
    custommodel.h \
    generator.h \
    integrator.h \
    linearalgebra.h \
    mathmodel.h \
    simplealgorithms.h
