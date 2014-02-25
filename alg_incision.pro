TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    alg_incision.c \
    main.c \
    table_trig.c

HEADERS += \
    alg_incision.h

# Specifies the head file directories that should be searched.
INCLUDEPATH += G:/opencv_mingw_hpbuild/install/include/
INCLUDEPATH += G:/opencv_mingw_hpbuild/install/include/opencv/
INCLUDEPATH += G:/opencv_mingw_hpbuild/install/include/opencv2/

# Contains the list of library files that should be linked with the project.
LIBS += G:/opencv_mingw_hpbuild/install/lib/libopencv_core244d.dll.a
LIBS += G:/opencv_mingw_hpbuild/install/lib/libopencv_highgui244d.dll.a
