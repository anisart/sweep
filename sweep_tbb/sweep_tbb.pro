TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp
INCLUDEPATH += "C:\tbb42_20140122oss\include"

LIBS += $$C:\tbb42_20140122oss\lib\ia32\vc10
LIBS += -ltbb
