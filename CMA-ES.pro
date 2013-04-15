TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    cmaes.cpp

HEADERS += \
    cmaes.h

# eigen3
INCLUDEPATH += /usr/local/include/eigen3

# boost
INCLUDEPATH += /usr/local/opt/boost/include
LIBS += -L/usr/local/opt/boost/lib/ \
  -lboost_chrono-mt \
  -lboost_date_time-mt \
  -lboost_filesystem-mt \
  -lboost_graph-mt \
  -lboost_iostreams-mt \
  -lboost_locale-mt \
  -lboost_math_c99-mt \
  -lboost_math_c99f-mt \
  -lboost_math_c99l-mt \
  -lboost_math_tr1-mt \
  -lboost_math_tr1f-mt \
  -lboost_math_tr1l-mt \
  -lboost_prg_exec_monitor-mt \
  -lboost_program_options-mt \
  -lboost_python-mt \
  -lboost_random-mt \
  -lboost_regex-mt \
  -lboost_serialization-mt \
  -lboost_signals-mt \
  -lboost_system-mt \
  -lboost_thread-mt \
  -lboost_timer-mt \
  -lboost_unit_test_framework-mt \
  -lboost_wave-mt \
  -lboost_wserialization-mt
