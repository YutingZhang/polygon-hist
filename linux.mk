SHARED_LIB=libzpolyhist.so
EXECUTABLE1=polyint-static-demo
EXECUTABLE2=zpolyhist-test
CFLAGS=$(P_CFLAGS) `pkg-config opencv --cflags` -fPIC -DBUILDING_DLL -fvisibility=hidden
EXE_LIBS=`pkg-config opencv --libs` $(P_LIBS)
LIB_LIBS=$(P_LIBS) -fPIC -fvisibility=hidden 

include noplatform.mk
