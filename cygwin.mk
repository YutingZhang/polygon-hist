SHARED_LIB_NAME=zpolyhist
SHARED_LIB=cyg$(SHARED_LIB_NAME).dll
SHARED_LIB_CYG=lib$(SHARED_LIB_NAME).dll.a
EXECUTABLE1=polyint-static-demo.exe
EXECUTABLE2=zpolyhist-test.exe
CFLAGS=$(P_CFLAGS) `pkg-config opencv --cflags` -fPIC -DBUILDING_DLL
EXE_LIBS=`pkg-config opencv --libs` $(P_LIBS)
LIB_LIBS=$(P_LIBS) -fPIC \
    -Wl,--out-implib=$(SHARED_LIB_CYG) \
    -Wl,--export-all-symbols \
    -Wl,--enable-auto-import

include noplatform.mk
