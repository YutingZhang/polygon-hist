SRC_DIR=polyint
OBJS1=polyint_main.o sliding_hist.o polyint.o poly_main.o zpolyhist.o
OBJS_LIB=sliding_hist.o polyint.o zpolyhist.o
OBJS2=zpolyhist_test.o
RM-F=rm -f

.PHONY: all sharelib clean veryclean rebuild

all: veryclean $(EXECUTABLE1) $(SHARED_LIB) $(EXECUTABLE2)

sharelib: veryclean $(SHARED_LIB)

clean:
	@$(RM-F) *.o
	@$(RM-F) *.d

veryclean: clean
	@$(RM-F) $(EXECUTABLE1) $(EXECUTABLE2) $(SHARED_LIB) $(SHARED_LIB_CYG)

rebuild: veryclean all

$(EXECUTABLE1): $(OBJS1)
	$(CC) $(OBJS1) -o $(EXECUTABLE1) $(EXE_LIBS)

$(SHARED_LIB): $(OBJS_LIB)
	$(CC) -shared $(OBJS_LIB) -o $(SHARED_LIB) $(LIB_LIBS) 

$(EXECUTABLE2): $(SHARED_LIB) $(OBJS2)
	$(CC) $(OBJS2) -o $(EXECUTABLE2) $(EXE_LIBS) -L. -lzpolyhist

polyint_main.o: $(SRC_DIR)/polyint_main.cpp $(SRC_DIR)/config.h $(SRC_DIR)/sliding_hist.h $(SRC_DIR)/hashspecial.h $(SRC_DIR)/zutils.h $(SRC_DIR)/polyint.h
	$(CC) $(CFLAGS) -c $< -o $@

sliding_hist.o: $(SRC_DIR)/sliding_hist.cpp $(SRC_DIR)/config.h $(SRC_DIR)/sliding_hist.h $(SRC_DIR)/zutils.h
	$(CC) $(CFLAGS) -c $< -o $@

polyint.o: $(SRC_DIR)/polyint.cpp $(SRC_DIR)/config.h $(SRC_DIR)/zutils.h $(SRC_DIR)/polyint.h
	$(CC) $(CFLAGS) -c $< -o $@

zpolyhist.o: $(SRC_DIR)/zpolyhist.cpp
	$(CC) $(CFLAGS) -c $< -o $@

poly_main.o: $(SRC_DIR)/poly_main.cpp
	$(CC) $(CFLAGS) -c $< -o $@

zpolyhist_test.o: $(SRC_DIR)/zpolyhist_test.cpp
	$(CC) $(CFLAGS) -c $< -o $@

