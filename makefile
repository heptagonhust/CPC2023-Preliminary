CC=sw9gcc
CXX=sw9g++

CFLAGS= -mslave -msimd -mieee
CXXFLAGS= -mhost -mieee -mftz -fpermissive
INCLUDE=-I.
INCLUDE+=-I$(CURDIR)/vector

EXE=pcg_solve

all: $(EXE)

$(EXE): main.o pcg.o vector_opt.o slave.o vector_slave.o
	$(CXX) -mhybrid -o $(EXE) $^ -L. -lpcg_solve

main.o:	main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

pcg.o: pcg.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

vector_opt.o: vector/vector_opt.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

slave.o: slave.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

vector_slave.o: vector/vector_slave.c
	$(CC) $(CFLAGS)  -c $< -o $@

clean:
	$(RM) $(EXE) *.o;
