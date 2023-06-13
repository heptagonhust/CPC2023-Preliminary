CC=sw9gcc
CXX=sw9g++

CFLAGS= -mslave -msimd -mieee 
CXXFLAGS= -mhost -mieee -mftz -fpermissive 
INCLUDE=-I.
INCLUDE+=-I./spmv

EXE=pcg_solve

all: $(EXE)

$(EXE): main.o pcg.o spmv_master.o slave.o spmv_slave.o
	$(CXX) -mhybrid -o $(EXE) $^ -L. -lpcg_solve

main.o:	main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ 

pcg.o: pcg.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@	

spmv_master.o: spmv/spmv_master.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@	

slave.o: slave.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

spmv_slave.o: spmv/spmv_slave.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

clean:
	$(RM) $(EXE) *.o;
