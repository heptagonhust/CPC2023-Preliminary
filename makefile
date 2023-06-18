CC  = sw9gcc
CXX = sw9g++

OBJECTS = main.o pcg.o vector_master.o spmv_master.o coo_matrix.o matrix_utils.o vector_slave.o spmv_slave.o

CFLAGS   = -mslave -msimd -mieee
CXXFLAGS = -mhost -mieee -mftz -fpermissive

INCLUDE = -I.
INCLUDE += -I./spmv
INCLUDE += -I./vector
INCLUDE += -I./matrix

EXE = pcg_solve

all: $(EXE)

$(EXE) : $(OBJECTS)
	$(CXX) -mhybrid -o $(EXE) $^ -L. -lpcg_solve libswperf.a

main.o:	main.cpp pcg.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

pcg.o: pcg.cpp matrix/coo_matrix.h 
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

spmv_master.o: spmv/spmv_master.cpp spmv/spmv_slave.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

vector_master.o: vector/vector_master.cpp vector/vector_def.h pcg_def.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

coo_matrix.o: matrix/coo_matrix.cpp spmv/spmv_def.h pcg.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

matrix_utils.o: matrix/matrix_utils.cpp pcg_def.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

spmv_slave.o: spmv/spmv_slave.c spmv/spmv_def.h slave_def.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

vector_slave.o: vector/vector_slave.c vector/vector_def.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

clean:
	$(RM) $(EXE) *.o;
