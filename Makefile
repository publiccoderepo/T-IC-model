INC_DIR = 
LIB_DIR = 
LIBS = 
CFLAGS = -std=c++11 -O3 -fopenmp 
CC = g++
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)

dtinfmax: $(OBJS)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) $? -o $@ $(LIBS)
.cpp.o:
	$(CC) $(CFLAGS) $(INC_DIR) -c $<
clean:
	-rm -f *.o
	-rm -f dtinfmax

