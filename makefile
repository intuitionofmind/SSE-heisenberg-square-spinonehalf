CXX = g++
CXXFLAGS = $(INC) -O2 -Wall -ansi -funroll-loops -std=gnu++0x 
OBJS = main.o auxi.o update.o

TARGET = Hesenberg 

all:$(TARGET)
$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)

main.o:main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)
auxi.o:auxi.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)
update.o:update.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIB)

clean:
	rm -f $(TARGET) *.o log *.dat
