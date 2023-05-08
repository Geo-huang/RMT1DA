
CXX              := g++
CXXFLAGS         := -O3  -w -g -Wfatal-errors

EIGEN_INCLUDE    := -I../../contrib/eigen


INCLUDE:= -I../../src/
INCLUDES:= $(wildcard ../../src/*.h)
SRCS:= $(wildcard ../../src/*.cpp)
OBJS:= $(patsubst %.cpp, %.o, $(SRCS))


# How to compile C++
%.o : %.cpp 
	@echo "RMT1DA is compiling C++ "$<"..."
	@$(CXX) $(CXXFLAGS) $(INCLUDE)  $(EIGEN_INCLUDE) -c $< -o $@


#--------------------------------------------------------------
# How to link C++ and Petsc
all: $(OBJS)

RMT1DA: $(OBJS) RMT1DA.o 
	 @$(CXX) $(CXXFLAGS) -o RMT1DA RMT1DA.o $(OBJS)
	 @chmod a+x RMT1DA

clean: 
	@rm -rf $(OBJS) RMT1DA RMT1DA.o  

# Check the compile and link options.	
echo:	
	@echo "C++ Compiler:\n$(CXX)\n"
	@echo "CXXFLAGS:    \n$(CXXFLAGS)\n"
	@echo "Header Files:\n$(INCLUDES)\n"
	@echo "Source Files:\n$(SRCS)\n"
	@echo "Object Files:\n$(OBJS)\n"
	@echo "INCLUDES:    \n$(INCLUDE_DIR)\n"
	
	
