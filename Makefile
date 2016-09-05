NAME = euler1d
PREFIX = build

INCLUDES = 

#Profile using gprof
PROFILE= #-pg

#CC = gcc-6
CFLAGS =
#CXX = g++-6
CPPFLAGS =  -std=c++11
LFLAGS =  #-lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lmkl_core -lpthread 
 
libsrcs =$(wildcard *.cpp)     
libobjst = $(libsrcs:.cpp=.o)
libobjs = $(foreach obj,$(libobjst),$(PREFIX)/$(obj))
clibsrcs =$(wildcard *.c)
clibobjst =$(clibsrcs:.c=.o)
clibobjs = $(foreach obj,$(clibobjst),$(PREFIX)/$(obj))

$(NAME): $(libobjs)
	$(CXX) $(LFLAGS) -o $(PREFIX)/$(NAME) $(libobjs) $(LIBS) $(PROFILE)

$(PREFIX)/%.o: %.cpp                    
	$(CXX) $(CPPFLAGS) $(EIGENDEF) -c -o $@ $<  $(INCLUDES) $(PROFILE)

$(PREFIX)/%.o: %.c
	$(CC)  $(CFLAGS) -c -o $@ $<  $(INCLUDES) $(PROFILE)

.PHONY : clean
clean:
	rm -f $(PREFIX)/$(NAME) $(libobjs) $(clibobjs)
	rm -f ./output/*

