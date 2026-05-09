CPPFLAGS = -g -Wall -DTEST=0

objects := $(patsubst %.cpp,%.o, $(wildcard *.cpp)) #string substitution all .cpp files in the directory by .o
headers := $(wildcard *.h)
test_objects := $(filter-out main.o, $(objects) )


all: CPPFLAGS := -g -Wall -DTEST=0
all: $(objects)   # implicit compilation of objects
	g++ -o all $(objects) -lsfml-graphics -lsfml-window -lsfml-system
	
test: CPPFLAGS := -g -Wall -DTEST=1 # target specific vaiable
test: $(test_objects)  # compile without main.cpp, main.cpp provided by doctest.h
	g++ -o test $(test_objects)


$(objects) :  $(headers) # recompile even if one header changes

clean:
	rm -f $(objects)
	rm -f all
	rm -f test
