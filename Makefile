# calls:
CC         = g++
CFLAGS     = -w -c -Wall -O3 -std=c++11
LDFLAGS    = 
EXECUTABLE = meganan

SOURCES    = previz.cpp skeleton.cpp motion.cpp displaySkeleton.cpp parse_stl.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o previz
