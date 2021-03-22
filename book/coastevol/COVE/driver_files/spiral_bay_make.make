# make with make -f spiral_bay_make.make

# COMPILER and LINKER MACROs
CC=g++
LD=g++

# COMPILER AND LINKER OPTION FLAGS MACRO
# -g option build tables for debugging
# -c option compile but do not try to link (yet)
# -Wall display all warning messages
# -pg does what?!
# -O3 is an optimisation flag, not good for debugging

CFLAGS= -g -c -O3 -pg  $(INCDIR)
LDFLAGS= -g -O3 -pg

# SOURCE FILES MACROS IN DEPENDENCY ORDER? SHOULDNT MATTER THANKS TO HEADERS
SOURCES = ../coastline.cpp ../waveclimate.cpp ./spiral_bay_driver.cpp

# LIBRARIES MACRO
LIBS   = -lm -lstdc++

# OBJECT FILES SAME NAME AS SOURCES MACRO
OBJECTS=$(SOURCES:.cpp=.o)

# EXECUTABLE MACRO
EXECUTABLE=cove

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
