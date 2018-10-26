#!gmake

TARGET_ARCH =

CC = gcc
LD = gcc
#LDFLAGS = -lm -g -pg -lgsl -lgslcblas
LDFLAGS = -lm 
LIBS = 
INC =

# an extra LDFLAG that is useful for segmentation faults and memory testing is -lefence

CFLAGS = $(LIBS) $(INC) -g -Wall
#CFLAGS = -g -Wall -O4 $(LIBS) $(INC)

# use this command to erase files.
RM = /bin/rm -f

# list of generated object files.
OBJS = wallpaper.o groups.o shapes.o fluke.o

# program executable file name.
PROG = wallpaper.exe

# top-level rule, to compile everything.
all: $(PROG) makefile

# linking rule remains the same as before.
$(PROG): $(OBJS)
	$(LD) -o $(PROG) $(OBJS) $(LDFLAGS)  

# now comes a meta-rule for compiling any "C" source file.
%.o: %.c
	$(CC) $(CFLAGS) -c $<

# rule for cleaning re-compilable files.
clean:
	$(RM) $(PROG) $(OBJS)
