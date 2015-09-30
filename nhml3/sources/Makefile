
# Make all executables by default:
# $ make
# 
# These targets each make only a single executable:
# $ make eval_nh
# $ make star_nh
# $ make shake_nh
#
# Delete all exectuables and object files:
# $ make clean



DEP_OBJECTS = simplex.o maseread.o compute_ml.o eigen.o
EXECUTABLES = eval_nh star_nh shake_nh

CFLAGS  = -c 

# ANSI c compiler and linker
CC = cc -O

# Unless a specific target is specified, make all executables.
all: $(EXECUTABLES)

# I think this says how to make .o files from .c files...
%.o : %.c
	$(CC) $(CFLAGS) $?


eval_nh: $(DEP_OBJECTS) eval_nh.o
	$(CC) -o $@ $(DEP_OBJECTS) eval_nh.o -lm

star_nh: $(DEP_OBJECTS) star_nh.o
	$(CC) -o $@ $(DEP_OBJECTS) star_nh.o -lm

shake_nh: $(DEP_OBJECTS) shake_nh.o
	$(CC) -o $@ $(DEP_OBJECTS) shake_nh.o -lm


# The force flag prevents complaints when the directory is already clean.
clean:
	rm --force *.o $(EXECUTABLES)
