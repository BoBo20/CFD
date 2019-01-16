TARGET	= dugks_2d_cavity
OBJ		= dugks_cavity.o
SRC		= dugks_cavity.c
CC		= gcc
CFLAGS	= -lm

all : $(TARGET)

$(TARGET) : $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(CFLAGS)

$(OBJ) : $(SRC)
	$(CC) -c $(SRC) -o $(OBJ)

run:
	./dugks_2d_cavity

.PHONY : clean
clean:
	rm $(OBJ)
