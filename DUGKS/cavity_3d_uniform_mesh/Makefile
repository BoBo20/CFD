TARGET	= dugks_3d_cavity
OBJ		= dugks_cavity.o
SRC		= dugks_cavity.c
CC		= gcc
CFLAGS	= -lm

all : $(TARGET)

$(TARGET) : $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(CFLAGS)

$(OBJ) : $(SRC)
	$(CC) -c $(SRC) -o $(OBJ) -O2

run:
	./$(TARGET)

.PHONY : clean
clean:
	rm -f $(TARGET) $(OBJ)
