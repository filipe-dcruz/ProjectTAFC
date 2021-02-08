
CC = g++
CFLAGS = -g -Wall
LDFLAGS = -I.
OBJFILES = species.o input.o initial.o pic.o
TARGET = pic1d

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
		rm -f $(OBJFILES) $(TARGET) *~
