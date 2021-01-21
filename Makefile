CC = gcc
CFLAGS = -I.
LDFLAGS =
OBJFILES = initial.o pic.o main.o
TARGET = pic1d

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
		rm -f $(OBJFILES) $(TARGET) *~
