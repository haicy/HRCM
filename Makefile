CC = g++
CFLAG = -O3 -Wall -std=c++0x
objects = compress.cpp decompress.cpp

lngc: $(objects)
	$(CC) -o hrcm $(objects) $(CFLAG)

clean:
	rm *.o
