CC = nvcc

.PHONY: all run

all: stable.cu
	$(CC) stable.cu -o stable

run: stable
	./stable
