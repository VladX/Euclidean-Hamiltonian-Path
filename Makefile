all:
	g++ -Wall -O3 -std=c++11 -march=native -mtune=native main.cpp -o christofides -lGL -lSDL2
clean:
	rm -f christofides