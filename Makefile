#######################################################################################################

# Mac OS X
#INCLUDE_PATH      = -I/usr/local/include/
#LIBRARY_PATH      = -L/usr/local/lib/
#OPENGL_LIBS       = -framework OpenGL -framework GLUT

# # Linux
INCLUDE_PATH      = -I/usr/local/include/
LIBRARY_PATH      = -L/usr/local/lib/
OPENGL_LIBS       = -lglut -lGL -lX11

# # Windows / Cygwin
# INCLUDE_PATH      = -I/usr/include/opengl
# LIBRARY_PATH      = -L/usr/lib/w32api
# OPENGL_LIBS       = -lglut32 -lopengl32

#######################################################################################################

CFLAGS = $(INCLUDE_PATH) -I./include
LFLAGS = $(LIBRARY_PATH)
LIBS = $(OPENGL_LIBS)

HPC_SDK=/opt/nvidia/hpc_sdk/Linux_x86_64/2020/compilers/bin
CC=$(HPC_SDK)/pgc++
#CC=g++ -O3
#FLAGS=-Minfo
#FLAGS=-fast -Minfo
FLAGS=-ta=tesla -Minfo
#FLAGS=-ta=multicore -fast -Minfo
#FLAGS=-Wall

all: fem

fem: main.o
	$(CC) $(FLAGS) $(LFLAGS) main.o $(LIBS) -o fem

main.o: main.cpp matrix_function.h
	$(CC) $(FLAGS) $(CFLAGS) -c main.cpp -o main.o

clean:
	rm -f fem main.o
