nozzle: av.o boundary.o data.o derivative.o eh.o initialise.o maccormack.o main.o memory.o roe.o schemes.o timestep.o
	gcc -Wall -lm -o nozzle av.o boundary.o data.o derivative.o eh.o initialise.o maccormack.o main.o memory.o roe.o schemes.o timestep.o

av.o: av.c main.h av.h
	gcc -Wall -c av.c

boundary.o: boundary.c main.h boundary.h
	gcc -Wall -c boundary.c

data.o: data.c main.h data.h
	gcc -Wall -c data.c

derivative.o: derivative.c main.h derivative.h
	gcc -Wall -c derivative.c

eh.o: eh.c main.h derivative.h eh.h
	gcc -Wall -c eh.c

initialise.o: initialise.c main.h initialise.h
	gcc -Wall -c initialise.c

maccormack.o: maccormack.c main.h av.h derivative.h maccormack.h
	gcc -Wall -c maccormack.c

main.o: main.c boundary.h data.h eh.h initialise.h maccormack.h memory.h roe.h timestep.h
	gcc -Wall -c main.c

memory.o: memory.c main.h memory.h
	gcc -Wall -c memory.c

roe.o: roe.c main.h roe.h schemes.h
	gcc -Wall -c roe.c

schemes.o: schemes.c main.h roe.h schemes.h
	gcc -Wall -c schemes.c

timestep.o: timestep.c main.h timestep.h
	gcc -Wall -c timestep.c
