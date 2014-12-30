# If you use the random number generator sprng() and the files are 
# installed in places the compiler searches, then you should not need
# to modify this file. But if the compiler complains 
# that it cannot find "sprng.h" or "libcmrg.a" or "-lcmrg", then 
# you will need to modify the path /home/tgk/sprng/include to the 
# path that contains the header file for sprng, and modify the path 
# /home/tgk/sprng/lib to the path that contains the executable file for sprng. 

lib1.o: src/public/lib1.c src/public/local.h src/public/variable.h; \
 g++ -I/home/tgk/sprng/include -L/home/tgk/sprng/lib -O3 -Wall  \
 src/public/lib1.c -c -o lib1.o

simul: src/public/simul.c lib1.o src/public/local.h src/public/variable.h; \
 g++ -I/home/tgk/sprng/include -L/home/tgk/sprng/lib -O3 -Wall  \
 src/public/simul.c lib1.o -lm -lcmrg -o simul


