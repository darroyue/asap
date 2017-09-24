CC = g++
SDSLFLAGS = -O3 -DNDEBUG -I ~/include
LFLAGS = -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64
CFLAGS = -c -std=c++11 -Wall -Wno-missing-braces $(SDSLFLAGS)
EXEC = main
OBJS = main.o
MWT = m_wt.o
BITVECTORS = normal_bv.o sd_bv.o rrr_bv.o

$(EXEC): $(OBJS) $(BITVECTORS)
	$(CC) $^ -o $@ $(LFLAGS) $(SDSLFLAGS)

mapper: mapper.o $(BITVECTORS)
	$(CC) $^ -o $@ $(LFLAGS) $(SDSLFLAGS)

%.o: src/%.cpp
	$(CC) $< $(CFLAGS)

%.o: src/m_wt/%.cpp
	$(CC) $< $(CFLAGS)

%.o: src/bit_vectors/%.cpp
	$(CC) $< $(CFLAGS)

clean:
	find . -name '*.o' -type f -delete

cleanall:
	find . -name '*.o' -type f -delete
	rm -f $(EXEC)

run:
	./$(EXEC) input.dat