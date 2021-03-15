CC      := gcc
CCFLAGS := -ggdb
LDFLAGS :=

TARGETS:= rsec
OBJ    := main.o RSErasureCode.o sha-256.o
DEPS   := sha-256.h RSErasureCode.h

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJ)

$(OBJ): %.o : %.c $(DEPS)
	$(CC) -c -o $@ $< $(CCFLAGS)

$(TARGETS): % : $(OBJ)
	$(CC) -o $@ $(LIBS) $^ $(CCFLAGS) $(LDFLAGS)

run: $(TARGETS)
	valgrind ./rsec
