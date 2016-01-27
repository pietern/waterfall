DEFS=-D_POSIX_C_SOURCE=200809
override CFLAGS+=-std=c99 -Wall -Wextra -Wpedantic -g -ggdb $(shell pkg-config --cflags gtk+-3.0) $(DEFS)
override LDFLAGS+=-lm -lpthread -lfftw3 $(shell pkg-config --libs gtk+-3.0)
MAIN=waterfall
MAIN_OBJS=main.o stream.o window.o wedge.o color.o
STREAM_OBJS=$(patsubst %.c,%.o,$(wildcard stream_*.c))
OBJS=$(MAIN_OBJS) $(STREAM_OBJS)

.PHONY: all
all: $(MAIN)

.PHONY: clean
clean:
	rm -f $(MAIN) $(OBJS)

.PHONY: dep
dep:
	$(CC) -MM $(OBJS:.o=.c) > Makefile.dep

-include Makefile.dep

$(MAIN): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<
