# ==========================================================
# Silicon Virovore
# Build System
# ==========================================================

CC = clang

CFLAGS = -O3 -Wall -fPIC

SRC = \
	c/engine.c \
	c/hydropathy.c \
	c/ga_loop.c \
	c/safety.c

LIB = c/libsafety.so

all: $(LIB)

$(LIB): $(SRC)
	$(CC) $(CFLAGS) -shared $(SRC) -o $(LIB)

clean:
	rm -f $(LIB)
	rm -rf build

run:
	python3 run_pipeline.py

rebuild: clean all