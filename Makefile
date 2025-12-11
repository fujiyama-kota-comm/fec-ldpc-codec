CC      = gcc
CFLAGS  = -O2 -Wall -std=gnu11 -Iinclude
LDFLAGS = -lm

# ============================================================
# Sources
# ============================================================
SRC = \
    src/ldpc_matrix.c \
    src/ldpc_encoder.c \
    src/ldpc_decoder.c

OBJ = $(SRC:.c=.o)

# Example program
GENE_HG_SRC = mains/gene_hg.c
GENE_HG_OBJ = $(GENE_HG_SRC:.c=.o)

# BER program
LDPC_BER_SRC = mains/ldpc_ber.c
LDPC_BER_OBJ = $(LDPC_BER_SRC:.c=.o)

# Output dir
BIN_DIR = bin

# OS-dependent executables
ifeq ($(OS),Windows_NT)
    GENE_HG_TARGET = $(BIN_DIR)/gene_hg.exe
    LDPC_BER_TARGET = $(BIN_DIR)/ldpc_ber.exe
    RUN_GENE_HG = $(GENE_HG_TARGET)
    RUN_LDPC_BER = $(LDPC_BER_TARGET)
else
    GENE_HG_TARGET = $(BIN_DIR)/gene_hg
    LDPC_BER_TARGET = $(BIN_DIR)/ldpc_ber
    RUN_GENE_HG = ./$(GENE_HG_TARGET)
    RUN_LDPC_BER = ./$(LDPC_BER_TARGET)
endif

# ============================================================
# Build rules
# ============================================================
all: $(GENE_HG_TARGET) $(LDPC_BER_TARGET)

# Create bin directory
$(BIN_DIR):
	@if [ ! -d "$(BIN_DIR)" ]; then \
		mkdir -p $(BIN_DIR) 2>/dev/null || mkdir $(BIN_DIR); \
	fi

$(GENE_HG_TARGET): $(BIN_DIR) $(OBJ) $(GENE_HG_OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(GENE_HG_OBJ) $(LDFLAGS)

$(LDPC_BER_TARGET): $(BIN_DIR) $(OBJ) $(LDPC_BER_OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDPC_BER_OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# ============================================================
# Run commands
# ============================================================
gene_hg: $(GENE_HG_TARGET)
	$(RUN_GENE_HG)

ldpc_ber: $(LDPC_BER_TARGET)
	$(RUN_LDPC_BER)

# ============================================================
# Clean
# ============================================================
clean:
	@echo "Cleaning object files..."
	rm -f $(OBJ) $(GENE_HG_OBJ) $(LDPC_BER_OBJ)

	@echo "Cleaning binaries..."
	@if [ -f "$(GENE_HG_TARGET)" ]; then rm -f "$(GENE_HG_TARGET)"; fi
	@if [ -f "$(LDPC_BER_TARGET)" ]; then rm -f "$(LDPC_BER_TARGET)"; fi

	@if [ -d "$(BIN_DIR)" ] && [ ! "$$(ls -A $(BIN_DIR))" ]; then \
		echo "Removing empty bin directory"; \
		rmdir $(BIN_DIR); \
	fi

.PHONY: all clean gene_hg ldpc_ber
