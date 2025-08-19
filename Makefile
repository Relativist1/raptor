# Makefile for RAPTOR Metal GPU Implementation
# For MacBook Pro M2 with Metal support

CC = clang
CFLAGS = -framework Foundation -framework Metal -framework MetalPerformanceShaders -framework MetalKit
TARGET = raptor_metal
SOURCE = raptor_metal.m

# Default target
all: $(TARGET)

# Build the Metal implementation
$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE)

# Clean build artifacts
clean:
	rm -f $(TARGET)

# Install (copy to /usr/local/bin)
install: $(TARGET)
	sudo cp $(TARGET) /usr/local/bin/

# Uninstall
uninstall:
	sudo rm -f /usr/local/bin/$(TARGET)

# Run the program
run: $(TARGET)
	./$(TARGET)

# Debug build with additional flags
debug: CFLAGS += -g -O0 -DDEBUG
debug: $(TARGET)

# Release build with optimization
release: CFLAGS += -O3 -DNDEBUG
release: $(TARGET)

# Profile build
profile: CFLAGS += -O2 -g -pg
profile: $(TARGET)

# Check Metal device capabilities
check-device: $(TARGET)
	./$(TARGET) --check-device

# Help
help:
	@echo "RAPTOR Metal GPU Implementation Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build the program (default)"
	@echo "  clean        - Remove build artifacts"
	@echo "  install      - Install to /usr/local/bin"
	@echo "  uninstall    - Remove from /usr/local/bin"
	@echo "  run          - Build and run the program"
	@echo "  debug        - Build with debug symbols"
	@echo "  release      - Build with optimization"
	@echo "  profile      - Build with profiling"
	@echo "  check-device - Check Metal device capabilities"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - macOS 10.15+ (Catalina)"
	@echo "  - Xcode Command Line Tools"
	@echo "  - Metal-capable GPU (M2 chip supported)"

.PHONY: all clean install uninstall run debug release profile check-device help