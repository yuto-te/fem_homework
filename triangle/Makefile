FC       = gfortran
FFLAGS   =
LDFLAGS  =
LIBS     =
INCLUDE  = -I ./include
MOD_DIR  = -J ./include
BIN_DIR  = ./bin
OBJ_DIR  = ./build
SRC_DIR  = ./src
BIN_LIST = a.out
SRC_LIST = main.f90 matrix.f90 input.f90 solver.f90 post.f90
SOURCES  = $(addprefix $(SRC_DIR)/, $(SRC_LIST))
OBJECTS  = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
TARGET   = $(addprefix $(BIN_DIR)/, $(BIN_LIST))
RM       = rm

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

run: $(TARGET)
	$(TARGET)

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(TARGET) ./include/*.mod
