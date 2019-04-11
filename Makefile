FC       = gfortran
INCLUDE  = -I ./include
MOD_DIR  = -J ./include
BIN_DIR  = ./bin
OBJ_DIR  = ./build
SRC_DIR  = ./src
SUFFIX = .f90
SOURCES:=$(wildcard $(SRC_DIR)/*$(SUFFIX))
OBJECTS:=$(addprefix $(OBJ_DIR)/, $(notdir $(SOURCES:$(SUFFIX)=.o)))
TARGETS:=$(addprefix $(BIN_DIR)/, $(notdir $(basename $(SOURCES))))

define OBJECT
	$(addprefix $(OBJ_DIR)/, $(notdir $(1)).o)
endef

define SOURCE
	$(addprefix $(SRC_DIR)/, $(notdir $(1))$(SUFFIX))
endef

define MAKEALL
$(1): $(call OBJECT, $(1))
	$(FC) -o $(1) $(call OBJECT, $(1))

$(addprefix $(OBJ_DIR)/, $(notdir $(1)).o): $(call SOURCE, $(1))
	$(FC) $(INCLUDE) $(MOD_DIR) -o $(call OBJECT, $(1)) -c $(call SOURCE, $(1))
endef

.PHONY: all
all: $(TARGETS)
$(foreach VAR,$(TARGETS),$(eval $(call MAKEALL,$(VAR))))

.PHONY:clean
clean:
	$(RM) $(OBJECTS) $(TARGETS) ./include/*.mod
