#Compilers
CC = gcc
ARM_CROSSCOMPILER = arm-linux-gnueabihf-gcc

#Libraries to use
LDFLAGS = -lm -w -lpthread

#Folders definition
SRCS = $(wildcard src/*.c) $(wildcard WFA2/*/*.c)
ECHO = echo

#Check if compilers are installed
EXECUTABLES = gcc arm-linux-gnueabihf-gcc
K := $(foreach exec,$(EXECUTABLES),\
        $(if $(shell which $(exec)),some string,$(error "$(exec) Not Installed")))


all: x86 arm success

x86:
	@ $(CC) -Iinclude -o micromapX86.elf $(SRCS) $(LDFLAGS)
	
arm:
	@ $(ARM_CROSSCOMPILER) -Iinclude -o micromapARM.elf $(SRCS) $(LDFLAGS) -static
	
success:
	@$(ECHO) "Build succeeded"	


help:
	@$(ECHO) ""	
	@$(ECHO) "Possible options:"
	@$(ECHO) ""
	@$(ECHO) "all      -   compile x86 and ARM"
	@$(ECHO) "clean    -   cleanup old .o and .elf"
	@$(ECHO) "x86 	 -   compile using \033[0;32m$(CC)\033[0m " 
	@$(ECHO) "arm 	 -   cross-compile using \033[0;32m$(ARM_CROSSCOMPILER)\033[0m "
	@$(ECHO) ""

clean:
	# @$(ECHO) "  Cleanup done"	
	@ rm -rf *.o *.elf | clear
