################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../datastructure/HashMap.cpp \
../datastructure/HashTable.cpp \
../datastructure/Matrices.cpp \
../datastructure/SparseMatrix.cpp \
../datastructure/UiSetPair.cpp \
../datastructure/Vector.cpp 

OBJS += \
./datastructure/HashMap.o \
./datastructure/HashTable.o \
./datastructure/Matrices.o \
./datastructure/SparseMatrix.o \
./datastructure/UiSetPair.o \
./datastructure/Vector.o 

CPP_DEPS += \
./datastructure/HashMap.d \
./datastructure/HashTable.d \
./datastructure/Matrices.d \
./datastructure/SparseMatrix.d \
./datastructure/UiSetPair.d \
./datastructure/Vector.d 


# Each subdirectory must supply rules for building sources it contributes
datastructure/%.o: ../datastructure/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Applications/Matlab/extern/include -I/Users/parastiwari/research/trnk/ipopt/src/include -I/Users/parastiwari/Ipopt3102/build/include/coin -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


