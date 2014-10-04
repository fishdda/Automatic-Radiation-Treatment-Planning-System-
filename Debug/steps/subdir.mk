################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../steps/StepIIImpl.cpp \
../steps/StepIImpl.cpp \
../steps/StepIVmpl.cpp \
../steps/StepImpl.cpp 

OBJS += \
./steps/StepIIImpl.o \
./steps/StepIImpl.o \
./steps/StepIVmpl.o \
./steps/StepImpl.o 

CPP_DEPS += \
./steps/StepIIImpl.d \
./steps/StepIImpl.d \
./steps/StepIVmpl.d \
./steps/StepImpl.d 


# Each subdirectory must supply rules for building sources it contributes
steps/%.o: ../steps/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Applications/Matlab/extern/include -I/Users/parastiwari/research/trnk/ipopt/src/include -I/Users/parastiwari/Ipopt3102/build/include/coin -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


