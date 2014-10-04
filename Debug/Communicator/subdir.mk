################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Communicator/DummyMatlabComm.cpp \
../Communicator/FileManager.cpp \
../Communicator/MatlabCommunicationManager.cpp \
../Communicator/OptimizationImpl.cpp \
../Communicator/PrioritizedMain.cpp \
../Communicator/PrioritizedNlp.cpp \
../Communicator/ProblemReps.cpp \
../Communicator/Result.cpp \
../Communicator/Util.cpp \
../Communicator/oneVoxel.cpp 

OBJS += \
./Communicator/DummyMatlabComm.o \
./Communicator/FileManager.o \
./Communicator/MatlabCommunicationManager.o \
./Communicator/OptimizationImpl.o \
./Communicator/PrioritizedMain.o \
./Communicator/PrioritizedNlp.o \
./Communicator/ProblemReps.o \
./Communicator/Result.o \
./Communicator/Util.o \
./Communicator/oneVoxel.o 

CPP_DEPS += \
./Communicator/DummyMatlabComm.d \
./Communicator/FileManager.d \
./Communicator/MatlabCommunicationManager.d \
./Communicator/OptimizationImpl.d \
./Communicator/PrioritizedMain.d \
./Communicator/PrioritizedNlp.d \
./Communicator/ProblemReps.d \
./Communicator/Result.d \
./Communicator/Util.d \
./Communicator/oneVoxel.d 


# Each subdirectory must supply rules for building sources it contributes
Communicator/%.o: ../Communicator/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Applications/Matlab/extern/include -I/Users/parastiwari/research/trnk/ipopt/src/include -I/Users/parastiwari/Ipopt3102/build/include/coin -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


