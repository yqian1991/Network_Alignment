################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/SemSim.cpp \
../src/Test.cpp \
../src/anemo.cpp \
../src/auxilliary.cpp \
../src/def.cpp \
../src/netaligner.cpp \
../src/networkblast.cpp \
../src/statistic.cpp 

OBJS += \
./src/SemSim.o \
./src/Test.o \
./src/anemo.o \
./src/auxilliary.o \
./src/def.o \
./src/netaligner.o \
./src/networkblast.o \
./src/statistic.o 

CPP_DEPS += \
./src/SemSim.d \
./src/Test.d \
./src/anemo.d \
./src/auxilliary.d \
./src/def.d \
./src/netaligner.d \
./src/networkblast.d \
./src/statistic.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


