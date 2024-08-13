@echo off
@echo Cleaning up directory:
set "curr_folder=%~dp0"
set "batch_folder=%~dp0/CUDA"
cd %batch_folder%
del *.exp *.lib *.o *.dll *.exe

@echo Making cuLUT.dll:
@REM I am not sure about what the "cudart" do
nvcc -c -o cuLUT.o cuLUT.cu -DDLL
nvcc -shared -o ../cuLUT.dll cuLUT.o -lcudart

@echo Making cuUnzip.dll:
nvcc -c -o cuUnzip.o cuUnzip.cu -DDLL
nvcc -shared -o ../cuUnzip.dll cuUnzip.o -lcudart

@echo (Note: objdump -- A good tool of gcc users to search for the addresses of functions)
@REM objdump -x cuLUT.dll

@echo Making main:
@echo (Note: the syntax of nvcc is different from gcc to link a lib)
@REM g++ main.cpp -o main.exe -L. -lparser
@REM nvcc main.cu -o main.exe -L. -llib_parser -rdc=true
set filename=main_gpu
nvcc main.cu -o ../%filename%.exe -rdc=true -std=c++17
cd %curr_folder%

@echo Cleaning up directory:
del *.exp *.lib *.o


@echo Testing: 
IF EXIST "example_data" (
    rd /s /q "example_data_unzip_curves"
    %filename% example_data 
)