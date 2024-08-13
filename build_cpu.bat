@ echo off
set /P LUT_SIZE="Enter look up table size (Default = 200): "
if "%LUT_SIZE%" == "" set "LUT_SIZE=200"

for %%f in (*.o *.lib *.exe) do (
    del "%%f"
)


set start_time=%TIME%
echo Compiling started at %start_time%

echo Building...
REM executable
cd CPU
g++ -std=c++20 -fconstexpr-ops-limit=100000000000 main.cpp -DJ_SIZE=%LUT_SIZE% -DEXT_SIZE=%LUT_SIZE% -o ../main_cpu.exe
cd ..

set end_time=%TIME%
for /f "tokens=1-3 delims=:." %%a in ("%start_time%") do (
    set /a "start_seconds=(((%%a*60)+%%b)*60)+%%c"
)
for /f "tokens=1-3 delims=:." %%a in ("%end_time%") do (
    set /a "end_seconds=(((%%a*60)+%%b)*60)+%%c"
)

echo Compiling finished at %end_time%

REM Calculate time difference in seconds
set /a "time_difference=end_seconds-start_seconds"
set /a "hours=time_difference/3600"
set /a "minutes=(time_difference%%3600)/60"
set /a "seconds=time_difference%%60"
set "time_difference=%hours% hour %minutes% minute %seconds% second."
echo Compiled in: %time_difference% 

echo Testing the program..
main_cpu.exe "example_data/Accessory colonization factor AcfD.txt" out.csv

