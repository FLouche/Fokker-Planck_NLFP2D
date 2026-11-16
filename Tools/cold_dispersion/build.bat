@echo off
rem Build cold_dispersion.exe with Intel ifx (oneAPI 2025.2) from a plain cmd prompt.
rem Uses the VS 2022 x64 linker environment plus the ifx directories directly
rem (oneAPI's setvars.bat fails when vswhere.exe is not on the PATH).
setlocal
cd /d "%~dp0"
set ONEAPI=C:\Program Files (x86)\Intel\oneAPI\compiler\2025.2
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
set "PATH=%ONEAPI%\bin;%PATH%"
set "LIB=%ONEAPI%\lib;%LIB%"
set "INCLUDE=%ONEAPI%\include;%INCLUDE%"
ifx /nologo /O2 /warn:all /exe:cold_dispersion.exe cold_dispersion.f90
if errorlevel 1 exit /b 1
del /q cold_dispersion.obj *.mod 2>nul
echo built %~dp0cold_dispersion.exe
