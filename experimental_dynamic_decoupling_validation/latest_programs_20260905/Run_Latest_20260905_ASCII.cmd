@echo off
setlocal
rem Run this launcher from the package directory without moving it.
rem Z: is only an ASCII alias used to avoid MATLAB/HDF5 Unicode-path issues.
set "PKG=%~dp0"
subst Z: "%PKG%" >nul 2>&1
if errorlevel 1 (
  echo Failed to create Z: mapping.
  exit /b 1
)
matlab -batch "cd('Z:\'); Check_Latest_Program_20260905(); Run_R5_Main_20251222();"
set "RC=%ERRORLEVEL%"
subst Z: /d >nul 2>&1
exit /b %RC%
