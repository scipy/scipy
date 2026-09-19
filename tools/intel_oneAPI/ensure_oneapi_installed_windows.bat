REM Install the Intel oneAPI HPC kit unless it is already present (e.g. restored from cache).
if not exist "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" (
  call tools\intel_oneAPI\install_intel_oneAPI_windows.bat %WINDOWS_HPC_URL% - %1
  exit /b %ERRORLEVEL%
)
