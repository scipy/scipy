REM Activate the Intel oneAPI and MSVC environments, then run the given command.
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
%*
exit /b %ERRORLEVEL%
