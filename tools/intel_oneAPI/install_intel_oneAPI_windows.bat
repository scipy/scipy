REM Reference - https://github.com/oneapi-src/oneapi-ci/blob/master/scripts/install_windows.bat
REM SPDX-FileCopyrightText: 2022 Intel Corporation
REM
REM SPDX-License-Identifier: MIT

setlocal enabledelayedexpansion

set URL=%~1
set COMPONENTS=%~2
set EXPECTED_SHA256=%~3

:: download installer from intel
curl.exe --output %TEMP%\webimage.exe --url %URL% --retry 5 --retry-delay 5

:: verify checksum if provided
if not "%EXPECTED_SHA256%"=="" (
  for /f "tokens=1" %%h in ('certutil -hashfile %TEMP%\webimage.exe SHA256 ^| findstr /v "hash"') do set FILE_HASH=%%h
  if /i not "!FILE_HASH!"=="%EXPECTED_SHA256%" (
    echo Checksum verification failed for %URL%
    echo Expected: %EXPECTED_SHA256%
    echo Got:      !FILE_HASH!
    del %TEMP%\webimage.exe
    exit /b 1
  )
)

start /b /wait %TEMP%\webimage.exe -s -x -f webimage_extracted
del %TEMP%\webimage.exe
if "%COMPONENTS%"=="" (
  webimage_extracted\bootstrapper.exe -s --action install --eula=accept -p=NEED_VS2017_INTEGRATION=0 -p=NEED_VS2019_INTEGRATION=0 -p=NEED_VS2022_INTEGRATION=0
) else (
  webimage_extracted\bootstrapper.exe -s --action install --components=%COMPONENTS% --eula=accept -p=NEED_VS2017_INTEGRATION=0 -p=NEED_VS2019_INTEGRATION=0 -p=NEED_VS2022_INTEGRATION=0
)
set installer_exit_code=%ERRORLEVEL%
rd /s/q "webimage_extracted"
exit /b %installer_exit_code%
