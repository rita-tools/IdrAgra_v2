@echo off
setlocal

rem Try current directory
for %%F in ("%~dp0IdrAgra*.exe") do (
    if exist "%%~fF" (
        "%%~fF"
        goto :done
    )
)

rem Try parent\release
for %%F in ("%~dp0..\release\IdrAgra*.exe") do (
    if exist "%%~fF" (
        "%%~fF"
        goto :done
    )
)

rem Try parent directory
for %%F in ("%~dp0..\IdrAgra*.exe") do (
    if exist "%%~fF" (
        "%%~fF"
        goto :done
    )
)

echo No IdrAgra executable was found.
echo Searched:
echo   %~dp0IdrAgra*.exe
echo   %~dp0..\release\IdrAgra*.exe
echo   %~dp0..\IdrAgra*.exe

:done
pause