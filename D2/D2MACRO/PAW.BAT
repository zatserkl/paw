@echo off
if %PAW% == N goto LOAD
	cls
        echo.
        echo PAW already loaded.
        echo.
        echo N No load
        echo Y Load again
        echo.
        choice /c:NY Choice an option 

        if errorlevel 2 goto LOAD
        if errorlevel 1 goto QUIT

:LOAD
SET PAW=Y

call c:\paw\pawenv

if not exist last.kumac goto RUN
move /Y last.kumac c:\temp\last_old.kumac

:RUN
%HOME%\go32 %HOME%\pawdos.exe %1 %2 %3 %4 %5 %6 %7 %8 %9
rem SET HOME=
SET PAW=N

:QUIT
