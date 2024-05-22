@if "%1" == "" goto NOTHING
echo *%1 >> t670.dat
copy t670.dat + %1
@echo off
del %1
goto END

:NOTHING
@echo File did not specified. Nothing to do!
@pause
@echo off

:END
