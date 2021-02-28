:: batch to check if MySQL service is running and run it if not
:: called from task scheduler every 5 minutes

@ECHO OFF
SETLOCAL ENABLEEXTENSIONS ENABLEDELAYEDEXPANSION

set EXE=mysqld.exe
:: check for running task, could check for service instead 
FOR /F %%x IN ('tasklist /NH /FI "IMAGENAME eq %EXE%"') DO IF %%x == %EXE% goto ProcessFound

goto ProcessNotFound

:ProcessFound
echo mysql running
goto END

:ProcessNotFound
echo mysql not running
:: start the MySQL windows service
NET start MySQL80
goto END

:END
ENDLOCAL
ECHO ON
@EXIT /B 1

