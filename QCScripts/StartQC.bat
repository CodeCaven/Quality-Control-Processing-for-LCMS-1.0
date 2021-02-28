:: batch to start server if it is not running
:: and open the browser at the QC address

@ECHO OFF
SETLOCAL ENABLEEXTENSIONS ENABLEDELAYEDEXPANSION

set EXE=mysqld.exe
:: check for running task, could check for service instead 
FOR /F %%x IN ('tasklist /NH /FI "IMAGENAME eq %EXE%"') DO IF %%x == %EXE% goto Node
goto NoDB

:Node
set EXE=node.exe
FOR /F %%x IN ('tasklist /NH /FI "IMAGENAME eq %EXE%"') DO IF %%x == %EXE% goto ProcessFound
goto ProcessNotFound

:ProcessFound
goto OPEN

:ProcessNotFound
cd C:\Users\sliggady\OneDrive\FITHonoursProject\AutomatedMSandLCPipeline\DEV\MPMF-Server\mpmf_qc_server
start npm start
::SET DEBUG= mbpf-server:* & npm run devstart
goto OPEN

:OPEN
cd C:\Program Files (x86)\Google\Chrome\Application
chrome "http://localhost:3000"
goto END

:NoDB
echo The MySQl database service is not running
echo Start the service (or wait for task scheduler to do it)
echo And Try Again
@pause
ENDLOCAL
ECHO ON
@EXIT /B 1

:END
ENDLOCAL
ECHO ON
@EXIT /B 1

