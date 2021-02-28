

@ECHO OFF
SETLOCAL ENABLEEXTENSIONS ENABLEDELAYEDEXPANSION
SET me=%~n0

:: logging
ECHO %me% >> D:\QCScripts\server_log.txt
ECHO %USERNAME% >> D:\QCScripts\server_log.txt
ECHO %DATE% >> D:\QCScripts\server_log.txt
ECHO %TIME% >> D:\QCScripts\server_log.txt

:: start server
cd D:\QualityControlServer\mpmf_qc_server  
npm start



:END
ENDLOCAL
ECHO ON
@EXIT /B 0
