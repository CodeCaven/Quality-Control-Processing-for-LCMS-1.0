:: Name:     check_clayton.bat
:: Purpose:  Update the database for Clayton metabolomics
:: Author:   Simon Caven
:: Revision: 17/2/2019

@ECHO OFF
SETLOCAL ENABLEEXTENSIONS ENABLEDELAYEDEXPANSION
SET me=%~n0

IF EXIST "D:\QCScripts\clayton.txt" (
    ECHO %me% >> D:\QCScripts\update_log.txt
    ECHO Script already running.. >> D:\QCScripts\update_log.txt
    ECHO %DATE% >> D:\QCScripts\update_log.txt
    ECHO %TIME% >> D:\QCScripts\update_log.txt
    
) ELSE (
    :: create temp file
    type NUL > clayton.txt
	
	:: logging
    ECHO %me% >> D:\QCScripts\update_log.txt
    ECHO Starting update.. >> D:\QCScripts\update_log.txt
    ECHO %DATE% >> D:\QCScripts\update_log.txt
    ECHO %TIME% >> D:\QCScripts\update_log.txt
    
    :: activate conda
    CALL C:\ProgramData\Miniconda3\Scripts\activate.bat

    :: run script
    cd D:\QualityControlDatabase\Metabolomics-Quality-Control-Pipeline
    python Metabolomics_Process_Raw_Files.py "\\storage.erc.monash.edu\Shares\R-MNHS-MBPF\Shared\Metabolomics\QC_runs\C1_Clayton" "\\storage.erc.monash.edu\Shares\R-MNHS-MBPF\Shared\Metabolomics\QC_runs\QC_outfiles" "Clayton"

    :: delete the temp file
    del "D:\QCScripts\clayton.txt" /q
      
)

:END
ENDLOCAL
ECHO ON
@EXIT /B 0
