:: Name:     check_clayton.bat
:: Purpose:  Update the database for Clayton proteomics
:: Author:   Simon Caven
:: Revision: 27/10/2019

@ECHO OFF
SETLOCAL ENABLEEXTENSIONS ENABLEDELAYEDEXPANSION
SET me=%~n0

IF EXIST "D:\mpmf_qc_scripts\parkville_metabolomics.txt" (
    ECHO %me% >> D:\mpmf_qc_scripts\update_log.txt
    ECHO Script already running.. >> D:\mpmf_qc_scripts\update_log.txt
    ECHO %DATE% >> D:\mpmf_qc_scripts\update_log.txt
    ECHO %TIME% >> D:\mpmf_qc_scripts\update_log.txt
    
) ELSE (
    :: create temp file
    type NUL > parkville_metabolomics.txt
	
	:: logging
    ECHO %me% >> D:\mpmf_qc_scripts\update_log.txt
    ECHO Starting update.. >> D:\mpmf_qc_scripts\update_log.txt
    ECHO %DATE% >> D:\mpmf_qc_scripts\update_log.txt
    ECHO %TIME% >> D:\mpmf_qc_scripts\update_log.txt
    
    :: activate conda
    CALL C:\ProgramData\Miniconda3\Scripts\activate.bat

    :: run script
    cd D:\Processing-Quality-Control-Pipeline
    python MPMF_Process_Raw_Files.py "\\storage.erc.monash.edu\Shares\R-MNHS-MBPF\Shared\Metabolomics\QC_runs" "\\storage.erc.monash.edu\Shares\R-MNHS-MBPF\Shared\Metabolomics\QC_runs\QC_outfiles" "parkville" "metabolomics"

    :: delete the temp file
    del "D:\mpmf_qc_scripts\parkville_metabolomics.txt" /q
      
)

:END
ENDLOCAL
ECHO ON
@EXIT /B 0
