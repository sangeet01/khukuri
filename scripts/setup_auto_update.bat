@echo off
REM Setup automatic daily updates via Windows Task Scheduler

echo Setting up automatic daily data updates...

set SCRIPT_DIR=%~dp0
set PROJECT_DIR=%SCRIPT_DIR%..
set PYTHON_SCRIPT=%PROJECT_DIR%\scripts\update_data.py

REM Create task to run daily at 2:00 AM
schtasks /create /tn "KhukuriDataUpdate" /tr "python %PYTHON_SCRIPT%" /sc daily /st 02:00 /f

if %errorlevel% equ 0 (
    echo.
    echo ✓ Task scheduled: Daily update at 2:00 AM
    echo   Task name: KhukuriDataUpdate
    echo.
    echo To view tasks: schtasks /query /tn KhukuriDataUpdate
    echo To remove: schtasks /delete /tn KhukuriDataUpdate
    echo.
    
    REM Create logs directory
    if not exist "%PROJECT_DIR%\logs" mkdir "%PROJECT_DIR%\logs"
    echo ✓ Setup complete!
) else (
    echo.
    echo ✗ Failed to create scheduled task
    echo   Run this script as Administrator
)

pause
