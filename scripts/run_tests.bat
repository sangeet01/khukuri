@echo off
REM Run Khukuri tests on Windows

echo Running Khukuri test suite...

pytest tests\ -v

echo.
echo Running with coverage...
pytest tests\ --cov=src --cov-report=html --cov-report=term

echo.
echo Tests complete! Coverage report in htmlcov\index.html
pause
