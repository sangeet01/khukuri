@echo off
REM Install Khukuri dependencies on Windows

echo Installing Khukuri Virtual Lab dependencies...

python -m pip install --upgrade pip
pip install -r requirements.txt
pip install pytest pytest-cov black flake8

echo.
echo Dependencies installed successfully!
echo Run 'pytest' to verify installation
pause
