@echo off
chcp 65001 >nul
setlocal
cd /d "%~dp0"

echo Installation des dépendances Python de l'interface Cavatappi Alpha V2...
set "PYTHON_CMD=python"
python --version >nul 2>&1
if errorlevel 1 (
    py -3 --version >nul 2>&1
    if errorlevel 1 (
        echo.
        echo Python est introuvable. Installez Python puis relancez ce fichier.
        pause
        exit /b 1
    )
    set "PYTHON_CMD=py -3"
)

%PYTHON_CMD% --version
%PYTHON_CMD% -m pip install --upgrade pip
if errorlevel 1 goto :install_error
%PYTHON_CMD% -m pip install -r "%~dp0requirements.txt"
if errorlevel 1 goto :install_error

echo.
echo Installation terminée.
pause
exit /b 0

:install_error
echo.
echo ERREUR : l'installation des dépendances a échoué.
pause
exit /b 1
