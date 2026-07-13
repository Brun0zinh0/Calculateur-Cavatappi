@echo off
chcp 65001 >nul
setlocal
cd /d "%~dp0"

echo Installation des dépendances Python du livrable Cavatappi...
python --version
if errorlevel 1 (
    echo.
    echo Python est introuvable. Installe Python puis relance ce fichier.
    pause
    exit /b 1
)

python -m pip install --upgrade pip
python -m pip install -r requirements.txt

echo.
echo Installation terminée.
pause
