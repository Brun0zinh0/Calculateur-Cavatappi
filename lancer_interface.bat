@echo off
chcp 65001 >nul
setlocal
cd /d "%~dp0"

echo Lancement de l'interface Cavatappi...
set "PYTHON_CMD=python"
python --version >nul 2>&1
if errorlevel 1 (
    py -3 --version >nul 2>&1
    if errorlevel 1 (
        echo.
        echo ERREUR : Python est introuvable. Lancez d'abord installer_dependances.bat.
        pause
        exit /b 1
    )
    set "PYTHON_CMD=py -3"
)

%PYTHON_CMD% -c "import sys; sys.path.insert(0, sys.argv[1]); import parametres as p; assert hasattr(p, 'MAXWELL_ANISOTROPY_LABELS'), 'parametres.py est incomplet'" "%~dp0"
if errorlevel 1 (
    echo.
    echo ERREUR : les fichiers du dossier Alpha V2 sont incomplets ou incompatibles.
    echo Verifiez que lancer_interface.bat, interface.py et parametres.py proviennent de la meme version.
    pause
    exit /b 1
)

set "APP_PORT="
for /f "delims=" %%P in ('%PYTHON_CMD% -c "import socket; s=socket.socket(); s.bind(('127.0.0.1', 0)); print(s.getsockname()[1]); s.close()"') do set "APP_PORT=%%P"
if not defined APP_PORT (
    echo.
    echo ERREUR : aucun port local libre n'a pu etre determine.
    pause
    exit /b 1
)

echo Adresse locale : http://localhost:%APP_PORT%
%PYTHON_CMD% -m streamlit run "%~dp0interface.py" --server.port %APP_PORT%

pause
