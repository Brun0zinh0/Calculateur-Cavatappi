@echo off
chcp 65001 >nul
setlocal
cd /d "%~dp0"

echo Lancement de l'interface Cavatappi...
python -m streamlit run "%~dp0livrable_interface.py"

pause
