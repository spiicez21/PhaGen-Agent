@echo off
echo Clearing Python cache files...
for /d /r . %%d in (__pycache__) do @if exist "%%d" rd /s /q "%%d"
echo Clearing pytest cache...
if exist ".pytest_cache" rd /s /q ".pytest_cache"
echo Done.
pause
