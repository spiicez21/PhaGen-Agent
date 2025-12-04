@echo off
rem Proxy OSRA calls through Ubuntu WSL so Windows Python can invoke the Linux binary.
wsl -d Ubuntu -- osra %*
