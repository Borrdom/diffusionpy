@echo off
chcp 65001 >nul &
for /f "delims=;" %%. in ('"prompt $H; & for %%. in (nul) do call"') do (set "bkspc=%%.")

set "ASCII219=██"

xcopy . %USERPROFILE%\diffusionpy\ /Y /E /D >NUL
xcopy %USERPROFILE%\diffusionpy\.matplotlib\ /Y /E /D %USERPROFILE%\.matplotlib >NUL
echo Loading, please wait. Loading for the first time might take one minute ...
<nul set /p "=%ASCII219%"||ver>nul

..\Python311\python.exe  -m pip list | findstr /i "diffusionpy" > nul
if %errorlevel% equ 0 (
    <nul set /p "=%ASCII219%"||ver>nul
) else (
    <nul set /p "=%ASCII219%"||ver>nul
    ..\Python311\python.exe -m pip install --user --quiet --no-warn-script-location -e %USERPROFILE%\diffusionpy\
)

<nul set /p "=%ASCII219%"||ver>nul

..\Python311\python.exe  -m pip list | findstr /i "jupyter" > nul
if %errorlevel% equ 0 (
    <nul set /p "=%ASCII219%"||ver>nul
) else (
    <nul set /p "=%ASCII219%"||ver>nul
    ..\Python311\python.exe -m pip install --user --quiet --no-warn-script-location jupyter
)
<nul set /p "=%ASCII219%"||ver>nul
echo :
echo Loading complete. Starting browser window. Please check out the folder example_notebooks on the left of the file explorer in the browser ...
cd %USERPROFILE%\diffusionpy\ 
%appdata%\Python\Python311\Scripts\jupyter-lab %USERPROFILE%\diffusionpy\  >NUL