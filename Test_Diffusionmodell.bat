xcopy "w:\user\01-Sadowski Group\Borrmann\diffusionpy\" %USERPROFILE%\diffusionpy\ /Y /E /D
xcopy %USERPROFILE%\diffusionpy\.matplotlib\ /Y /E /D %USERPROFILE%\.matplotlib
@echo off
echo 
echo 
echo The source code was copied to %USERPROFILE%\diffusionpy\ .
echo 
echo Your version uses a python installation on w:\user\"01-Sadowski Group"\Borrmann\Python311
echo 
echo If you want to install this package on your local Python Version check out the installation instruction on https://github.com/Borrdom/diffusionpy/
echo 
w:\user\"01-Sadowski Group"\Borrmann\Python311\python.exe -m pip install --user -e %USERPROFILE%\diffusionpy\

w:\user\"01-Sadowski Group"\Borrmann\Python311\python.exe -m pip install --user jupyter
%appdata%\Python\Python311\Scripts\jupyter-lab