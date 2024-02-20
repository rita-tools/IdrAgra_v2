# Readme #
## Project description ##
### IdrAgra - Idrologia Agraria (in Italian) ###
IdrAgra is a distributed-parameter agro-hydrological model that allows the simulation of the irrigation water distribution in agricultural areas and the estimation of the hydrological balance on a daily basis.
## Contact ##
prof. Claudio Gandolfi

email: claudio.gandolfi@unimi.it

## Contributors ##
* Anna Borghi %AB%
* Enrico Antonio Chiaradia %EAC% enrico.chiaradia@unimi.it
* Rachele Riva %RR%

## Conventions ##
Please refer to the following note for reading the code and apply edits.
### Variable names ####
* t_* = soil water content, theta, dimensionless  [m^3/m^3]
* h_* = water height/depth specific volume always in millimeter [mm]
* v_* = water volume in cubic meter [m^3]
* q_* = discharge in cubic meter per second [m^3/s]

## Compile instruction ##
### Windows ###
#### Step 1: install MSYS2 (once) ####
MSYS2 is a collection of tools and libraries providing you with an easy-to-use build system.
1. download the installer from https://www.msys2.org
2. run the installer, enter the desired installation folder, leave other settings as default
3. when installation finished, run MSYS2
4. type the following commands in the MSYS2 console:
```Shell
pacman -Syu
pacman -Su
pacman -S --needed base-devel mingw-w64-x86_64-toolchain
pacman -S mingw-w64-x86_64-python3-pip mingw-w64-x86_64-python3-setuptools
/mingw64/bin/pip3 install fortran-language-server
pacman -S mingw-w64-x86_64-lapack
pacman -S git unzip zsh 
```
> [!NOTE]
> After step 5, the terminal will be closed and you have to reopen it. Step 7 will install a complete building toolchain, other steps install useful packages
#### Step 2: install Visual Studio Code (once) ####
Install Visual Studio Code, VSC, (not Visual Studio!) from https://code.visualstudio.com/
#### Step 3: install Visual Studio Code Extensions (once) ####
Launch VSC. In the side bar (left), click the “Extensions” icon (in the middle). The following extensions of VSC should be (already) installed:
- C/C++ extension
- Makefile extension 
- Code Runner
- Modern Fortran
- Data Preview (to view and analyze csv files, optional)
#### Step 4: edit settings (once) ####
In the side bar (left), click the "Manage" icon (bottom) , and select "Settings".
In the settings, search for "code runner" and set "on" the following three items:
- Code runner: run in terminal -> ON
- Code runner: save all file before run -> ON
- Code runner: save file before run -> ON
- Code runner: ignore selection -> true

#### Step 5: get the code from github ####
Open a new command terminal (Start menu -> cmd)

Go to the desired root folder (where the project folder will be saved):
```Shell
cd YOUR-ROOT-PATH
```

Clone your repo with the following git command:

```Shell
git clone https://github.com/rita-tools/IdrAgra_v2
```

Change your terminal into that new subdirectory:

```Shell
cd IdrAgra_v2
```
Open Visual Studio Code  and select File -> Open Folder

#### Step 6: compile the executable ####

In the terminal, type the following:

```Shell
make cleanall
make
```
A brand-new idragra executable should be created!

#### Step 6bis: alternative, run in debug mode ####

Go to the launch.json file under .vscode folder and replace the path with that to your idragra dataset (where idragra_parameters.txt exists):
```json
...
"cwd": "YOUR-PATH-TO-IDRAGRA-DATASET",
...
make
```

On the side bar (left), select Run and Debug. Then, click on the green triangle on the top. If necessary, select the debugger (gdb) Launch
