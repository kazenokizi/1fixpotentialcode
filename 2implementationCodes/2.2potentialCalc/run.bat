call "C:\Program Files (x86)\Intel\Compiler\11.0\066\fortran\Bin\ifortvars.bat" ia32
ifort main.for -o main.exe
main.exe
del main.obj
del main.exe