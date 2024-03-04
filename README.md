
you will need to download

1. VScode
2. Julia
3. julia ext in VScode
4. this package from github (manually download or clone whichever) put it in your user folder. the same on that has .vscode

then:

5. open FoilSim in vscode then shift enter first line or run the file, it should error
   
   5.5 Might happen automatically but at the bottom left of the VScode window select the folder sFoil as the environment if there are any issues
   
7. in julia kernel that will appear at the bottom type:
   
   ] activate "path"/sFoil  (etc. type folder path; this folder is called the environment)
   
   ] instantiate
   
8. now the file is runnable, I like to shift-enter each line and run blocks for testing but ctrl-enter should work to run the entire file



* include("") adds the files functions as if they were in the same file
