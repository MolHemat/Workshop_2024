# Workshop 2024 bash tutorial
Commands to be entered are highlighted as below.
```
command
``` 
Press "Enter" after you type the command  
In order to understand the usage of a command type ```command --help```

## 1. Login to your account
Open MobaXterm ("Click on the Windows icon on the lower left corner, search for MobaXterm, click on the MobaXterm icon")
Type the follwing command to login to your account (username@ipaddress)
```
ssh ciw12@10.100.75.81
```
Enter password
## 2. Find your current location (home directory)
```
pwd
```
## 3. List the files present in the current directory
```
ls 
```
## 4. Make a temporary directory for output
```
mkdir temp
```
## 5. Enter into the directory made in above step
```
cd temp
```
## 6. Make an empty file
```
touch myfile.txt
``` 
## 7. Add a line to your empty file
```
echo "this is my first line" >> myfile.txt
```
## 8. Add another line by opening the file in a text editor (notepad, wordpad etc)
Go to the left handside panel on MobaXterm, double click on the "temp" directory.  
Right click on the myfile.txt, => Open with , select Wordpad / Notepad  
Add a new line saying "This line was added using wordpad"  
Save it by pressing "Ctrl + S"  
Close the file

## 9. Check the contents of your file
```
cat myfile.txt
```
## 10. Count the occurences of the word "line" in your file
```
grep -c 'line' myfile.txt
```
## 11. How many lines are present in your file ?
```
wc myfile.txt
```
To understand the output type 
```
wc --help
```
## 12. Find the number of reads in a fastq.gz file
### Location of the fastq.gz file
```
source path_files.dat
echo ${sequences}
```
The above ```source ```command helps to assign values to the variables.  
In our case, the file path_files.dat contains the location of files / programs assigned to different variables  
Check the details by opening the file using wordpad or ```cat``` command  
Find the location of input files. Check if files are present by doing a list  
```
ls `path to sequences/file name`
```
### Count the number of lines using wc 
We will use ```zcat``` instead of ```cat``` to open the file since the file has .gz extension.  
It stands for a compressed file
```
zcat `path to sequences/file name` | wc 
```
The above value denotes the total number of lines.  
Total no. of reads in a fastq.gz file = Total number of lines / 4  
You can perform the above calculation in linux using ```bc``` command
```
echo (Total no. of lines / 4) | bc
```