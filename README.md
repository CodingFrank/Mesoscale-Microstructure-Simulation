PPC Group Project

To compile:
 - make all

To run a job, use sbatch command
 - sbatch job1.sh

To convert binary data file into a matrix of integer file, use hexdump command,
with N being the y dimension size.
 - hexdump -v -e 'N/4 "%10d "' -e '"\n"' outputfile > outputfile_ascii

To plot the .dat file, use the matlab function recolor.m and plotSnapshots.m

