****************************************************************************                             
			    HRCM
  (Hybrid Referential Compression Method )

     https://github.com/haicy/HRCM

          Copyright (C) 2019                  
****************************************************************************

1. Introduction

1.1 HRCM is implemented with C++.

1.2 For details of methods, please refer to the original paper.

1.3 This tool is implemented with linux operating system, you can run the tool on Linux operating system, other operating system is not guaranteed.

1.4 Open source tool 7-zip is used in HRCM to provide the implementation of PPMD encoding. 

****************************************************************************

2. Compile

2.1 Download: 
You can find and download source code and test datasets at:
https://github.com/haicy/HRCM

2.2 System requirements:
A C/C++ compiler and the STL library. In the Makefile, we use g++ as the default compiler command.

2.3 Compiling: 　
Run "make" in the HRCM source code directory, the executable file hrcm is in the current directory.
Run "chmod 777 7za" to change the access permission of the 7-zip executables.

****************************************************************************

3 Use

3.1 Usage
hrcm {compress | decompress}  -r {ref-file-path}{ [-t] {tar-file-path}|[-f] {filename} [percent]}
     {compress | decompress} is mode,  choose one of them according to requirement, required
     -r is the reference, the {ref-file-path} followed, required
     -t is the target, a single to-be-compressed file path {tar-file-path} followed, optional
     -f is the alternative option of -t, a set of to-be-compressed file paths included in {filename}, optional
     [percent] is the percentage of the second-level matching, default is 10, means 10% of sequences will be used for the second-level matching, optional when -f, illegal when -t

3.2 Output:
1 compressed file named filename.7z 
2 decompressed file named filename.fasta

****************************************************************************

4 Example

4.1 compress and decompress hg18_chr22.fa using hg17_chr22.fa as reference
    ./hrcm compress -r hg17_chr22.fa -t hg18_chr22.fa
    output: hg18_chr22.7z
	
    ./hrcm decompress -r hg17_chr22.fa -t hg18_chr22.7z
    output: hg18_chr22.fasta
	
4.2 compress and decompress hg17_chr22.fa and hg18_chr22.fa, using hg13_chr22.fa as reference. The to-be-compressed paths are written in chr22.txt
    ./hrcm compress -r hg13_chr22.fa -f chr22.txt 100
    output: chr22.7z


    ./hrcm decompress -r hg13_chr22.fa -f chr22.txt 100
    output: hg17_chr22.fasta 
	    hg18_chr22.fasta
    The decomprssed file path is specified in chr22.txt.

4.3 check the difference between original file and decompressed file
    diff hg17_chr22.fa hg17_chr22.fasta
    diff hg18_chr22.fa hg18_chr22.fasta

***************************************************************************
