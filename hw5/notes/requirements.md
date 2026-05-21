Hello everyone,

For homework 5, the declarations/returns of the functions are clearly laid out in the homework pdf. Here I will address the structure of your submission.

Your zipped file should contain TWO folders:

1) Matrices

This folder should contain your solution to Problem 1 with all required header/c/cpp files/makefile. Your makefile should create an executable "matrices" that depends on your main and other files. Note that this is essentially the Matrices folder provided to you (without the extra .txt files).

2) MatricesGSL

This folder should contain your solution to Problem 2 with all required header/c/cpp files/makefile. Your makefile should create an executable "matrices" that depends on your main and other files. Note that this is essentially the MatricesGSL folder provided to you (without the extra .txt files).

If you encountered a module loading issue when compiling the second question in HW 5. The common solution is to load other libraries that your scripts need. For example:

module load GSL

module load R

Please make sure that your makefiles have the commands for "all" and "clean" appropriately as shown in class/given to you. Before running your code, I will "make clean" and "make all" on the cluster. PLEASE make sure that your homework compiles AND runs ON THE CLUSTER. You will get no credit if I cannot compile your code on the cluster by running "make clean" then "make all". Please read that last sentence again.
