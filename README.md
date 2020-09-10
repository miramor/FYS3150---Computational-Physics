# FYS3150
FYS3150 – Computational Physics .  Filer til undervisning og oppgaver. Prosjekter.

## Projects
### Project 1:
Thomas algorithm and LU decomposition to solve second order differential equation. Writing the results to files for different n values and visualizing results.

**How to use:**

Create a directory called **ResultsComputation** in project 1 folder, this houses all the results.
Run make all for useSpecial = true once and then again with  useSpecial = false. This does all the calculations and writes the results to the folder mentioned earlier. (Can comment out the part which handles LU calculations, 2nd part of main, if already ran main.cpp once to avoid excessive calculations). Adjust repeat paramater if needed (low = fast and less precision *** high = slow and more precice time). 

When all results for the different n values are written. Then run visualize.py if you wanna see tables which shows things such as error and results for all methods for *n = 10*. Also creates some plots which is placed in a directoray to visualize the solution vs exact curves.

