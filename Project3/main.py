import sys
import os
import matplotlib.pyplot as plt

try:
    if sys.argv[1] == "test":
        os.system("echo compiling...")
        os.system("c++ test.exe test.cpp") #compile codes

        os.system("echo executing...")
        os.system("./test.exe ")
    else:
        os.system("echo compiling...")
        os.system("make") #compile codes

        os.system("echo executing...")
        os.system("./output.exe")

except:
    os.system("echo compiling...")
    os.system("make") #compile codes

    os.system("echo executing...")
    os.system("./output.exe")
