import os

dir = "test" + str(str(1)) + "/"

if not os.path.exists(dir):
    os.mkdir(dir)
    print("Create directory")
else:
    print("Make directory")
