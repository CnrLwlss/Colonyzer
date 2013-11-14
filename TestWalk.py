import os
wlk=os.walk(os.getcwd())

for dirname, dirnames, filenames in os.walk(os.getcwd()):
    if not any([x in dirname for x in ['Output_Images','.git']]):
        for filename in filenames:
            print os.path.join(dirname, filename)
