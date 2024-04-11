from glob import glob
import os, shutil

# clear out old fits
folder = 'newFitData/test_simple/fitdata'

for filename in os.listdir(folder):

    file_path = os.path.join(folder, filename)

    try:
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)

    except Exception as e:
        print('Failed to delete %s. Reason: %s' % (file_path, e))

# delete the executable

exe = 'predictStructureQvary'

if os.path.isfile(exe):
	os.remove(exe)
# print(glob('*'))
