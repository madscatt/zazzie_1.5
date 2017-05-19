import os,glob,locale,re,shutil
import platform
import datetime

def back_up_existing_module_folder(runname,module):
    folders = glob.glob(os.path.join(runname,module+'*'))
    if os.path.join(runname,module) not in folders:
        return
    max_num = 0
    for folder in folders:
        folder_name = os.path.basename(folder)
        if re.compile('^'+module+'_\d+$').match(folder_name):
            num = locale.atoi(folder_name[len(module)+1:])
            if num>max_num:
                max_num = num
    shutil.move(os.path.join(runname,module), os.path.join(runname,module+'_%d'%(max_num+1)))

    return


import os

def creation_date(path_to_file):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    if platform.system() == 'Windows':
        return os.path.getctime(path_to_file)
    else:
        stat = os.stat(path_to_file)
        try:
            return stat.st_birthtime
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return stat.st_mtime

def modification_date(filename):
    t = os.path.getmtime(filename)
    return datetime.datetime.fromtimestamp(t).strftime("%Y_%b_%d_%H:%M:%S")

if __name__ == '__main__':

    path = './'
    runname = 'run_fred'
    module = 'happy_feet'
    back_up_existing_module_folder(runname,module)


