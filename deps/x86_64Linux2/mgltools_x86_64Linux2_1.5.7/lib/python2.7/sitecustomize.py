mglroot = '/disk1/ariane/vscode/enzymetk/software/x86_64Linux2/mgltools_x86_64Linux2_1.5.7'
# specify mglroot here
import sys, os
path = os.path.join(mglroot, "MGLToolsPckgs")
sys.path.append(path)

from os import getenv
if getenv('MGLPYTHONPATH'):
    sys.path.insert(0, getenv('MGLPYTHONPATH'))
    
from Support.path import setSysPath
setSysPath(path)
#sys.path.insert(0,os.path.abspath('.'))
