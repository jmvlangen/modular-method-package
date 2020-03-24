load('~/Documents/SageFiles/load.sage')

import traceback

try:
    load('<path>')
except Exception as e:
    print("")
    print(">>>> Exception <<<<<")
    print("")
    print(traceback.format_exc())
    
    
