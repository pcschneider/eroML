#!/bin/bash

hn=$(hostname)
echo "Hostname=${hn}"


if [ $hn == 'thinkpad-T490s' ]; then
        echo "Exporting 'PYTHONPATH=/home/majestix/hdd/Xdata/ero/'."
        export PYTHONPATH=/home/majestix/hdd/Xdata/ero/
#         exit 0;
elif [ $hn == 'pc13' ]; then
        echo "Exporting 'PYTHONPATH=/hs/pc13/data/stgh325/Xdata/ero'."
        export PYTHONPATH=/hs/pc13/data/stgh325/Xdata/ero
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hs/pc13/data/stgh325/python/lib/

else 
    echo "I don't know the path for hostname=${hn}"
fi


# export PYTHONPATH=/home/majestix/hdd/Xdata/ero/
# export PYTHONPATH=/hs/pc13/data/stgh325/Xdata/ero
