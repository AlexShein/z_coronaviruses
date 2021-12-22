import tempfile
import os
import pandas as pd
import subprocess

EXECUTABLE = '/Users/alexandershein/Code/z_coronaviruses/zhunt3-alan'


def alanrun(query: str, windowsize: int = 8, minsize: int = 6, maxsize: int = 8):
    fd, temp = tempfile.mkstemp()
    os.close(fd)

    with open(temp, 'w') as stream:
        stream.write(query)

    subprocess.run([EXECUTABLE, str(windowsize), str(minsize), str(maxsize), temp], check=True)
    with open(temp + ".Z-SCORE", 'r') as stream:
        df = pd.read_csv(
            stream,
            names=['Start', 'End', 'nu-1', 'nu-2', 'nu-3', 'Z-Score', 'Sequence', 'Conformation'],
            skiprows=1,
            sep='\s+',
        )
    os.remove(temp)
    os.remove(temp + ".Z-SCORE")
    return df
