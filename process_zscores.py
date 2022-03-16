import logging
import os
import subprocess
import tempfile
from typing import Dict
import fire
import pandas as pd
from Bio import SeqIO

logging.basicConfig(level=logging.DEBUG,)
logger = logging.getLogger('process_zscores')

EXECUTABLE = '/Users/alexandershein/Code/z_coronaviruses/zhunt3-alan'
SEQUENCES_FILE = '/Users/alexandershein/Code/z_coronaviruses/corr_analysis/sequences.fasta'
RESULT_PATH = '/Users/alexandershein/Code/z_coronaviruses/corr_analysis/zscores/'


def compute_and_save_zscores(
    sequences: Dict[str, str], windowsize: int = 8, minsize: int = 6, maxsize: int = 8
):
    processes = []
    temp_filenames = []
    current = 1
    total = len(sequences)
    for accession, sequence in sequences.items():
        fd, temp = tempfile.mkstemp()
        os.close(fd)
        temp_filenames.append(temp)

        with open(temp, 'w') as stream:
            stream.write(sequence)
        logger.info(f'Running zhunt for {current} out of {total} sequences')
        processes.append(
            subprocess.Popen(
                [EXECUTABLE, str(windowsize), str(minsize), str(maxsize), temp],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        )
        current += 1

    logger.info('Waiting for subprocesses to finish')
    for process in processes:
        # Making sure that all background processes finished
        process.wait()

    logger.info('Processing zhunt outputs')
    for i, accession in enumerate(sequences.keys()):
        temp = temp_filenames[i]
        with open(temp + ".Z-SCORE", 'r') as stream:
            df = pd.read_csv(
                stream,
                names=['Start', 'End', 'nu-1', 'nu-2', 'nu-3', 'Z-Score', 'Sequence', 'Conformation'],
                skiprows=1,
                sep='\s+',
            )
        os.remove(temp)
        os.remove(temp + ".Z-SCORE")
        df.to_csv(os.path.join(RESULT_PATH, f'{accession}.csv'))

    logger.info('Done!')


def read_sequences(path: str) -> Dict[str, str]:
    return {record.id: str(record.seq) for record in SeqIO.parse(path, "fasta")}


def main(sequences_file: str):
    sequences = read_sequences(sequences_file)
    compute_and_save_zscores(sequences)


if __name__ == '__main__':
    fire.Fire(main)
