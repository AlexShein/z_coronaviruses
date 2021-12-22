from typing import List, Tuple
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
import logging

ZSCORE_THRESHOLDS = [50, 150, 300]  # Arbitrary
SEQUENCES_FILE = '/Users/alexandershein/Code/z_coronaviruses/evolutionary_zscores/61_coronaviruses.fasta'

logging.basicConfig(level=logging.DEBUG,)
logger = logging.getLogger('get_coverage')


def generate_coverage(*args):
    sequence_id, sequence_length = args
    pass


def get_sequences_stats(path: str) -> List[Tuple[str, int]]:
    return [(record.id, len(record.seq)) for record in SeqIO.parse(path, "fasta")]


if __name__ == '__main__':
    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(generate_coverage, get_sequences_stats(SEQUENCES_FILE))

