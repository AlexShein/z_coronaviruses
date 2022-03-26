import logging
import multiprocessing as mp
from itertools import chain
from typing import List, Tuple
import os

import fire
import pandas as pd
from Bio import SeqIO

OUT_FILE = f'corr_analysis/data/combined_zscores_{35}_seq.csv'
ZSCORES_PATH = '/Users/alexandershein/Code/z_coronaviruses/corr_analysis/zscores/'
SEQUENCES_FILE = '/Users/alexandershein/Code/z_coronaviruses/corr_analysis/data/aligned_35_sequences.fasta'

logging.basicConfig(
    format='%(asctime)s - %(message)s', level=logging.DEBUG,
)
logger = logging.getLogger('get_coverage')


def generate_coverage(*args) -> List[dict]:
    results = []
    sequence_id, _ = args[0]
    logger.info(f'Processing {sequence_id}')
    logger.info(f'Reading {os.path.join(ZSCORES_PATH, f"{sequence_id}.csv")}')
    zscores_df = pd.read_csv(os.path.join(ZSCORES_PATH, f'{sequence_id}.csv'), header=0,)
    results = zscores_df.apply(
        lambda row: {
            "Start": row['Start'],
            "End": row['End'],
            "Score": row['Z-Score'],
            "Sequence_id": sequence_id,
        },
        axis=1,
    ).tolist()

    logger.info(f'Done processing {sequence_id}')
    return results


def get_sequences_info(path: str) -> List[Tuple[str, int]]:
    res = [(record.id, len(record.seq)) for record in SeqIO.parse(path, "fasta")]
    logger.info(f'Found {len(res)} sequences to process')
    return res


def main():
    sequences_file = SEQUENCES_FILE
    out_file = OUT_FILE
    logger.info('Combininng Zscores into 1 file')
    with mp.Pool(4) as pool:
        results = chain(*pool.map(generate_coverage, get_sequences_info(sequences_file)))
        logger.info('Processing is finished, saving results to file')
    res_df = pd.DataFrame(results)
    res_df.to_csv(out_file)
    logger.info(f'Done!\n{out_file}')


if __name__ == '__main__':
    fire.Fire(main)
