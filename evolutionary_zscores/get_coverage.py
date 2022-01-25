import logging
import multiprocessing as mp
from itertools import chain
from typing import List, Tuple

import fire
import pandas as pd
import pyranges as pr
from Bio import SeqIO

ZSCORE_THRESHOLDS = [50, 150, 300]  # Arbitrary
SEQUENCES_FILE = '/Users/alexandershein/Code/z_coronaviruses/evolutionary_zscores/61_coronaviruses.fasta'
STEP_COVERAGE = 10
OUT_FILE_POSITIONS = 'zscores_by_position.csv'
OUT_FILE_MERGED = 'zscores_merged.csv'

logging.basicConfig(
    format='%(asctime)s - %(message)s', level=logging.DEBUG,
)
logger = logging.getLogger('get_coverage')


def generate_coverage(*args) -> List[dict]:
    results = []
    sequence_id, sequence_length = args[0]
    logger.info(f'Processing {sequence_id}')
    df = pd.read_csv(f'zscores/{sequence_id}.csv', header=0,).iloc[:, 1:]
    df['Chromosome'] = 'chr1'
    for threshold in ZSCORE_THRESHOLDS:
        ranges = pr.PyRanges(
            df[df['Z-Score'] >= threshold], chromosomes='Chromosome', starts='Start', ends='End'
        )
        for coord in range(0, sequence_length, STEP_COVERAGE):
            intersections = ranges.intersect(
                pr.from_dict({"Chromosome": ['chr1'], "Start": [coord], "End": [coord + STEP_COVERAGE]})
            )
            intersections_df = intersections.as_df()
            mean_score = 0 if intersections_df.empty else intersections_df['Z-Score'].mean()
            median_score = 0 if intersections_df.empty else intersections_df['Z-Score'].median()
            results.append(
                {
                    "Start": coord,
                    "End": coord + STEP_COVERAGE,
                    "Start_r": coord / sequence_length,
                    "End_r": (coord + STEP_COVERAGE) / sequence_length,
                    'Mean_Z-Score': mean_score,
                    'Median_Z-Score': median_score,
                    "Threshold": threshold,
                    "Sequence_id": sequence_id,
                }
            )
    logger.info(f'Done processing {sequence_id}')
    return results


def merge_intersecting_ranges(*args) -> List[dict]:
    results = []
    sequence_id, _ = args[0]
    logger.info(f'Processing {sequence_id}')
    df = pd.read_csv(f'zscores/{sequence_id}.csv', header=0,).iloc[:, 1:]
    df['Chromosome'] = 'chr1'
    for threshold in ZSCORE_THRESHOLDS:
        ranges_df = (
            pr.PyRanges(df[df['Z-Score'] >= threshold], chromosomes='Chromosome', starts='Start', ends='End')
            .merge()
            .as_df()
        )
        ranges_df['Threshold'] = threshold
        ranges_df['Sequence_id'] = sequence_id

        results.extend(ranges_df.to_dict(orient='records'))

    logger.info(f'Done processing {sequence_id}')
    return results


def get_sequences_info(path: str) -> List[Tuple[str, int]]:
    res = [(record.id, len(record.seq)) for record in SeqIO.parse(path, "fasta")]
    logger.info(f'Found {len(res)} sequences to process')
    return res


def main(sequences_file: str, out_file: str, merge: bool = False):
    func = generate_coverage if not merge else merge_intersecting_ranges
    if merge:
        logger.info('Running merge')
    else:
        logger.info('Running coverage')
    with mp.Pool(4) as pool:
        results = chain(*pool.map(func, get_sequences_info(sequences_file)))
        logger.info('Processing is finished, saving results to file')
        res_df = pd.DataFrame(results)
        res_df.to_csv(out_file)
        logger.info(f'Done!\n{out_file}')


if __name__ == '__main__':
    fire.Fire(main)
