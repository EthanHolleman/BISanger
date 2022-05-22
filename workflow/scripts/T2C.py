# Measure relative T / C signal at potiental conversion sites
# from aligned Sanger traces.

import pandas as pd
import json
import argparse

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('A', help='Path to Tracy json formated alignment file.')
    parser.add_argument('S', help='Sample name for aligned file.')
    parser.add_argument('T', help='Sample treatment')
    parser.add_argument('R', help='Name of reference sequence.')
    parser.add_argument('O', help='Output path. Output is always tsv formated.')

    return parser.parse_args()


def read_tracy_alignment(filepath):

    with open(filepath) as handle:
        return json.load(handle)


def make_TC_table(tracy_dict, sample_name, ref_name, treatment):

    # get list of all basecall positions in order
    basecall_index = sorted(tracy_dict['basecalls'].keys(), key=lambda k: int(k))

    TC_table = []

    assert len(tracy_dict['refalign']) == len(tracy_dict['altalign'])
    site_counter = 0
    for i in range(len(tracy_dict['refalign'])):
        # if reference is a C (could have been converted interogate this site)
        if tracy_dict['refalign'][i].upper() == 'C':
            # get signal for each base at this position in the aligned sequence
            signal = {base: float(tracy_dict[f'peak{base}'][i]) for base in ['A', 'T', 'G', 'C']}
            total_signal = sum(signal.values())
            T_to_C = signal['T'] / (signal['T'] + signal['C'])
            TC_row = {
                'total_signal': total_signal,
                'T_to_C': T_to_C,
                'read_index': i,
                'converted_site_index': site_counter,
                'sample_name': sample_name,
                'treatment': treatment,
                'ref_name': ref_name
            }
            TC_table.append(TC_row)
            site_counter += 1
    
    # convert list of dicts to a DataFrame to easy export
    return pd.DataFrame(TC_table)


def main():
    
    args = get_args()
    tracy_dict = read_tracy_alignment(args.A)
    TC_table = make_TC_table(tracy_dict, args.S, args.R, args.T)
    TC_table.to_csv(args.O, sep='\t', index=False)


if __name__ == '__main__':
    main()