# Measure relative T / C signal at potiental conversion sites
# from aligned Sanger traces.


import pandas as pd
import json
import argparse

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--A', help='Path to Tracy json formated alignment file.')
    parser.add_argument('--S', help='Sample name for aligned file.')
    parser.add_argument('--T', help='Sample treatment')
    parser.add_argument('--R', help='Name of reference sequence.')
    parser.add_argument('--B', help='Boolean. Has the sample been treated with bisulfute?')
    parser.add_argument('--O', help='Output path. Output is always tsv formated.')
    parser.add_argument('--C', help='Topo state of template when treated')
    parser.add_argument('--E', help='Expected read length')

    return parser.parse_args()


def read_tracy_alignment(filepath):

    with open(filepath) as handle:
        return json.load(handle)


def make_TC_table(tracy_dict, sample_name, ref_name, treatment, bisulfite, topo, expected_read_length):

    # get list of all basecall positions in order
    basecall_index = tracy_dict['gappedTrace']['basecallPos']

    TC_table = []

    assert len(basecall_index), len(tracy_dict['altalign'].replace('-', ''))
    print(len(basecall_index), len(tracy_dict['altalign'].replace('-', '')))  # these are the same


    site_counter = 0
    ungapped_steps = 0
    for i in range(len(tracy_dict['altalign'])):
        # check to make sure not a gap in the alt align
        if tracy_dict['altalign'][i] != '-' and i < int(expected_read_length):
            if tracy_dict['refalign'][i].upper() == 'C':
                # check the RFU values at this position
                signal = {}

                for each_base in ['A', 'T', 'G', 'C']:
                    call_index = int(basecall_index[ungapped_steps]) - 1
                    base_rfu = tracy_dict['gappedTrace'][f'peak{each_base}'][call_index]
                    signal[each_base] = base_rfu

                total_signal = sum(signal.values())

                try:
                    T_to_C = signal['T'] / (signal['T'] + signal['C'])
                except ZeroDivisionError:
                    print('Zero divison error', signal)
                    T_to_C = 0
                TC_row = {
                    'total_signal': total_signal,
                    'T_to_C': T_to_C,
                    'T': signal['T'],
                    'C': signal['C'],
                    'refBase': tracy_dict['refalign'][i].upper(),
                    'altBase': tracy_dict['altalign'][i].upper(),
                    'read_index': i,
                    'converted_site_index': site_counter,
                    'sample_name': sample_name,
                    'treatment': treatment,
                    'ref_name': ref_name,
                    'bisulfite': bisulfite,
                    'topo_state': topo
                }
                TC_table.append(TC_row)
                site_counter += 1

            ungapped_steps += 1
    
    # convert list of dicts to a DataFrame to easy export
    return pd.DataFrame(TC_table)


def main():
    
    args = get_args()
    print('Got args')
    tracy_dict = read_tracy_alignment(args.A)
    print('read dict')
    TC_table = make_TC_table(tracy_dict, args.S, args.R, args.T, args.B, args.C, args.E)
    TC_table.to_csv(args.O, sep='\t', index=False)
    print('wrote table')


if __name__ == '__main__':
    main()