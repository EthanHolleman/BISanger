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
    parser.add_argument('--Z', help='Template table path')

    return parser.parse_args()


def read_tracy_alignment(filepath):

    with open(filepath) as handle:
        return json.load(handle)
    


def make_template_table(tracy_dict, template_name):
    template_table = []
    for i in range(len(tracy_dict['altalign'])):
        

            row = {
                'template_base': tracy_dict['refalign'][i].upper(),
                'index': i,
                'pos': i+1,
                'template_name': template_name,
                'value': 0
            }

            if tracy_dict['refalign'][i].upper() == 'C':
                row['value'] = 1
            
            template_table.append(row)


    return pd.DataFrame(template_table)



def make_TC_table(tracy_dict, sample_name, ref_name, treatment, bisulfite, topo, expected_read_length):

    TC_table = []

    # get list of all basecall positions in order
    basecalls = tracy_dict['gappedTrace']['basecalls']
    refalign, altalign = tracy_dict['refalign'], tracy_dict['altalign']

    # basecall specify position using base 1 indexing

    # create list of basecall positions sorted by integer value. These are
    # the keys of the basecall dictionary this also appears to be base 1

    basecall_pos_raw = sorted(basecalls.keys(), key=lambda x: int(x))
    for raw_call in basecall_pos_raw:
        call_info = basecalls[str(raw_call)].split(':')
        print(call_info)
        seq_pos = call_info[0]

        total_signal = 0
        signal_dict = {base: 0 for base in ['A', 'T', 'G', 'C']}
        refbase = '-'
        altbase = '-'
        
        if seq_pos != '-':  # is *not* a gap
            seq_pos = int(seq_pos)

            if seq_pos >= int(expected_read_length):
                break
            else:
                for base in signal_dict.keys():
                    signal_dict[base] = tracy_dict['gappedTrace'][f'peak{base}'][int(raw_call)-1]
                    total_signal = sum(signal_dict.values())

                    try:
                        T_to_C = signal_dict['T'] / (signal_dict['C'] + signal_dict['T'])
                    except ZeroDivisionError:
                        T_to_C = 0
                
                refbase = refalign[seq_pos-1]
                altbase = altalign[seq_pos-1]
        
        else:
            seq_pos = 0

        TC_row = {
                    'total_signal': total_signal,
                    'T_to_C': T_to_C,
                    'T': signal_dict['T'],
                    'C': signal_dict['C'],
                    'refBase': refbase,
                    'altBase': altbase,
                    'read_index': seq_pos-1,
                    'sample_name': sample_name,
                    'treatment': treatment,
                    'ref_name': ref_name,
                    'bisulfite': bisulfite,
                    'topo_state': topo
                }
        
        TC_table.append(TC_row)

    return pd.DataFrame(TC_table)










    # TC_table = []

    # assert len(basecall_index), len(tracy_dict['altalign'].replace('-', ''))
    # print(len(basecall_index), len(tracy_dict['altalign'].replace('-', '')))  # these are the same


    # site_counter = 0
    # gaps = 0
    # for i in range(len(tracy_dict['altalign'])):
    #     # check to make sure not a gap in the alt align
    #     if tracy_dict['altalign'][i] == '-':
    #         gaps += 1
    #     elif i < int(expected_read_length):
    #         if tracy_dict['refalign'][i].upper() == 'C':
    #             # check the RFU values at this position
    #             signal = {}

    #             for each_base in ['A', 'T', 'G', 'C']:
    #                 call_index = int(basecall_index[i]) - 1
    #                 base_rfu = tracy_dict['gappedTrace'][f'peak{each_base}'][call_index]
    #                 signal[each_base] = base_rfu

    #             total_signal = sum(signal.values())

    #             try:
    #                 T_to_C = signal['T'] / (signal['T'] + signal['C'])
    #             except ZeroDivisionError:
    #                 print('Zero divison error', signal)
    #                 T_to_C = 0
    #             TC_row = {
    #                 'total_signal': total_signal,
    #                 'T_to_C': T_to_C,
    #                 'T': signal['T'],
    #                 'C': signal['C'],
    #                 'refBase': tracy_dict['refalign'][i].upper(),
    #                 'altBase': tracy_dict['altalign'][i].upper(),
    #                 'read_index': i,
    #                 'converted_site_index': site_counter,
    #                 'sample_name': sample_name,
    #                 'treatment': treatment,
    #                 'ref_name': ref_name,
    #                 'bisulfite': bisulfite,
    #                 'topo_state': topo
    #             }
    #             TC_table.append(TC_row)
    #             site_counter += 1
    

def main():
    
    args = get_args()
    print('Got args')
    tracy_dict = read_tracy_alignment(args.A)
    print('read dict')
    TC_table = make_TC_table(tracy_dict, args.S, args.R, args.T, args.B, args.C, args.E)
    template_table = make_template_table(tracy_dict, args.R)
    TC_table.to_csv(args.O, sep='\t', index=False)
    template_table.to_csv(args.Z, sep='\t', index=False)
    print('wrote table')


if __name__ == '__main__':
    main()