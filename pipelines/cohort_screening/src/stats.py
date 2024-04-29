# useful for getting some stats on the positions file after converting to df


import argparse
import json
import pandas as pd


def stats(args):
    with open(args.json_file, "r") as f:
        jdata = json.load(f)
        fstring = str(jdata)
        print('Number of transcripts: ' + str(fstring.count('transcripts')))

    df = pd.DataFrame(jdata)
    print('Number of variants: ' + str(df.variants.count()))
    print('other stats:')
    print(df.describe(include='all'))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get some stats on this json.')
    parser.add_argument('--json_file', type=str, help='path to json', required=True)

    args = parser.parse_args()

    stats(args)