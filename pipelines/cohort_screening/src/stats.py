import argparse
import json
import pandas as pd


def stats(args):
    with open(args.json_file) as f:
        jdata = json.load(f)

    df = pd.DataFrame(jdata)
    print('Number of variants: ' + str(df.variants.count()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get some stats on this json.')
    parser.add_argument('--json_file', type=str, help='path to json', required=True)

    args = parser.parse_args()

    stats(args)