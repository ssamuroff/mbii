import argparse
import numpy as np
import yaml


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--config', '-c', type=str, action='store')
parser.add_argument('--ctype', '-p', type=str, action='store')
args = parser.parse_args()

options = yaml.load(open(args.config))

if (args.ctype.lower()=='eta'):
	from mbii.pipeline import calculate_eta as fns
	print 'Will calculate %s correlation.'%args.ctype.lower()
if (args.ctype.lower()=='xiee'):
	from mbii.pipeline import calculate_xiee as fns
	print 'Will calculate %s correlation.'%args.ctype.lower()

fns.compute(options)


