import argparse
import numpy as np
import yaml


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--config', '-c', type=str, action='store')
args = parser.parse_args()

options = yaml.load(open(args.config))

correlations = options['2pt']['ctypes'].split()
mode = options['2pt']['mode'].lower()

for correlation in correlations:
	print('Processing %s'%correlation )
	if (mode=='errors'):
		exec('from mbii.pipeline.twopoint import calculate_%s_errors as fns'%correlation)
	else:
		exec('from mbii.pipeline.twopoint import calculate_%s as fns'%correlation)
	fns.compute(options)


