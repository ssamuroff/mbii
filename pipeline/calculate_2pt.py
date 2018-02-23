import argparse
import numpy as np
import yaml


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--config', '-c', type=str, action='store')
args = parser.parse_args()

options = yaml.load(open(args.config))

correlations = options['2pt']['ctypes'].split()

for correlation in correlations:
	print 'Processing %s'%correlation 
	exec('from mbii.pipeline.twopoint import calculate_%s as fns'%correlation)
	fns.compute(options)


