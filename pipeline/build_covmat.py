import argparse
import numpy as np
import yaml
import mbii.symmetrise_lib as lib
from mbii.pipeline.twopoint import calculate_covmat as cov

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--config', '-c', type=str, action='store')
parser.add_argument('--mpi', action='store_true')
args = parser.parse_args()

options = yaml.load(open(args.config))
binning = lib.parse_binning(options)

correlations = options['covariance']['ctypes'].split()
snapshots = options['covariance']['snapshots'].split()
bins = [binning[c] for c in correlations]

cov.compute(options, binning, snapshots, mpi=args.mpi)


