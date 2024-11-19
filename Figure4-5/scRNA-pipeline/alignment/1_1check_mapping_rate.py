#! /usr/bin/python
""" check overall statistics for all log files from star aligner"""

from os import listdir
from os.path import isfile
import re

OUTPUT_PATH_STAR = '/data/nier/cycle/RNA-seq/singlecell/mouse/esc_npc/mesc/'

if __name__ == "__main__":
	regex = "(?<=Uniquely mapped reads \% \|\t)[0-9]+\.[0-9]+"
	regex = re.compile(regex)
	
	stats = {}
	for folder in listdir(OUTPUT_PATH_STAR):
		log_file = "{0}/{1}/Log.final.out"\
			   .format(OUTPUT_PATH_STAR, folder)
	
		if not isfile(log_file):
			continue

		stats[folder] = regex.findall(
			open(log_file, 'r').read())[0]
	
	f_csv = open(OUTPUT_PATH_STAR + 'aligner_unique_read.csv', 'w')

	for key in stats:
		f_csv.write('{0},{1}\n'.format(key,stats[key]))


