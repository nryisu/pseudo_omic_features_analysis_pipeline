#!/usr/bin/python
import sys
from fnmatch import fnmatch

from os import popen
from os import listdir
from os import mkdir
from os.path import isdir
from os.path import isfile 
from os.path import getsize

#from sys import agrv
from sys import stdout

from distutils.dir_util import mkpath


OUTPUT_PATH = '/data/nier/cycle/RNA-seq/singlecell/mouse/esc_npc/mesc/'
FASTQ_PATH = '/data/nier/cycle/RNA-seq/singlecell/mouse/esc_npc/mesc/'
STAR_INDEX_PATH = '/data/reference/mouse/mm10_STAR/RNA'
PATH_SOFTWARE = '/data/biosoft/software/STAR-2.7.1a/bin/Linux_x86_64/STAR'
PATTERN = ""
GTF = "/data/reference/mouse/mm10_STAR/RNA/Mus_musculus.GRCm38.93.gtf"
#STAR_INDEX_READ_LENGTH = 76
#THREADS = 1
#sleep(2 * random()

if __name__ == "__main__":
	for fil in listdir(FASTQ_PATH):
		if isfile(FASTQ_PATH + fil):
			continue
		if PATTERN and not fnmatch(fil, PATTERN):
			continue

		#print "====> file to be aligned:", fil

		if not isdir(OUTPUT_PATH + fil):
			mkdir(OUTPUT_PATH + fil)

		#if isfile(OUTPUT_PATH + fil + "Aligned.sortedByCoord.out.bam") \
		#	and getsize(OUTPUT_PATH + fil + "Aligned.sortedByCoored.out.bam"):
		#		print 'bam file result already exists for:{0}\nskoipping...' \
		#		.format(fil)
		#		continue

		fastq_str_1 = ""
		fastq_str_2 = ""
		for fastq_fil in listdir(FASTQ_PATH + fil):
		#	print fastq_fil
			if fnmatch(fastq_fil,"*_1.fastq.gz"):
				fastq_str_1 += "{0}{1}/{2}".format(FASTQ_PATH, fil, fastq_fil)
			if fnmatch(fastq_fil,"*_2.fastq.gz"):
				fastq_str_2 += "{0}{1}/{2}".format(FASTQ_PATH, fil, fastq_fil)
				#print(fastq_str)
			#if not fastq_str:
			#	print 'no fastq file found for:{0}!\n'.format(fil)
				#print(OUTPUT_PATH + fil + "/")
		cmd = "{0} --readFilesIn {1} {2} --readFilesCommand zcat --runThreadN 4 " \
		      " --outFileNamePrefix {3} --genomeDir {4} --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate" \
		      .format(PATH_SOFTWARE,
			      fastq_str_1,
			      fastq_str_2,
			      OUTPUT_PATH + fil + "/",
			      STAR_INDEX_PATH
		      )
		#print(cmd)
		res = popen(cmd)
		c = res.read(1)
		
		while c:
			stdout.write(c)
			stdout.flush()
			c = res.read(1)
