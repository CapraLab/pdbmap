#!/usr/bin/python27
import sys,gzip

unp_file = sys.argv[1]

with gzip.open(unp_file,'r') as fin:
	data_buffer = []
	confirmed_human=False
	for line in fin:
		row = line.strip().split('   ')
		if row[0] == '//':
			# Found new ID, if the previous entry was human,
			# flush the buffer
			if confirmed_human:
				for dat in data_buffer:
					print dat,
			# Wait for confirmation that next entry is human
			confirmed_human = False
			# Clear the data buffer for the next entry
			data_buffer = []
		elif row[0] == 'OS' and row[1] == 'Homo sapiens (Human).':
			# The current entry is human, flush when finished
			confirmed_human = True
		# Store the row in the data buffer in case it is
		# human and needs to be printed
		data_buffer.append(line)
