#import group specific data from "inputfile_group_10.txt"

file=open("inputfile_group_10.txt", 'r')
lines=file.readlines()
interval = False

for line in lines:
	if len(line.strip())>0:
		line=line.split()
		if line[0]=='L':
			L=float(line[2])
		if line[0]=='hz':
			hz=float(line[2])
		if line[0]=='k':
			k=float(line[2])
		if line[0]=='c':
			c=float(line[2])
		if interval:
			if line[0].strip() == ']':
				interval = False
				elements_modify = []
				for interval in elements_intervals:
					#print(interval.strip(","))
					elements_modify += range(int(interval.split("-")[0]), int(interval.split("-")[1].strip(",")) + 1)
			else:
				elements_intervals += line
		if line[0]=='elements_to_be_modified':
			interval = True
			elements_intervals = []
		if line[0]=='q(y=0)':
			bc_bottom=float(line[2])
		if line[0]=='T(y=L)':
			bc_top=float(line[2])

