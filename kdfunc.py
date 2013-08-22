def pwmdicter(pwms):
	"""Return of dict of pwms and ChIP experiments mapped to TFs."""
	pwmlist = {}
	for line in pwms:
		if '#' in line:
			continue
		if line.strip().split()[0] in pwmlist.keys():
			pwmlist[line.strip().split()[0]].append(line.strip().split()[1])
		if line.strip().split()[0] not in pwmlist.keys():
			pwmlist[line.strip().split()[0]] = [line.strip().split()[1]]
	return(pwmlist)

def genereader(thisfile):
	"""Return array of genes from file split by '|'."""
	currfile = open(thisfile,'r')
	return(currfile.readline().strip().split('|'))

def pathmaker(path):
	"""Make a new directory unless it already exists."""
	import os
	try:
		os.makedirs(path)
	except OSError:
		print 'Warning: ' + path + ' already exists!'

def bedmaker(bedrecord,tssrecord):
	"""Returns bed-like record of all binding within 'w' kb of TSSs."""
	sortrecord = bedrecord.sort()
	intersectbed = tssrecord.window(sortrecord, w=windowsize)
	setdiffbed = tssrecord.window(sortrecord, w=windowsize, v=True)
	reportbed = intersectbed.cat(setdiffbed,postmerge=False,force_truncate=False)
	return(reportbed)

if __name__ == '__main__':
	func()