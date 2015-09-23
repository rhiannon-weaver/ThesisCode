#!/bin/env/python
import os, commands, string

rowids = commands.getoutput("cat datasets/LargeNetworksIndex.txt")
rowids = rowids.splitlines()

pydir = "/Users/rweaver/Research/analysis/Conficker/summaries-and-subsets/"
dir = "/Users/rweaver/Research/analysis/Conficker/spikestatemodel/code/parameterfiles/datasets/SingleNets/"

for rowid in rowids:
	nrows = str(int(rowid)+1)
	fileext = dir + "Network" + rowid
	filename = fileext + "Index.txt"
        s = not os.system("test -f " + filename)
        if not s:
		indfile= file(filename, "w")
		indfile.write(rowid)
		indfile.close()
	s = not os.system("test -f " + fileext +"FL.py")
	if not s:
		os.system("python " + pydir + "getsubsetIPFL.py " + filename + " " + fileext + " " + nrows)
