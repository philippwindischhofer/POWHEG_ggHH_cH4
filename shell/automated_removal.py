import os
import numpy

#########################################################
##                                                      #
## REMOVAL script for NLO.top files from unstable runs  #
##                                                      #
#########################################################
## To be used if instabilities are present in the pro-  #
## duced histograms pwg-****-NLO.top,                   #
## leading to unphysical spikes in the histograms       #
#########################################################

def listfiles(d):

	top = []
	for f in os.listdir(d):
	  if f.find("pwg-")!=-1 and f.find("-NLO.top")!=-1 and f.find(".outlier") == -1:
	    top += [f]
	return top


path = '.'
hists = ['sigtot-000']  # ['sigtot-000', 'HH-m-000'] or any histogram
                        # you'd like to check for spikes

files = listfiles(path)
data = [[ [] for f in files  ] for h in hists]
baddata = []
biggest = [ [] for h in hists]

print "-----------------------------------------------\n",\
      "               ################                \n", \
      "          ... Starting the search ...          \n\n\n"

for i, h in enumerate(hists):
  for j, f in enumerate(files):
    handle = open(f, 'r')

    tmp = 0.
    flag = False
    for line in handle:
      if flag and line != '\n':
        buf = line
        y = float(buf.split(' ')[3])
        dy = float(buf.split(' ')[4])
        if y > tmp: tmp = y
        data[i][j] += [ [y, dy] ]
      if line.find(h) != -1: flag = True
      if (line.find('#') != -1 or line == '\n') and line.find(h) == -1: flag = False 
    biggest[i] += [tmp]

for i, h in enumerate(hists):
  data[i] = map(list, zip(*data[i]))

#baddata = []

## As a simple trick, check if the relative error in bins with high-statistics
## (i.e. with a total bin weight of at least WMIN * (largest bin weight through all bins))
## is larger than a given threshold ERRMAX. Move the corresponding event pwgevents.lhe to pwgevents.lhe.outlier

## WMIN = 0.2
## ERRMAX = 0.1

WMIN   = 0.2
ERRMAX = 0.1

for i, h in enumerate(hists):
  for j, bin in enumerate(data[i]):
     for k in range(len(bin)):
       if data[i][j][k][1] / data[i][j][k][0] > ERRMAX and data[i][j][k][0] / biggest[i][k] > WMIN:
         baddata += [k]

baddata = list(set(baddata))

print "* * * Found: \033[91m " + str(len(baddata)) + "\033[0m bad files...\n"
for i in range(len(baddata)):
  print "   " + files[baddata[i]]

for i in range(len(baddata)): os.system("mv " + files[baddata[i]] + " " + files[baddata[i]] + ".outlier")
print "\n\n* * *\033[92m CLEAN ! \n \033[0m"

