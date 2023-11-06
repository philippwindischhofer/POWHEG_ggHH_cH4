#! /usr/bin/env python
##
##  Author(s):
##   Simone Alioli, Gionata Luisoni
##

import math
import copy
import os,sys,logging,time
import re
import types
import numpy as np
import cStringIO
import collections as cl
from collections import Iterable

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib import ticker as mticker
from matplotlib import rc




class TDHist:
  """ A pair of dictionaries containing the histograms """

  def fromfile(self,f):
    title=''
    while True:
      line=f.readline()
      if line=='':
        return
      if re.search('#',line)!=None:
        title=re.sub('#',' ',line).strip()
        title=re.sub('  *','_',title)
        title=title.split()[0]
        sio=cStringIO.StringIO()
        while True:
          last_pos = f.tell()
          line=f.readline().replace('D','E')
          if re.search("#", line)!=None or not line:
            f.seek(last_pos)
            break
          if len(line)>2:
            sio.write(line)
          if sio.tell()>0:
            sio.seek(0)
            self.hdata[title]=np.loadtxt(sio,comments='(',ndmin=2)
            self.rdata[title]=np.transpose(self.hdata[title])
        sio.close()

  def tofile(self,f):
    for title,histo in self.hdata.items():
      title=re.sub('_',' ',title)
      f.write('# {0}\n'.format(title))
      for row in histo:
        f.write('{: e}  {: e}  {: e}  {: e}\n'.format(*row).replace('e','D'))
      f.write('\n\n')


  def __init__(self,f=None):
    self.hdata=cl.OrderedDict()
    self.rdata=cl.OrderedDict()
    if f:
      self.fromfile(f)



def genplots(outdir,prefix,files,rslist):

  if len(files) != len(rslist):
    logging.error("Something is wring with the number of analyzed files.")
    sys.exit(2)

  # Open output file:
  genplotsfile = open(outdir+'/genplots.gp', 'w')

  for tit,hist in rslist[0].hdata.items():

    title=tit.split('index')[0][:-1]
    index=tit.split('index')[1]
    index=re.sub('_','',index)
    filename=prefix+'-'+title+'_index_'+str(index)
    genplotsfile.write('\n')
    genplotsfile.write('reset\n')
    genplotsfile.write('set term post eps enhanced color dashed "Helvetica" 24\n')
    genplotsfile.write('set output "%s.eps"\n'%filename)
    genplotsfile.write('iindex = "{}"\n'.format(title))

    for f in enumerate(files):
      plottitle=f[1].split('/')[-1]
      genplotsfile.write('title{0}  = "{1}"\n'.format(str(f[0]+1),plottitle))
      genplotsfile.write('infile{0} = "{1}"\n'.format(str(f[0]+1),f[1]))

    for i in range(len(files)):
      genplotsfile.write('allfiles{0} = "< pastegnudata -i ".iindex." ".infile1." ".infile{1}." "\n'.format(str(i+1),str(i+1)))

    genplotsfile.write('\n')
    genplotsfile.write('set multiplot')
    genplotsfile.write('\n\n')
    genplotsfile.write('set origin 0.0,0.3\n')
    genplotsfile.write('set size 1.0,0.7\n')
    genplotsfile.write('set lmargin at screen 0.2\n')
    genplotsfile.write('set rmargin at screen 0.95\n')
    genplotsfile.write('set tmargin at screen 0.95\n')
    genplotsfile.write('set bmargin at screen 0.4\n')
    genplotsfile.write('\n')
    genplotsfile.write('# set log x\n')
    genplotsfile.write('set log y\n')
    genplotsfile.write('\n')
    genplotsfile.write('set datafile fortran\n')
    genplotsfile.write('set style data xyerrorlines\n')
    genplotsfile.write('\n')
    genplotsfile.write('set key font ",15"\n')
    genplotsfile.write('set xtics format ""\n')
    genplotsfile.write('set format y "10^{%L}"\n')
    genplotsfile.write('set xrange [] writeback\n')
    genplotsfile.write('\n')
    genplotsfile.write('plot \\')
    genplotsfile.write('\n')
    for f in range(len(files)-1):
      genplotsfile.write('infile{0} index iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title title{0},\\'.format(str(f+1)))
      genplotsfile.write('\n')
    genplotsfile.write('infile{0} index iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title title{0}\n'.format(str(len(files))))

    genplotsfile.write('set tmargin at screen 0.4\n')
    genplotsfile.write('set bmargin at screen 0.1\n')
    genplotsfile.write('\n')
    genplotsfile.write('unset log y\n')
    genplotsfile.write('set yrange [0.5:1.5]\n')
    genplotsfile.write('set xrange restore\n')
    genplotsfile.write('set format x\n')
    genplotsfile.write('set format y\n')
    genplotsfile.write('set key off\n')
    genplotsfile.write('\n')

    if len(files) > 1:
      genplotsfile.write('plot \\')
      genplotsfile.write('\n')
      genplotsfile.write('infile1 index iindex using ($1+$2)/2:($3/$3):1:2:(1-($4/$3)*2**0.5):(1+($4/$3)*2**0.5) title title1,\\')
      genplotsfile.write('\n')
      for i in range(1,len(files)):
        genplotsfile.write('allfiles{0} using ($1+$2)/2:($7/$3):1:2:($7/$3)*(1-(($4/$3)**2+($8/$7)**2)**0.5):($7/$3)*(1+(($4/$3)**2+($8/$7)**2)**0.5) title title{1},\\'.format(str(i+1),str(i+1)))
        genplotsfile.write('\n')
      genplotsfile.write('1 with lines linestyle -1\n')

    else:
      genplotsfile.write('allfiles using ($1+$2)/2:($3/$3):1:2:(1-($4/$3)*2**0.5):(1+($4/$3)*2**0.5) title title1,\\')
      genplotsfile.write('\n')
      genplotsfile.write('plot 1 with lines linestyle -1\n')

    genplotsfile.write('unset multiplot\n')

  genplotsfile.close()

  # Call gnuplot
  logging.info('Generating plots..')
  pwd = os.getcwd()
  os.chdir(outdir)
  os.system('gnuplotsplit.sh genplots.gp')
  os.chdir(pwd)




####################
if __name__ == '__main__':

  ## Find number of cpu cores
  import multiprocessing
  try:
    numcores = multiprocessing.cpu_count()
  except NotImplementedError:
    numcores = 1

  ## Argument parser
  from argparse import ArgumentParser

  parser = ArgumentParser(description='Manipulates POWHEG gnuplot output files.')
  parser.add_argument('inputfiles', nargs='+', metavar='pwgNLO.top', help='One or more analyzer input files. When more than one input file is given, a comparison plot between them is created.')
  parser.add_argument("-o", "--output", dest="OUTPUT", default="plots", help="name of directory for output plots [plots]")
  parser.add_argument("-f", "--format", dest="FORMAT", default="genplots", help="output format for plots [genplots]")
  parser.add_argument("-s","--show", action='store_true', dest="SHOW", default=False, help="show results in default web browser [NO]")
  parser.add_argument("-w","--overwrite", action='store_true', dest="OVERWRITE", default=False, help="overwrite existing plots [NO]")
  parser.add_argument("-n", "-j", "--num-threads",action="store", dest='NUM_THREADS', type=int, default=numcores, help="max number of threads to be used [%s]" % numcores)
  parser.add_argument("-a", "--alpha", action="store", dest='ALPHA_BAND', type=float, default=0.2, help="transparency of error bands [0.2]")
  parser.add_argument("-v", "--verbose", action="store_const", const=logging.DEBUG, dest="LOGLEVEL", default=logging.INFO, help="print debug (very verbose) messages [NO]")
  parser.add_argument("-q", "--quiet", action="store_const", const=logging.WARNING, dest="LOGLEVEL", default=logging.INFO, help="be very quiet [NO]")
  parser.add_argument("-p", "--prefix", dest="PREFIX", default="xxx", help="prefix to histogram file names")

  args = parser.parse_args()
  logging.basicConfig(level=args.LOGLEVEL, format="%(message)s")
  filenames = args.inputfiles
  outdir    = args.OUTPUT
  logging.debug('Saving output to directory %s'%outdir)
  refFile = filenames[0]

  logging.info('Loading %s file(s):',len(filenames))
  logging.info('%s',filenames)

  # Output format:
  # default: genplots style
  outform = args.FORMAT
  prefix  = args.PREFIX

  # define colors
  if (len(filenames) < 8):
     color_cycle = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
  else:
     colormap = plt.cm.rainbow
     color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(filenames))]

  if os.path.exists(outdir) and os.path.isdir(outdir):
    if args.OVERWRITE==False:
      logging.error("Output directory "+outdir+" already existing.\n")
      logging.error("Call again with options -w or --overwrite to overwrite it.\n")
      sys.exit(2)
    else:
      import shutil
      shutil.rmtree(outdir)
      os.makedirs(outdir)
  else:
    if os.path.exists(outdir):
      logging.error(outdir+" exists but is not a valid directory.\n")
      logging.error("Please check your inputs.\n")
      sys.exit(2)
    else:
      os.makedirs(outdir)
  os.system("chmod 755 "+ outdir)

  rs = []
  analyzerfiles= []

  # Read in files
  for c,rfile in enumerate(filenames):

    response=rfile

    analyzerfiles.append(os.path.abspath(rfile))
    logging.info('Processing %s',rfile)
    file=open(rfile,'r')
    rs.append(TDHist(file))
    file.close()


  # Switch for different output formats:

  if outform == "genplots":
    genplots(outdir,prefix,analyzerfiles,rs)

  elif outform == "write":

    for f in enumerate(analyzerfiles):
      filename='new'+f[1].split('/')[-1]
      file=open(outdir+'/'+filename,'w')
      rs[f[0]].tofile(file)
      file.close()

  elif outform == "rebcum":

    # rebin={'b1-pt-001-050_index_176':4,
    #        'b1-pt-100-050_index_230':4,
    #        'b1-pt-500-050_index_284':4,
    #        'bb-pt-001-050_index_214':4,
    #        'bb-pt-100-050_index_268':4,
    #        'bb-pt-500-050_index_322':4,
    #        'W-pt-100-050_index_223':4,
    #        'W-pt-500-050_index_277':4,
    #        'j1-pt-100-050_index_242':4,
    #        'j1-pt-500-050_index_296':4,
    #        'lept-pt-001-050_index_172':4,
    #        'lept-pt-100-050_index_226':4,
    #        'lept-pt-500-050_index_280':4,
    #        'm-bb-pt-001-050_index_479':4,
    #        'm-bb-pt-100-050_index_518':4,
    #        'm-bb-pt-500-050_index_557':4,
    #        'm-Wbb-m-001-050_index_472':4,
    #        'm-Wbb-m-100-050_index_511':4,
    #        'm-Wbb-m-500-050_index_550':4,
    #        'm-Wbb-ptratio-001-050_index_468':2,
    #        'm-Wbb-ptratio-100-050_index_507':2,
    #        'm-Wbb-ptratio-500-050_index_546':2,
    #        'm-W-pt-100-050_index_487':4,
    #        'm-W-pt-500-050_index_526':4,
    #        'm-j1-pt-100-050_index_494':4,
    #        'm-j1-pt-500-050_index_533':4,
    #        'm-W-bb-dphi-001-050_index_475':2,
    #        'm-W-bb-dphi-100-050_index_514':2,
    #        'm-W-bb-dphi-500-050_index_553':2,
    #        'm-W-bb-dr-001-050_index_476':2,
    #        'm-W-bb-dr-100-050_index_515':2,
    #        'm-W-bb-dr-500-050_index_554':2,
    #        'st-j1-pt-100-000_index_613':2,
    #        'st-j1-pt-100-500_index_706':2,
    #        'st-j1-pt-100-1e3_index_799':2,
    #        'st-j1-y-100-500_index_704':2,
    #        'st-j1-y-100-1e3_index_797':2,
    #        'st-j1-y-500-500_index_735':2,
    #        'st-j1-y-500-1e3_index_828':2,
    #        'st-Wb-pt-100-000_index_594':4,
    #        'st-Wb-pt-500-000_index_625':4,
    #        'st-Wb-y-100-500_index_685':2,
    #        'st-Wb-y-100-1e3_index_778':2,
    #        'st-Wb-y-500-500_index_716':2,
    #        'st-Wb-y-500-1e3_index_809':2,
    #        'W-bb-dphi-001-050_index_208':2,
    #        'W-bb-dphi-100-050_index_262':2,
    #        'W-bb-dphi-500-050_index_316':2,
    #        'W-bb-dr-001-050_index_209':2,
    #        'W-bb-dr-100-050_index_263':2,
    #        'W-bb-dr-500-050_index_317':2,
    #        'm-Wbb-ptratio-500-000_index_429':2,
    #        'm-Wbb-ptratio-500-050_index_546':2}

    # cumulative={'rebin4-m-bb-pt-001-050_index_4791':'m-sigtot_Wbb-001-050_index_444',
    #             'rebin4-m-bb-pt-100-050_index_5181':'m-sigtot_Wbb-100-050_index_483',
    #             'rebin4-m-bb-pt-500-050_index_5571':'m-sigtot_Wbb-500-050_index_522',
    #             'rebin4-m-Wbb-m-001-050_index_4721':'m-sigtot_Wbb-001-050_index_444',
    #             'rebin4-m-Wbb-m-100-050_index_5111':'m-sigtot_Wbb-100-050_index_483',
    #             'rebin4-m-Wbb-m-500-050_index_5501':'m-sigtot_Wbb-500-050_index_522'}

    lumi = 20.0e6

    ##################################################################################

    if len(rebin) > 0:
      for f in rs:
        for h in rebin.items():
          # Add here check on the existence of the histogram!
          m     = f.hdata[h[0]].shape[0]
          step  = h[1]
          title = 'rebin'+str(h[1])+'-'+h[0].split('index')[0]
          index = h[0].split('index')[1]
          title = title+'index'+str(index)+'1'
          sio   = cStringIO.StringIO()
          #print f.hdata[h[0]]
          for m in range(step-1,m,step):
            area = 0
            err2 = 0
            for i in range(step):
              area = area + f.hdata[h[0]][m-(step-1)+i][2]*(f.hdata[h[0]][m-(step-1)+i][1]-f.hdata[h[0]][m-(step-1)+i][0])
              err2 = err2 + f.hdata[h[0]][m-(step-1)+i][3]**2*(f.hdata[h[0]][m-(step-1)+i][1]-f.hdata[h[0]][m-(step-1)+i][0])
            area = area/(f.hdata[h[0]][m][1]-f.hdata[h[0]][m-(step-1)][0])
            err2 = err2/(f.hdata[h[0]][m][1]-f.hdata[h[0]][m-(step-1)][0])
            line=str(f.hdata[h[0]][m-(step-1)][0])+' '+str(f.hdata[h[0]][m][1])+' '+str(area)+' '+str(math.sqrt(err2))+'\n'
            sio.write(line)
          sio.seek(0)
          f.hdata[title]=np.loadtxt(sio,comments='(',ndmin=2)
          f.rdata[title]=np.transpose(f.hdata[title])

  #   for f in enumerate(analyzerfiles):
  #     filename='new-'+f[1].split('/')[-1]
  #     file=open(outdir+'/'+filename,'w')
  #     rs[f[0]].tofile(file)
  #     file.close()

  # elif outform == "cumulative":

  # print rs[0].hdata.keys()

    if len(cumulative) > 0:
      for f in rs:
        for h in cumulative.items():
          mmax     = f.hdata[h[0]].shape[0]
          totxs = f.hdata[h[1]][0][2]
          errxs = f.hdata[h[1]][0][3]
          title = 'cumul-'+h[0].split('index')[0]
          index = h[0].split('index')[1]
          title = title+'index'+str(index)+'2'
          sio   = cStringIO.StringIO()
          for m in range(len(f.hdata[h[0]])):
            area = totxs
            err2 = 0
            for j in range(m):
              area = area - f.hdata[h[0]][j][2]*(f.hdata[h[0]][j][1]-f.hdata[h[0]][j][0])
            for j in range(m,mmax):
              err2 = err2 + f.hdata[h[0]][j][3]**2

            line=str(f.hdata[h[0]][m][0])+' '+str(f.hdata[h[0]][m][1])+' '+str(area*lumi)+' '+str(math.sqrt(err2)*lumi/(mmax-m+1))+'\n'
            sio.write(line)
          sio.seek(0)
          f.hdata[title]=np.loadtxt(sio,comments='(',ndmin=2)
          f.rdata[title]=np.transpose(f.hdata[title])

    for f in enumerate(analyzerfiles):
      filename='new-'+f[1].split('/')[-1]
      file=open(outdir+'/'+filename,'w')
      rs[f[0]].tofile(file)
      file.close()



#  for tit,hist in rs[0].hdata.items():
#    print tit
#    print hist



# import sys, os
# import math
# from optparse import OptionParser
# from os.path  import basename, splitext


# # parser of command line inputs:
# parser = OptionParser()

# parser.add_option("-i", "--input", dest="infilename",
#                   help="name of input file", metavar="INFILE")
# parser.add_option("-o", "--output", dest="outfilename",
#                   help="name of output file", metavar="OUTFILE")
# parser.add_option("-v", "--verbosity", dest="verbosity",
#                   help="sets verbosity of output", metavar="VERBOS")
# parser.add_option("-p", "--process", dest="process",
#                   help="name of the process in histograms", metavar="PROC")


# # verbosity ...
# (options, args) = parser.parse_args()

# if options.verbosity is not None:
#     verb = int(options.verbosity)
# else:
#     verb = 0

#     if options.process is not None:
#         process = options.process
#     else:
#         process = ''

# #####################################################

# filepath = os.path.realpath(options.infilename)

# print "----------------------------------"
# print "Input file is: %s" % filepath
# print "Output file is: %s" % options.outfilename
# print "----------------------------------"

# inputfile  = open(options.infilename,  'r')
# outputfile = open(options.outfilename, 'w')

# lin_n = 0
# his_n = 0
# x_min = 0
# x_max = 0

# his_n = 0

# skip_observables = []

# cumulated_observables = {}
# #    'rebin4_m-bb-pt-001-050':'sigma_rebin4_m-bb-pt-001-050'
# #}

# rebin_observables = {
# #    'bb-pt-001-050':4
# #    'bb-pt-050-050':4
#     'm-j1-pt-050-050':4
# }

# total_xs_observable="sigtot Wbb-001-050"

# lines    = inputfile.readlines()
# prevline = ""
# name     = ""
# skip     = False
# readxs   = False

# cumul_store = False
# cumul_name  = ""
# cumul_dict  = {}
# cumul_histo = {}
# cumul_i     = 0
# cumul_y     = 0
# cumul_dy2   = 0

# rebin = False
# rebin_size = 1
# rebin_dict  = {}
# rebin_histo = {}
# rebin_i     = 0
# rebin_y     = 0
# rebin_dy2   = 0
# newbin = {}


# for line in lines:

#     if line.startswith("#"):
#         #print "Title line found: " + line
#         his_n += 1

#         #--> Observable:
#         histo = line.rsplit('#', 1)[1].rsplit('index',1)[0].strip()
#         if histo == total_xs_observable:
#             readxs = True

#         #--> Observables to skip:
#         if histo in skip_observables:
#             print "Skip    : " + histo
#             skip = True

#         if histo in cumulated_observables:
#             cumul_name  = cumulated_observables[histo]
#             print "Cumulate: " + histo + " --> " + cumul_name
#             cumul_store = True
#             cumul_i = 0

#         #--> Observables to rebin:
#         if histo in rebin_observables:
#             rebin_size = rebin_observables[histo]
#             print "Rebin   : " + histo + " --> rebin" + str(rebin_size) + "_" + histo
#             rebin = True


#     #--> Write copy:
#     if skip is False:
#         outputfile.write(line)

#     if len(line.strip()) == 0 and skip is True:
#         skip = False

#     if len(line.strip()) == 0 and cumul_store is True:
#         cumul_store = False

#         #outtmp = open('tmp.log', 'w')

#         dx   = cumul_dict[0][1]-cumul_dict[0][0]
#         area = 0
#         cumul_area = 0

#         # Compute overflow:
#         for i in cumul_dict:
#             area = area + cumul_dict[i][2]*dx
#         overflow = totalxs - area

#         # Conversion to no. of events by multiplication with luminosity:
#         # --> assume 20 ab-1 (units are pb):
#         lumi = 20.0e6

#         # Fill cumulative histogram:
#         cumul_entries = len(cumul_dict)

#         for entry in cumul_dict:
#             cumul_y   = overflow
#             cumul_dy2 = 0
#             for i in range(entry,cumul_entries):
#                 cumul_y  = cumul_y + cumul_dict[i][2]*dx*lumi
#                 cumul_dy2 = cumul_dy2 + (cumul_dict[i][3]*lumi)**2

#             cumul_histo[entry]=[cumul_dict[entry][0],
#                                 cumul_dict[entry][1],
#                                 cumul_y,
#                                 math.sqrt(cumul_dy2)]
#             cumul_area = cumul_area + cumul_y

#         # Write it to file:
#         outputfile.write("\n\n")
#         outputfile.write("# " + cumul_name + "\n")
#         for entry in cumul_histo:
#             outputfile.write(('%.9e' % cumul_histo[entry][0]).ljust(5).replace('e','D')  + '\t' +
#                              ('%.9e' % cumul_histo[entry][1]).ljust(5).replace('e','D')  + '\t' +
#                              ('%.9e' % cumul_histo[entry][2]).ljust(10).replace('e','D') + '\t' +
#                              ('%.9e' % cumul_histo[entry][3]).ljust(10).replace('e','D') + '\n')

#         outputfile.write("\n\n")
#         cumul_dict.clear()
#         cumul_histo.clear()
#         #outtmp.close()


#     if len(line.strip()) == 0 and rebin is True:
#         rebin = False

#         #outtmp = open('tmp.log', 'w')
#         dx       = rebin_dict[0][1]-rebin_dict[0][0]
#         dx_rebin = rebin_dict[0+rebin_size-1][1]-rebin_dict[0][0]

#         # Fill cumulative histogram:
#         rebin_entries = len(rebin_dict)
#         entry = 0

#         while entry < rebin_entries:
#             rebin_y   = 0
#             rebin_dy2 = 0

#             for i in range(entry,entry+rebin_size):
#                 rebin_y  = rebin_y + rebin_dict[i][2]*dx
#                 rebin_dy2 = rebin_dy2 + (rebin_dict[i][3])**2

#             rebin_histo[entry]=[rebin_dict[entry][0],
#                                 rebin_dict[entry+rebin_size-1][1],
#                                 rebin_y/dx_rebin,
#                                 math.sqrt(rebin_dy2)/rebin_size]
#             entry = entry + rebin_size


#         # Write it to file:
#         outputfile.write("\n\n")
#         outputfile.write("# rebin" + str(rebin_size) + "_" + histo + "\n")
#         for entry in range(0,len(rebin_dict)): #rebin_histo:
#             if entry in rebin_histo:
#                 outputfile.write(('%.9e' % rebin_histo[entry][0]).ljust(5).replace('e','D')  + '\t' +
#                                  ('%.9e' % rebin_histo[entry][1]).ljust(5).replace('e','D')  + '\t' +
#                                  ('%.9e' % rebin_histo[entry][2]).ljust(10).replace('e','D') + '\t' +
#                                  ('%.9e' % rebin_histo[entry][3]).ljust(10).replace('e','D') + '\n')
#         outputfile.write("\n\n")
#         rebin_dict.clear()
#         rebin_histo.clear()
#         #outtmp.close()


#     if line.startswith(" 0.50000000") and readxs is True:
#         totalxs = float(line.rsplit(' ')[3].strip().replace('D','E'))
#         readxs = False

#     if not line.startswith("#") and cumul_store is True:
#         cumul_dict[cumul_i] = [float(line.split()[0].replace('D','E')),
#                                float(line.split()[1].replace('D','E')),
#                                float(line.split()[2].replace('D','E')),
#                                float(line.split()[3].replace('D','E'))]
#         cumul_i = cumul_i + 1

#     if not line.startswith("#") and rebin is True:
#         rebin_dict[rebin_i] = [float(line.split()[0].replace('D','E')),
#                                float(line.split()[1].replace('D','E')),
#                                float(line.split()[2].replace('D','E')),
#                                float(line.split()[3].replace('D','E'))]
#         rebin_i = rebin_i + 1


# inputfile.close()
# outputfile.close()

# print "----------------------------------"
# print "%i HISTOGRAM FOUND" % int(his_n)
# print "----------------------------------"
