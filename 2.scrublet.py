#Run in Python 3.8.6
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
if len(sys.argv)<2 :
	print('Error:sample name not entered')
	sys.exit()
else:
  sample=sys.argv[1]
  input_dir = './'+sample+'/filtered_feature_bc_matrix'
  os.system('gunzip '+ input_dir+ '/*.gz')
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
  genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1))
  scrub = scr.Scrublet(counts_matrix)
  doublet_scores, predicted_doublets = scrub.scrub_doublets()
  scrub.plot_histogram()
  plt.savefig(input_dir+"/doublet_histogram.png")
  np.savetxt(input_dir+'/doublet_scores.txt',doublet_scores)
  np.savetxt(input_dir+'/predicted_doublets.txt',predicted_doublets)
