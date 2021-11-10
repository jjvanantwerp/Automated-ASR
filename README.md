# Automated-ASR
This is a python code that automates the workflow of Ancestral Sequence Reconstruction.
It's written in Python3, and the command lines may be different depnding on the executables for IQTree, MAFFT, and CD-Hit in your system.
The current version (1.0) should also produce libraries of degenerate DNA sequences from each ancestral sequence.
When downloading this software to use for your own analysis, it will be nessecary to change a few things:
  Paste in the sequence you wish to analyze at the top, and call it the variable "sequence"
  Rename the directory to somehthing that is memorable and that is not already a directory in the same workspace as this software. Its current form WILL overwrite existing directories.
  Be sure that CD-Hit, MAFFT, and IQ-Tree are installed, and their executables have been added to the $PATH such that the executable lines will run.
      CD-Hit: https://github.com/weizhongli/cdhit 
      MAFFT: https://mafft.cbrc.jp/alignment/software/
      IQ-TREE: http://www.iqtree.org/#download
  The Biopython module will need to be installed.
  If you have any questions, feel free to reach out - my email is in the header of the python file.
