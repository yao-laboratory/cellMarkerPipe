
from pipeline import *
import os
#import pipeline.block
from pipeline import block

# input data directory: 10X
data_dir = "../../data/Zeisel/10x"
# output data directory
work_dir = "./output"
if not os.path.isdir(work_dir):
    os.mkdir(work_dir)
os.chdir(work_dir)
# 1st step: preprocess
block.preprocess(data_dir=data_dir, work_dir=work_dir, nvariable=2000, Cluster=False, max_RNA = 2500, max_mt = 5)
# 2nd step: selection
block.selection(work_dir=work_dir, data_dir ="", method="de")
# 3rd step: evaluation
block.evaluation(work_dir, nPCA=10)
