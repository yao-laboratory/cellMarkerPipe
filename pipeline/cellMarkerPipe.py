import os
import argparse
#import pipeline.block
from pipeline import block

__version__ = '0.0.0'

def main():

    parser = argparse.ArgumentParser(
                    prog = 'cellMarkerPipe',
                    description = 'Find marker genes for single cell datasets')
    subparsers = parser.add_subparsers(help='help for subcommand: preprocess, selection, evaluation', dest="command")
    parser.add_argument('--version', action='version', version=__version__)
    
    parser_a = subparsers.add_parser('preprocess', help='Preprocess the 10x data')
    parser_a.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir",required=True)
    parser_a.add_argument('-10xd','--10xdir', type=str, help='10x data directory',dest="datadir",required=True)
    parser_a.add_argument('-nvb','--nvariable', type=int, default=2000,  help='The number of highly variable genes selected for later selection',dest="nvariable")
    parser_a.add_argument('-maR','--maxRNA', type=int, default=2500, help='max number of RNA for each cell',dest="maxRNA")
    parser_a.add_argument('-miR','--minRNA', type=int, default=5, help='min number of RNA for each cell',dest="minRNA")
    parser_a.add_argument('-mam','--maxmt', type=int, default=5, help='max number of MT genes for each cell',dest="maxmt")
    parser_a.add_argument('-np', '--nPCA', type=int, default=10, help='The number of PCA chosen for re-cluster',dest="npca")
    parser_a.add_argument('-res', '--resolution', type=float, default=0.5, help='The solution value used in FindClusters for re-cluster',dest="resolution")
    parser_a.add_argument('-alg', '--algorithm', type=int, default=1, help='The algorithm chosen in FindClusters for re-cluster',dest="algorithm")
    parser_a.add_argument('--cluster', dest='cluster',  action='store_true', help="Do cluster")
    parser_a.add_argument('--no-cluster', dest='cluster', action='store_false', help="Do not do cluster")
    parser_a.set_defaults(cluster=False)
    parser_a.add_argument('--know-marker', dest='know', action='store_true', help="Have a file named Known_marker.csv")
    parser_a.add_argument('--no-know-marker', dest='know', action='store_false', help="Do not have a file named Known_marker.csv")
    parser_a.set_defaults(know=False)
    parser_a.add_argument('--keep-known-marker', dest='keep', action='store_true', help="Keep the knwon markers during screening")
    parser_a.add_argument('--no-keep-known-marker', dest='keep', action='store_false', help="Do not keep the known markers during screening")
    parser_a.set_defaults(keep=False)
    
    parser_b = subparsers.add_parser('selection', help='Select marker genes')
    parser_b.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir",required=True)
    parser_b.add_argument('-10xd','--10xdir', type=str, help='10x data directory',dest="datadir",required=True)
    parser_b.add_argument('-m','--method', type=str, help='Method used for selection',dest="method",required=True)

    parser_c = subparsers.add_parser('evaluation', help='Evaluate selected marker genes')
    parser_c.add_argument('-wd', '--workdir', type=str, help='Working directory',dest="workdir",required=True)
    parser_c.add_argument('-np', '--nPCA', type=int, default=10, help='The number of PCA chosen for re-cluster',dest="npca")
    parser_c.add_argument('-res', '--resolution', type=float, default=0.5, help='The solution value used in FindClusters for re-cluster',dest="resolution")
    parser_c.add_argument('-alg', '--algorithm', type=int, default=1, help='The algorithm chosen in FindClusters for re-cluster',dest="algorithm")
    parser_c.add_argument('--know-marker', dest='know', action='store_true', help="Have a file named Known_marker.csv")
    parser_c.set_defaults(keep=False)
    parser_c.add_argument('-10xd','--10xdir', type=str, help='10x data directory',dest="datadir",required=True)
    


    args = parser.parse_args()
    print(args) 


    if args.command == "preprocess":
        print("Command: preprocess the 10x data ...")
        if not os.path.isdir(args.workdir):
            os.mkdir(args.workdir)
        os.chdir(args.workdir)
        block.preprocess(work_dir=args.workdir, data_dir=args.datadir, nvariable=args.nvariable, max_RNA = args.maxRNA, min_RNA = args.minRNA, max_mt = args.maxmt, nPCA=args.npca, resolution=args.resolution, algorithm=args.algorithm, Cluster=args.cluster, know_marker=args.know, keep_known=args.keep)
    
    elif args.command == "selection":
        print("Command: select the marker genes...")
        block.selection(work_dir=args.workdir, data_dir =args.datadir, method=args.method)

    elif args.command == "evaluation":
        print("Command: evaluate the selected marker genes...")
        block.evaluation(work_dir=args.workdir, nPCA=args.npca, resolution=args.resolution, algorithm=args.algorithm, know_marker=args.know, data_dir=args.datadir)


if __name__ == "__main__":
    main()
