#/bin/bash
cd testsuit/test
cellMarkerPipe preprocess -wd ./ -10xd ../../data/Zeisel/10x
cellMarkerPipe selection -wd ./ -10xd ../../data/Zeisel/10x -m de
cellMarkerPipe evaluation -wd ./ -np 10
