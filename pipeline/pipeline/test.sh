#!/bin/bash

start=`date +%s`

sleep 5

end=`date +%s`

runtime=$((end-start))

echo $runtime second
