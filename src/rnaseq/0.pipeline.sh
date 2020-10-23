#!/bin/bash
conda activate format

./1.download_data.sh
./2.do_align.sh
./3.do_shotgun.sh
./4.do_join.sh
./5.do_dge.sh
./6.do_annotate.sh
