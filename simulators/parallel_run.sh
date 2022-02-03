#!/bin/bash

#simple example showing how to run 10 versions in parallel
#the sleep command is only here temporarily since the sims currently
#can't set their own output directory
for i in {1..10}; do ./mt_dgg_simulator >> test$i & sleep 2; done 
