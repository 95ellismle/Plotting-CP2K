# Plotting-CP2K
Will plot the parameters that are outputted by CP2K when doing a CTMQC run.

The code is very `hacked together` and is a bit of a mess. It will be improved/updated periodically as I progress in my PhD.

To run it run the main.py file (e.g. python main.py). The code works with python2.7+ (including python3).

This will both plot and analyse the data giving some nice info such as the best and worst replicas, the coupling used, various bits from the input file, drifts (norm and energy) etc... 

# TODO
* Move all plotting parameters from PLOT to main.py

* Fix the average coupling calculation (there is a unit mmix up and it would be good if only the 'settled' bits of the graph are used in the calculation)
