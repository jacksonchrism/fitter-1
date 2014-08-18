
 git clone kkamdin/fitter.git somewhere. I'm going to pretend that you do it in ~/fitter so the rest of my directions will make sense

Setting up RAT to run in "aged concentrator mode"
1. clone the rat-tracking-aging rat branch
2. vim src/physics/PMT/AgedConcentratorModel.cc
3. line 119: change path to ~/fitter/AGED_CONCENTRATOR_PARAMS.ratdb (in case you're wondering why this ratdb file isn't in the data directory, it's because you can't write there during run time, so in order to change the fit parameters for each iteration and still pass them to rat, we write to a ratdb file outsdie of the rat directory) *although this had changed, I should update this
4. source env.sh
5. scons

Setting up fitter to run properly 
1. in Makefile: change /home/kate to your home
2. in wwfitter.cpp and MyMinisim.cc change the various paths to match the appropriate path on your machine (sorry, I realize this is not very convenient)
3. source the rat-tracking-aging env.sh file
4. make
5. ./runfitter


