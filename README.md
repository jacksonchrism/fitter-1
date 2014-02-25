
 git clone kkamdin/fitter.git somewhere. I'm going to pretend that you do it in ~/fitter so the rest of my directions will make sense

Setting up RAT to run in "aged concentrator mode"
1. in the directory of your choosing: git clone git@github.com:kkamdin/rat.git
2. cd rat
3. checkout rat-con
(I guess you can get the rat-con branch directly without the master, but i find this the easiest)
3. vim src/physics/PMT/AgedConcentratorModel.cc
4. line 119: change path to ~/fitter/AGED_CONCENTRATOR_PARAMS.ratdb (in case you're wondering why this ratdb file isn't in the data directory, it's because you can't write there during run time, so in order to change the fit parameters for each iteration and still pass them to rat, we write to a ratdb file outsdie of the rat directory)
5. source env.sh
6. scons

Setting up fitter to run properly 
1. in Makefile: change /home/kate to your home
2. make
3. in AGED_CONCENTRATOR_PARAMS.ratdb: change p0 to say, 0.4
4. ./wwfitter > logs/mylog.txt 2>&1 (is how i usually run it)

This will do a one parameter fit (p1 is fixed to 0.05, and p0 floats) to the angular response in MyHistos_pbomb800cm_386nm_constDegredation0p3_diffuseconst0p05.root
