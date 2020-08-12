adp_database_fit


# Note: This is test version. (u and w function and ADP_code.init are not good !  Please, you would develop them)


# Debian 10.0 or Ubuntu 18.04 LTS
## Install (cif2cell-informal)
1. sudo apt install -y git python python-setuptools python-dev gcc
2. git clone https://github.com/by-student-2017/cif2cell-informal.git
3. cd cif2cell-informal
4. tar zxvf PyCifRW-3.3.tar.gz
5. cd PyCifRW-3.3
6. sudo python setup.py install
7. cd ..
8. sudo python setup.py install


## Install (lammps)
1. cd ~
2. sudo apt update
3. sudo apt install -y gcc g++ build-essential gfortran libopenblas-dev libfftw3-dev libopenmpi-dev wget
4. wget https://lammps.sandia.gov/tars/lammps-stable.tar.gz
5. tar zxvf lammps-stable.tar.gz
6. cd lammps-3Mar20
7. mkdir build ; cd build 
8. cmake -D BUILD_MPI=on -D PKG_USER-MEAMC=on -D PKG_MANYBODY=on -D PKG_MC=on ../cmake
9. cmake --build .


## Environment settings (lammps)
1. echo ' ' >> ~/.bashrc
2. echo '# lammps environment settings' >> ~/.bashrc
3. echo 'export PATH=$PATH:$HOME/lammps-3Mar20/build' >> ~/.bashrc
4. bash


## Install (PWscf)
1. cd ~
2. sudo apt update
3. sudo apt install -y gcc g++ build-essential gfortran libopenblas-dev libfftw3-dev libopenmpi-dev wget
4. wget https://github.com/QEF/q-e/archive/qe-6.4.1.tar.gz
5. tar zxvf qe-6.4.1.tar.gz
6. cd q-e-qe-6.4.1
7. wget https://github.com/QEF/q-e/releases/download/qe-6.4.1/backports-6.4.1.diff
8. patch -p1 --merge < backports-6.4.1.diff
9. ./configure
10. make pw
11. sudo make install


## set fitting 
1. cd ~
2. sudo apt update
3. sudo apt install -y git python-pip python-scipy csh gfortran gnuplot
4. git clone https://github.com/by-student-2017/adp_database_fit.git
5. cd ~/adp_database_fit
6. gfortran -o Zhou04_ADP_1 Zhou04_create_adp.f
7. gfortran -o plot_adp plot_adp.f


## fit parameters by other methods
1. cd ~/adp_database_fit
2. cp ADP.input_temp ADP.input
3. sed -i 's/Xx/Cu/g' ADP.input
4. python nm_v4_debian.py


## fit parameters by genetic algorithm
1. cd ~/adp_database_fit
2. pip install -U deap --user
3. cp ADP.input_temp ADP.input
4. sed -i 's/Xx/Cu/g' ADP.input
5. python ga_v4_debian.py


## fit parameters by baysian method
1. cd ~/adp_database_fit
2. pip install bayesian-optimization==1.1.0
3. cp ADP.input_temp ADP.input
4. sed -i 's/Xx/Cu/g' ADP.input
5. python baysian_v4_debian.py


## plot F, rho or z2r
1. gnuplot png_adp.gp


# Google Colaboratory
## Install (cif2cell-informal)
	!apt update
	!apt install -y git python python-setuptools python-dev gcc
	%cd /content
	!git clone https://github.com/by-student-2017/cif2cell-informal.git
	%cd cif2cell-informal
	!tar zxvf PyCifRW-3.3.tar.gz
	%cd PyCifRW-3.3
	!python2 setup.py install
	%cd /content/cif2cell-informal
	!python2 setup.py install


## Install (lammps)
	!apt update
	!apt install -y gcc g++ build-essential gfortran libopenblas-dev libfftw3-dev libopenmpi-dev wget
	%cd /content
	!wget https://lammps.sandia.gov/tars/lammps-stable.tar.gz
	!tar zxvf lammps-stable.tar.gz
	%cd lammps-3Mar20
	!mkdir build
	%cd build 
	!cmake -D BUILD_MPI=on -D PKG_USER-MEAMC=on -D PKG_MANYBODY=on -D PKG_MC=on ../cmake
	!cmake --build .
	import os
	os.environ['PATH'] = "/content/lammps-3Mar20/build:"+os.environ['PATH']


## Install (PWscf)
	!apt update
	!apt install -y gcc g++ build-essential gfortran libopenblas-dev libfftw3-dev libopenmpi-dev wget
	%cd /content
	!wget https://github.com/QEF/q-e/archive/qe-6.4.1.tar.gz
	!tar zxvf qe-6.4.1.tar.gz
	%cd q-e-qe-6.4.1
	!wget https://github.com/QEF/q-e/releases/download/qe-6.4.1/backports-6.4.1.diff
	!patch -p1 --merge < backports-6.4.1.diff
	!./configure
	!make pw
	import os
	os.environ['PATH'] = "/content/q-e-qe-6.4.1/bin:"+os.environ['PATH']


## set fitting
	!apt update
	!apt install -y git python-pip python-scipy csh gfortran gnuplot
	%cd /content
	!git clone https://github.com/by-student-2017/adp_database_fit.git
	%cd /content/adp_database_fit
	!gfortran -o Zhou04_ADP_1 Zhou04_create_adp.f
	!gfortran -o plot_adp plot_adp.f


## fit parameters by other methods
	%cd /content/adp_database_fit
	import os
	os.environ["OMP_NUM_THREADS"] = "1,1"
	os.environ["MKL_NUM_THREADS"] = "1"
	!cp ADP.input_temp ADP.input
	!sed -i 's/Xx/Cu/g' ADP.input
	!python2 nm_v4_gc.py


## fit parameters by genetic algorithm
	!pip2 install -U deap --user
	%cd /content/adp_database_fit
	import os
	os.environ["OMP_NUM_THREADS"] = "1,1"
	os.environ["MKL_NUM_THREADS"] = "1"
	!cp ADP.input_temp ADP.input
	!sed -i 's/Xx/Cu/g' ADP.input
	!python2 ga_v4_gc.py


## fit parameters by baysian method
	!pip2 install bayesian-optimization==1.1.0
	%cd /content/adp_database_fit
	import os
	os.environ["OMP_NUM_THREADS"] = "1,1"
	os.environ["MKL_NUM_THREADS"] = "1"
	!cp ADP.input_temp ADP.input
	!sed -i 's/Xx/Cu/g' ADP.input
	!python2 baysian_v4_gc.py


## plot F, rho or z2r
	!gnuplot png_adp.gp


