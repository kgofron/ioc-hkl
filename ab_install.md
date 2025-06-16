# Installation steps for EPICS, IOCs (hkl, motorSim) 
## Assumptions
-Ubuntu

## dependencies
sudo apt-get install build-essential

## Directories
```bash
mkdir ~/epics/{modules,IOCs,support,GUI}
```

## EPICS
download EPICS base https://epics-controls.org/
```bash
cd ~/epics
mv ~/Downloads/base-7.0.9.tar.gz .
tar -xzf base-7.0.9.tar.gz
mv base-7.0.9/ base/
cd base/
make -j4
```
in ~/.bashrc, append to last line 
```
export PATH=$PATH:/epics/base/bin/linux-x86_64
export EPICS_CA_AUTO_ADDR_LIST=YES
```

## GUI
download Phoebus https://controlssoftware.sns.ornl.gov/css\_phoebus/
```bash
cd ~/epics/GUI
mv ~/Downloads/pheobus-linux.zip .
unzip phoebus-linux
```
run with ./phoebus-4.7.4-SNAPSHOT/phoebus.sh

## hkl
### download
```bash
cd ~/epics/support
git clone https://repo.or.cz/hkl.git
```
### install dependencies
```bash
sudo apt install gtk-doc-tools autoconf libgtkmm-3.0-dev libyaml-dev gettext autopoint gobject-introspection libtool autoconf-archive debhelper gnuplot-nox gobject-introspection gtk-doc-tools libbullet-dev libg3d-dev libg3d-plugins libgirepository1.0-dev libgl-dev libgsl-dev libgtk-3-dev libgtkglext1-dev libhdf5-dev python3-gi python3-pip elpa-htmlize dvipng libhdf5-dev povray asymptote libhdf5-dev libcglm-dev libinih-dev cabal-install
```
### configure and build
```bash
cd ~/epics/support/hkl
./autogen.sh
./configure --enable-introspection --disable-binoculars
make -j4
sudo make install
```
/usr/local/lib should now have libhkl.so, among others
```bash
export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/usr/local/lib
export GI_TYPELIB_PATH=/usr/local/lib/girepository-1.0 
```

## hkl IOC
### download
```bash
cd /epics/iocs
git clone https://github.com/hkl-projects/ioc-hkl.git
```
### check python executable, install numpy
```bash
which python3
```
should return
/usr/bin/python3
if on <=ubuntu 20.04
```bash
pip3 install numpy
```
if on >=ubuntu 24.04
```bash
sudo apt install python3-numpy
```
### set paths
if epics dir in ~
in RELEASE, set EPICS\_BASE=/home/ab/epics/base (needs to be absolute path)
if epics dir in /
in RELEASE, set EPICS\_BASE=/epics/base
### build
```bash
cd ioc-hkl
make
```

## motorSim IOC
### dependencies
```bash 
sudo apt install libtirpc-dev re2c flex
```

### download
```bash
cd ~/epics/iocs
git clone https://github.com/hkl-projects/motorMotorSim
```
```bash
cd ~/epics/modules
git clone https://github.com/epics-modules/motor.git
git clone https://github.com/epics-modules/asyn.git
git clone https://github.com/epics-modules/sequencer.git seq #renames dir to seq
```
### configuration files
                            motorMotorSim
For motorSim ioc, in ~/epics/iocs/motorMotorSim/configure/
```bash
cp EXAMPLE_CONFIG_SITE.local CONFIG_SITE.local
cp EXAMPLE_RELEASE.local RELEASE.local
```
in ~/epics/iocs/motorMotorSim/configure/RELEASE.local
```
MOTOR=/home/ab/epics/modules/motor
SEQ=/home/ab/epics/modules/seq
MOTOR_MOTORSIM=/home/ab/epics/iocs/motorMotorSim
EPICS_BASE=/home/ab/epics/base
```
MOTOR and SEQ lines should be added at the top of the file, these will not be here by default

in ~/epics/iocs/motorMotorSim/configure/CONFIG\_SITE.local, uncomment BUILD\_IOCS = YES
module: MOTOR MODULE
For motor module, in ~/epics/modules/motor/configure/RELEASE, change:
```
SUPPORT=/home/ab/epics/modules
ASYN=$(SUPPORT)/asyn
EPICS_BASE=/home/ab/epics/base
```
comment out all submodules except motorSim in ~/epics/motor/modules/Makefile
module: ASYN
in ~/epics/modules/asyn/configure/RELEASE
```
SUPPORT=/home/ab/epics/modules
EPICS_BASE=/home/ab/epics/base
```
module: SEQ
in ~/epics/modules/seq/configure/RELEASE
```
EPICS_BASE=/home/ab/epics/base
```
### build
```bash
cd ~/epics/modules/asyn
make -sj
cd ~/epics/modules/motor
make -sj
cd ~/epics/modules/seq
make -sj
cd ~/epics/iocs/motorMotorSim
make -sj
```

### external package for elastic peak intensities
```bash
cd /epics/support
git clone https://gitlab.com/soleil-data-treatment/soleil-software-projects/cif2hkl.git
make
make test
sudo make install
```

or with apt
```bash
sudo apt install cif2hkl
```
then move executable to /usr... #TODO



### run motorSim
run with ./st.cmd in ~/epics/iocs/motorMotorSim/iocs/motorSimIOC/iocBoot/iocMotorSim
can find .bob files in ~/epics/modules/motor/motorApp/op/bob/autoconvert
for example, motor8x.bob

### couple PVs between motorsim and hkl-ioc
copy motor/.../op/.../bob/.../motor8x.bob into op dir of hkl-ioc
adjust macros in phoebus to P => "IOC:motorsimtest_ab:"
add FLNK to records in order to bridge the two iocs (motorMotorSim and ioc-hkl)

change /epics/iocs/motorMotorSim/iocs/motorSimIOC/iocBoot/iocMotorSim/motor.substitutions to:
```
file "$(MOTOR)/db/basic_asyn_motor.db"
{
pattern
{P,     N,  M,        DTYP,         PORT,       ADDR,  DESC,          EGU,      DIR,  VELO,  VBAS,  ACCL,  BDST,  BVEL,  BACC,  MRES,  PREC,  DHLM,  DLLM,  INIT, FLNK}
{IOC:motorsimtest_ab:,  1,  "m$(N)",  "asynMotor",  motorSim1,  0,     "motor $(N)",  degrees,  Pos,  1,     .1,    .2,    0,     1,     .2,    0.01,  5,     100,   -100,  "", "HB3:ioc-hkl:omega_e4c.PROC"}
{IOC:motorsimtest_ab:,  2,  "m$(N)",  "asynMotor",  motorSim1,  1,     "motor $(N)",  degrees,  Pos,  1,     .1,    .2,    0,     1,     .2,    0.01,  5,     100,   -100,  "", "HB3:ioc-hkl:chi_e4c.PROC"}
{IOC:motorsimtest_ab:,  3,  "m$(N)",  "asynMotor",  motorSim1,  2,     "motor $(N)",  degrees,  Pos,  1,     .1,    .2,    0,     1,     .2,    0.01,  5,     100,   -100,  "", "HB3:ioc-hkl:phi_e4c.PROC"}
...
```

then add FLNK field in:

/epics/modules/motor/db/basic_asyn_motor.db
```
record(motor,"$(P)$(M)")
{
	field(DESC,"$(DESC)")
	field(DTYP,"$(DTYP)")
	field(DIR,"$(DIR)")
	field(VELO,"$(VELO)")
	field(VBAS,"$(VBAS)")
	field(ACCL,"$(ACCL)")
	field(ACCU,"$(ACCU=0)")
	field(BDST,"$(BDST)")
	field(BVEL,"$(BVEL)")
	field(BACC,"$(BACC)")
	field(OUT,"@asyn($(PORT),$(ADDR))")
	field(MRES,"$(MRES)")
	field(PREC,"$(PREC)")
	field(EGU,"$(EGU)")
	field(DHLM,"$(DHLM)")
	field(DLLM,"$(DLLM)")
	field(INIT,"$(INIT)")
	field(FLNK,"$(FLNK)")
	field(TWV,"1")
}

```


get this warning if exclude .PROC in motor_substitution (I guess this is needed when forward linking across iocs)

WARNING: Forward-link uses Channel Access without pointing to PROC field
    IOC:motorsimtest_ab:m1.FLNK => HB3:ioc-hkl:omega_e4c
WARNING: Forward-link uses Channel Access without pointing to PROC field
    IOC:motorsimtest_ab:m2.FLNK => HB3:ioc-hkl:chi_e4c

to see details of individual PV
(in ioc-hkl epics shell)
```
dbpr HB3:ioc-hkl:omega_e4c, 10
```
PINI: YES
SCAN: Passive

test with increment motors in phoebus and look to e4c (for now) results
```
caget HB3:ioc-hkl:omega_e4c
caget IOC:motorsimtest_ab:m1

camonitor also helps to check if a PV changes, like
camonitor HB3:ioc-hkl:omega_e4c
```


