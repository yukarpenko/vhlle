*******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            by Iurii Karpenko                                                *
*  For the detailed description please refer to:                              *
*  Comput. Phys. Commun. 185 (2014), 3016    [arXiv:1312.4160]                *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appears in all the copies. If you decide to make modifications   *
*  to the code, please contact the authors, especially if you plan to publish *
*  the results obtained with such modified code. Any publication of results   *
*  obtained using this code must include the reference to:                    *
*  Comput. Phys. Commun. 185 (2014), 3016    [arXiv:1312.4160]                *
*                                                                             *
*******************************************************************************


 This repository contains the source code of vHLLE.

The following branches may be particularly useful:
+ cpc_final branch corresponds to the code published in Comput. Phys. Commun.
  This branch is more useful to get familiarized with the structure of the code
  and reproduce the published test results.
+ main branch (RECOMMENDED) contains the most recent, feature-rich code
  used for the actual event-by-event simulations, and includes:
  - different initial state modules (classes)
  - Cornelius module to reconstruct the elements of particlization hypersurface
    (please cite P. Huovinen, H. Petersen, Eur. Phys. J. A (2012) 48: 171,
     https://doi.org/10.1140/epja/i2012-12171-9  when using this branch)
  - improved scheme for ideal-viscous splitting
  - cross terms for shear/bulk viscosity according to
    J. Bernhard et al, Phys.Rev. C 94, 024907, G.S. Denicol et al, Phys.Rev. C 90, 024912
    See cross_terms branch for more details
  - dynamic grid resize for CPU time-saving simulation
  - and many more improvements.

 0. BRIEF DESCRIPTION OF THE PACKAGE

 The package contains the following files:
 Makefile
 src/      : source subdirectory
  |- inc.h : predefined constants
  |- rmn.h, rmn.cpp : transformation procedures from conserved variables to
  |                   primitive variables and back
  |- cll.h, cll.cpp : Cell class which stores and manipulates the properties of
  |                   the individual hydro cell
  |- eos.h, eos.cpp : EoS class, equation of state information
  |- fld.h, fld.cpp : Fluid class, which contains fluid (3D array of fluid cells)
  |                   as a whole and performs auxiliary actions on a fluid
  |- ic.h, ic.cpp : initial conditions for fluid-dynamical evolution from
  |                 optical Glauber or Gubser flow
  |- s95p.h, s95p.cpp: manages tabular initial conditions and EoS 's95' (for optional use)
  |- trancoeff.h, trancoeff.cpp : Trancoeff class containing transport
  |                               coefficients: (T-dependent) shear, bulk
  |                               viscosities and corresponding relaxation times
  |- hdo.h, hdo.cpp : Hydro class, which contains the algorithms of hydrodynamic solution
  \- main.cpp : main() function, initializes all the objects and calls evolution
                loop, also takes care of reading the parameters of the
                simulation from parameter file
 params/  : subdirectory with parameter files, see their description below
 eos/Laine_nf3.dat : equation of state, latticeQCD inspired results
                     from M. Laine and Y. Schroder, PRD 73, 085009
 ic/      : sample tabulated initial conditions, used for 3D run
 hydro_output/  : subdirectory containing hydrodynamic output (generated by the
                  code) and Gnuplot scripts to make some of the plots presented
                  in the paper. See the description in Sect. II below.
 doc/     : Latex drafts documenting some of the improvements from the version
            published in Comput.Phys.Commun.


 I. GETTING and BUILDING vHLLE on Linux

 1) The following software must be installed in order to compile the code:
 make, g++, binutils, GSL

 To install those packages,
 -> on Ubuntu run:
 sudo apt install make g++ binutils gsl-bin
 -> on Fedora run:
 sudo yum install make gcc binutils gsl

 ROOT is not required to compile and run the code.
 As of commit 589b7cb , C++17 is required (because of calls to std::filesystem),
 therefore make sure `root-config --cflags` returns '-std=c++17' among the options.

 Optionally, to run Gnuplot scripts provided in the program package one has to
 install Gnuplot and awk (gawk).

 2) clone the main vHLLE repository:
 git clone https://github.com/yukarpenko/vhlle.git
 cd vhlle
 git checkout main
 mkdir obj
 make

 3) clone another repository, which contains the equation of state,
    sample initial state tables and related sample parameter files:
 cd ..
 git clone https://github.com/yukarpenko/vhlle_params.git
 cd vhlle_params

 4) Copy the eos/, ic/ and params/ subdirectories from vhlle_params/
 into the directory containing the vhlle code, vhlle/, e.g. by executing:
 mv ic eos ../vhlle/
 mv params/* ../vhlle/params/

 II. RUNNING vHLLE

 1. to run vHLLE, cd into the vhlle subdirectory and type:
 ./hlle_visc -params <parameter-file> [-system <system>] -ISinput <IS-file> -outputDir <output-directory>

<parameter-file> is the file name (including relative path) of a parameter file.
<system> is an optional command-line parameter, which is mandatory for GLISSANDO
 and TrENTo initial state options, and specifies a setup for a particular collision
 energy. The standard values are: RHIC200, LHC276, which corresponds to sqrt(s)=200 GeV
 RHIC and sqrt(s)=2760 GeV LHC energies, respectively. See the implementation of
 IcGlissando::IcGlissando() in src/icGlissando.cpp for more details.
<IS-file> is the file name of an initial state table. This file has to be supplied
 for all but the simplest initial state options (optical Glauber or Gubser).
<output-directory> is the directory to write the output files into. If the directory
 does not exist, it will be created. If not specified a 'data' subfolder in the current directory
 will be created or used if already existing.

 A sample command:
 ./hlle_visc -params params/glissRHIC/gliss2RHIC.20-30 -system RHIC200 -ISinput ic/glissando/sources.RHIC.20-30.dat -outputDir output/RHIC.20-30
 executes vHLLE with an averaged initial state table pre-computed from GLISSANDO code, and corresponding to 20-30% central Au-Au collisions at sqrt(s)=200 GeV.

Other sample parameter files can be found in a code package linked to the reference
publication [Comput. Phys. Commun. 185 (2014), 3016]. The parameter files are located
in params/ subdirectory and cover several of the results published in the ref. publication:

 1) comparison to VISH2+1 code in Section 4.3 ("Matter expansion in heavy ion
  collisions"). The corresponding parameter files are:
params/song2DCPC.*
The filename suffixes correspond to:
 b0 = impact parameter 0 fm
 b7 = impact parameter 7 fm
 etas008 = viscous hydro simulation with eta/s=0.08
 ideal = ideal hydro simulation

 In order to reproduce the comparison plots (Figures 9, 10) one should run the
 code with all 4 input files. This will create subdirectories "song2D.*" with
 hydrodynamic output for each parameter set. Then Gnuplot scripts hydro_output/*.plot
 will read and parse the hydrodynamic output and draw the corresponding plots.
 Then run Gnuplot with "hydro_output" as working directory, and load one of the
 scripts in Gnuplot prompt:

 cd hydro_output
 gnuplot
 gnuplot> load 'radFlow.plot'  # which creates "vradSong.eps" postscript containing Figure 9,
 gnuplot> load 'epsilonp.plot'  # which creates "epsilonpSong.eps" postscript containing Figure 10

 2)numerical solution for ideal Gubser flow (end of Section 4.1).
  The corresponding parameter file is:
params/gubserCPC

 To produce the plots run Gnuplot with "hydro_output" as working directory. In Gnuplot prompt, type:
 gnuplot> load 'gubser.plot'
 which creates "gubserEps.eps" and "gubserVx.eps", which are Figures 5 and 6,
 respectively.

 3)3D hydro simulation described in Section 4.4 ("Energy conservation").
  The corresponding parameter file is:
params/3DCPC

 Note that this simulation uses full 3D grid. With the parameters provided, the
 program consumes about 1900 Mbytes of RAM. Therefore please make sure that the
 host machine has enough RAM to run this simulation.

 The purpose of this simulation is to check the energy/entropy conservation
 during 3D hydrodynamic expansion of a closed system. The corresponding results
 are presented in Table 1 of the paper. To reproduce the results, please watch
 the console output of the program, which contains total energy and total entropy
 after each timestep.
 Note that by default the code uses the equation of state 'p=epsilon/3'.
 To reproduce the numbers in the middle column of the Table 1 corresponding to
 'Laine' EoS, please recompile the code with  'EXTRA_FLAGS   = -D TABLE' in Makefile.
 The typical console output of the program is the following:

-----
.....
step= 0  dtau= 0.05

calcTotals: E =        1402.31  Efull =        1402.33
           Px =              0      S =        6869.23
.....
------
where
      E : total energy of the system [GeV], computed from ideal part of the
          energy-momentum system only,
  Efull : total energy of system including viscous terms in energy-momentum tensor,
     Px : x-component of the total momentum
      S : total entropy including viscous terms


 III. STRUCTURE OF THE OUTPUT FILES

 Each hydrodynamic simulation records its results in a separate directory, which
 is specified by the value of 'outputDir' parameter in the corresponding parameter
 file. Below we denote the cells used for the output as  cell(ix,iy,iz), where
 ix=[0...nx-1]  is x-coordinate of cell on 3D grid
 iy=[0...ny-1]  is y-coordinate of cell on 3D grid
 iz=[0...nz-1]  is z-coordinate of cell on 3D grid

 As a result of the simulation the output directory contains the following files:
 outx.dat    :   distributions in cells along X direction at every timestep.
                 The cells used are cell(ix,ny/2,nz/2),  ix=0...nx-1
 outy.dat    :   distributions in cells along Y direction at every timestep.
                 The cells used are cell(nx/2,iy,nz/2),  iy=0...ny-1
 outz.dat    :   distributions in cells along Z direction at every timestep.
                 The cells used are cell(nx/2,ny/2,iz),  iz=0...nz-1
 outdiag.dat :   distributions in cells along diagonal direction in XY plane
                 at every timestep. The cells used are cell(ix,ix,iz), ix=0...nx-1
                 Works reliably if nx=ny
 out.aniz.dat:   grid-integrated quantities $\epsilon_p$, $\epsilon_p'$
                 (see Subsection 4.3 "Matter expansion in ...")

 The format of the columns in the output files is:

 outx.dat:
 t  x  vx  vy  eps  nb  T  mub  [pi*10]  Pi  cut_flag

 outy.dat:
 t  y  vy  vx  eps  nb  T  mub  [pi*10]  Pi  cut_flag

 outz.dat:
 t  z  vz  vx  eps  nb  T  mub  [pi*10]  Pi  cut_flag

 outdiag.dat:
 t  sqrt(x*x+y*y)  vx  vy  eps  nb  T  mub  [pi*10]  Pi  cut_flag

 out.aniz.dat:
 t  <vt>  <epsilon_p>  <epsilon_p'>

 where
 t : proper time [fm/c]
 x : x coordinate [fm]
 y : y coordinate [fm]
 z : rapidity
 vx : x-component of 3-velocity
 vy : x-component of 3-velocity
 vz : longitudinal flow rapidity
 eps : energy density in fluid rest frame [GeV/fm^3]
 nb : baryon density in fluid rest frame [1/fm^3]
 T  : temperature [GeV]
 mub : baryon chemical potential [GeV]
 [pi*10] : $\pi^{mu\nu}$ components: $\pi^{\tau\tau}$, $\pi^{\tau x}$,
           $\pi^{\tau y}$, $\pi^{\tau\eta}$, $\pi^{xx}$, $\pi^{xy}$, $\pi^{x\eta}$,
           $\pi^{yy}$, $\pi^{y\eta}$, $\pi^{\eta\eta}$. All the components
           correspond to $\tilde{\pi^{\mu\nu}}$, see Section 2. Units are [GeV/fm^3]
 Pi : bulk pressure [GeV/fm^3]
 cut_flag : viscous corrections are cut by cut_flag factor in the cell
            (not cut if cut_flag = 1.0)
