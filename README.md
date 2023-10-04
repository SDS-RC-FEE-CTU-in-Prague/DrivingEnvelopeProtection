# DrivingEnvelopeProtection
Implementation of the driving envelope protection algorithm in Simulink with testing environment in IPG Carmaker.

This is a source code of our implementation of the driving envelope protection scheme, described in a paper, which is under review process. Driving envelope protection is implemented in the core as a model predictive controller. As a solver for the optimal control problem, we used qpOASES version 3.2. As an Implementation environment, we used Matlab/Simulink R2021a. As a test environment, we used IPG CarMaker 10.1.

This project also contains the setting of the used car and test environments in CarMaker. We added a technical paper with a short description of provided ride tests. Video results of the experiments are published on our YouTube channel.

The folder 'DEP in C' also contains in implementation of the same code in C/C++. To run this code make sure that you have installed the library qpOASES from https://github.com/coin-or/qpOASES 

If you want to contact authors, please use the following emails:
_denis.efremov@fel.cvut.cz_,
_tomas.hanis@fel.cvut.cz_.
