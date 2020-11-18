# Wireless Power Transfer
Algorithms for design of Wireless Power Transfer power pads.
A series of functions developed in Matlab to help designing power pads (or couplers) used in inductive power transfer applications.

If you are not familiar to this branch of Electrical Engineering, "wireless power tranfer" refers to tranfering power between two or more physically separated electrical circuits:
- One is usualy called the emitter or the primary (which is the same name used in transformers);
- The other ones are called receivers or secondaries (again, a reference to transformer circuit);

The theory and descriptions about these algorithms can be found here: **Exhaustive algorithms applied to the design of inductive power transfer couplers**, R.C. Fernandesm A.A. Oliveira Jr., Cambridge University Press, August 2015. DOI: https://doi.org/10.1017/wpt.2015.7

These algorithms are called "exhaustive" for one simple reason: they make dozens, maybe hundreds of finite element simulations, adjusting the pad geometry until it outputs the target power. In simple terms, users will input design parameters (distance between couplers and switching frequency, for example) and select a coupler configuration. An instance of the Finite Element Method Magnetics (FEMM) software will be opened (must be installed previously) and from that point on the algorithms will literally draw the chosen geometry. As long as the desired output parameters (such as power or voltage) are not met, the design dimensions will be slightly changed and re-simulated. The work ends if the software finds a suitable geometry that solves your coupler's target outputs. So, this solves the problem of trial and error approach.

The coupler configurations currently available are:
- The Circular-Circular;
- The Double D (DD);
- The Double D Quadrature (DDQ).



