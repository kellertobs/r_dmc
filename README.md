**********************************************************************
*****  R_DMC USER MANUAL  ********************************************
**********************************************************************

SOFTWARE PACKAGE: Reactive Disequilibrium Multi-Component method for
computing partial melting and magmatic evolution models in the
mantle.

COPYRIGHT (C) 2015
	  Tobias Keller [tobias.keller@earth.ox.ac.uk] 
	  University of Oxford
	  FOALAB
	  Earth Sciences
	  South Parks Road
	  Oxford, OX1 3AN
	  United Kingdom

LICENSE: This program is free software: you can redistribute it and/or modify 
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

GET STARTED: Use provided example user control routines and
calibration file to plot phase diagrams and compute reactive
equilibration scenarios for a four-component compositional space
calibrated for hydrated and carbonated mantle melting. Add new
calibrations by duplicating and modifying the calibration file. Any
number of components can be used with this model, as long as
calibration parameters and initial concentration are specified for
each component.

MODEL DESCRIPTION: For details on the reactive model refer to:
"Effects of volatiles on melt production and reactive flow in the
mantle",  Tobias Keller & Richard F. Katz, (Journal, Year).

PACKAGE CONTENTS: This distribution contains the following Matlab
routines to compute thermodynamic equilibrium and reactive
equilibration in a multi-component compositional space:

- usr/R_DMC_USER_EX1.m
   Example user control routine set up to plot various phase diagram
   readings for hydrated and carbonated mantle melting.

- usr/R_DMC_USER_EX2.m
   Example user control routine set up to compute reactive
   equilibration for hydrated and carbonated mantle melting.

- cal/R_DMC_CAL_MANTLE.m
   Example calibration file for volatile-rich mantle melting using five
   components, and Simon's law pressure dependence.
   
- cal/R_DMC_CAL_JPET16_REF.m
   Example calibration file for hydrated and carbonated mantle melting
   models in a four-component compositional space as used for reference 
   models in "Effects of volatiles on melt production and reactive flow 
   in the mantle", Tobias Keller & Richard F. Katz, (J Pet, 2016).

- cal/R_DMC_CAL_JPET16_HiR.m
   Example calibration file for hydrated and carbonated mantle melting
   models in a four-component compositional space as used for high-R 
   calibration in "Effects of volatiles on melt production and reactive 
   flow in the mantle", Tobias Keller & Richard F. Katz, (J Pet, 2016).

- cal/R_DMC_CAL_JPET16_LoR.m
   Example calibration file for hydrated and carbonated mantle melting
   models in a four-component compositional space as used for low-R 
   calibration in "Effects of volatiles on melt production and reactive 
   flow in the mantle", Tobias Keller & Richard F. Katz, (J Pet, 2016).

- cal/R_DMC_CAL_MANTLE.m
   Example calibration file for hydrated and carbonated basaltic mantle 
   melting models in a six-component compositional space (for experimentation 
   only, not published).

- src/R_DMC_Equilibrium.m
   Routine computing thermodynamic equilibrium in a multi-component
   compositional space defined in terms of P-dependent component
   melting points, and P.T-dependent component distribution
   coefficients. Routine computes solidus and liquidus temperature, and
   equilibrium melt fraction and phase compositions for a given
   pressure, temperature and bulk composition.

- src/R_DMC_ReactionRates.m
   Routine computing linear kinetic equilibration rates given 
   the disequilibrated melt fraction and phase compositions in a
   multi-component system, along with the thermodynamic equilibrium
   state towards which reactions are driven.

- src/R_DMC_PhaseDiagrams.m
   Routine plotting various phase diagram readings for a calibrated
   multi-component compositional system.

- src/R_DMC_ReactiveEquilibration.m
   Routine computing reactive equilibration for a calibrated
   multi-component compositional system. The routine solves a coupled
   system of ODEs for the time evolution of temperature, pressure,
   melt fraction, and phase compositions using linear kinetic
   equilibration rates. User specified rates of de/compression,
   heating/cooling, and de/compaction are applied to disequilibrate 
   the system and drive reactions.

DISCLAIMER: This distribution of R_DMC was developed and tested
on MATLAB 2014b. Running these routines on older versions of
MATLAB may lead to compatibility issues.