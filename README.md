---Description---
Package for pulse propagation with parametric gain conversion coupling towards other wavelength using a intermediary decaying wave (vibration wave). Each wavelength is considered an individual pulse (a single wavelength "line"). Developed for Raman gain in gases. Can easily be adapted to Brillouin gain. Can be generalised to any parametric gain with system and waveguide that does not include intermediary wave. 

---Instructions---
The MBE package contains 1 main class (Model.m) and sub-packages for type of waveguide (+Fiber), gain medium (+Medium) and system (+System). The Data sub-package is for internal usage only.

In a first step, the main class (Model) must be loaded. Ex:
model = MBE.Model;

A Fiber, Medium and System class must be assigned as follows:
model.syst = MBE.System.system;
model.fiber = MBE.Fiber.fiber;
model.medium = MBE.Medium.medium; 

Note, in the above example are shown generic parent classes. These do not contain physically valid functions. 
Specific child classes should be selected instead (see choices below).

---Examples---
After, all parameters for mesh, fields, fiber, medium and system must be set. View examples in the following scripts:

ExecuteMBE_back.m
ExecuteMBE_forward.m (used in Loranger et al. (2024), Opt. Express 32(5), 7622-7632)
ExecuteMBE_FBG.m
ExecuteMBE_Brillouin.m (not physically validated)
ExecuteMBE_MM.m (not tested, development in progress)

---Systems---
Available systems are:
General_Fwrd_Raman	: Generalized multi-mode forward Raman propagation. Has not been fully tested.
SM_Fwrd_Raman 		: Basic single-mode forward Raman propagation. A simplified version of General_Fwrd_Raman for faster computation.
SM_Fwrd_Raman_2sect	: Same as SM_Fwrd_Raman, but with 2 sections of fibers
SM_Back_R_Raman		: Same as SM_Fwrd_Raman, but with possible backward-propagating fields from fiber-end reflection.
SM_Back_FBG_Raman	: Same as SM_Back_R_Raman, but with additional reflection anywhere by a Bragg grating.
SM_Brillouin		: Same as SM_Fwrd_Raman, but adapted to Brillouin gain (each line couples to a backward line only). 
			  Has not been physically validated.

---Fibers---
Ideal_S			: Single-mode anti-resonant hollow-core fiber waveguide. Includes AR-HCF dispersion calculation. 
			  Overlap between all lines is consiered equal and ideal (=1) except for first (lowest) line. 
			  Hence, only 2 overlap values are considered: 1) between the lowest line, its own Q and its pump; 
			  2) between the lowest line, its pump's Q and its pump.
Ideal_MM		: Same as Ideal_S, but generalised to multi-modes. Cannot be used with SM_... systems. Can input inter-modal overlap.
			  Has not been fully tested.

---Medium---
gas_H2Xe		: Raman gain with H2 and added Xe buffer.
gas_H2			: Same as H2Xe, but without Xe.
gas_SBS			: Generalised Brillouin gain. Has not been physically validated, values should be revised.

---Data---
Fields			: Field class, contains starting pulse parameter and results.
Mesh			: Mesh class. Input grid parameters and used to generate grid.
cmap_biegert		: data for color-code display
