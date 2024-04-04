# THAMES

## CHANGES

| Date                     | Description                         |
|:-------------------------|:-------------------------------------|
| Apr 04 12:13:14 2024 | Changed naming convention of img files, with time string now in minutes. |
| Mar 19 14:49:02 2024 | Better handling of color scheme for visualization of microstructure |
| Mar 19 13:41:50 2024 | Very minor changes to output of pore size distribution in Lattice.cc |
| Mar 19 13:33:18 2024 | viz program now creates one master xyz file with time stamps for Ovito |
| Mar 17 22:08:20 2024 | Removed createmic from main branch |
| Mar 15 14:09:23 2024 | Updated test file input |
| Mar 15 12:55:17 2024 | Merged some improvements made by Florin Nita |
| Mar 14 18:10:02 2024 | Added new createmic program to make THAMES input microstructures directly |
| Feb 22 09:15:56 2024 | Added Cemdata18 database files from EMPA to the repository; updated build instructions |
| Jan 28 15:27:45 2024 | Minor changes to visualization programs |
| Jan 20 15:17:23 2024 | Added StandardKineticModel class for dissolution of salts. Compiles and runs |
| Jan 18 10:16:14 2024 | Cleaned up some debugging output |
| Jan 17 00:58:36 2024 | Debugged behavior for adding growth sites |
| Dec 29 19:21:12 2023 | Formatting uniformity and added viz to build |
| Dec 23 22:34:19 2023 | Streamlined adjustment
of microstructure volumes and output chemical shrinkage in Microstructure csv file |
| Dec 22 13:51:41 2023 | Total microstructure volume now forced constant and
capillary porosity modified to keep it that way; actual volume still tracked |
| Dec 20 16:43:17 2023 | Modified pozzolanic rate law; all solid phases are increasing in volume fraction (needs debug) |
| Dec 19 22:51:08 2023 | Fixed retrieval of saturation index from solution object for pozzolanic kinetic step |
| Dec 19 17:18:17 2023 | Debugging of pozzolanic kinetics; code compiles but does not run correctly |
| Dec 18 23:36:31 2023 | Modified pozzolanic rate equations for diffusion and updated chemistry.xsd |
| Dec 18 16:01:39 2023 | Added LOI influences on clinker component reaction rates |
| Dec 18 14:09:33 2023 | Fixed increment of IC moles for kinetic models and added hydroxyl activity term to pozzolanic model |
| Dec 17 15:11:00 2023 | Added two diffusion rate constants for pozzolanic reactions, akin to the k2 and k3 parameters in the Parrot-Killoh model |
| Dec 15 16:08:51 2023 | Updating kinetic models to have initial specific surface areas for each phase; compile not checked |
| Dec 14 16:28:39 2023 | Fixed run-time errors in KineticController that were not assigning water or solid masses initially |
| Dec 08 09:50:54 2023 | Added more in-depth queries for Rd values and ICs; compiles error-free |
| Dec 07 17:40:56 2023 | Completed moved of Rd partitioning from KineticController to ChemicalSystem; compiles error-free |
| Dec 07 16:17:19 2023 | Moved Rd partitioning from KineticController to ChemicalSystem; compile not checked |
| Nov 11 22:42:13 2023 | Compile-time errors fixed; need to complete Pozzolanic model |
| Nov 10 17:48:17 2023 | Fixed several compile-time errors; some need to still be fixed |
| Nov 03 17:14:04 2023 | Making KineticController handle the models; still no compile check |
| Nov 02 01:34:03 2023 | Completed ParrotKillohModel.cc; still no compile check |
| Nov 01 11:38:14 2023 | Completed ParrotKillohModel.h; still need ParrotKillohModel.cc |
| Nov 01 11:00:08 2023 | Completed KineticModel and KineticController; compile not checked |
| Oct 29 21:13:00 2023 | Setting up multiple kinetic models; compile not checked |
| Oct 27 17:15:10 2023 | Distinguishing pozzolanic kinetics from Parrot-Killoh kinetics; compile not checked |
| Oct 27 16:48:36 2023 | Created branch pozzolan for kinetically modeling pozzolan glasses |
| Aug 10 17:02:31 2023 | Fixed spelling error in Interface.cc |
| Jul 31 12:46:53 2023 | Improved calculation of internal RH using Kelvin equation |
| Jul 30 23:28:36 2023 | Merge DCPorosity branch with master, added output of enthalpy |
| Jun  3 22:55:14 2023 | Specify and read gas phase composition.  Does not do anything with it yet |
| May 23 10:37:38 2023 | Fixed bug in handling of spaces in path names |
| May 23 08:05:26 2023 | Allow path and file names that include spaces |
| May 15 15:41:03 2023 | Updated test file input after rebase |
| May 15 10:52:34 2023 | Updated tag to 2.6.1; updated install instructions |
| May 13 14:05:20 2023 | Fixed minor calculation error in w/c effect of the PK model; cleaned up test cases |
| May  6 16:37:11 2023 | Fix behavior when GEM_run fail happens on first try |
| Apr 25 21:45:58 2023 | Fixed error that mistakenly "corrected" capillary pore volume for its porosity |
| Apr 30 20:59:06 2023 | Fuller implementation of global DEBUG preprocessor define |
| Apr 27 17:47:31 2023 | Removed debug as class variables and now have a global DEBUG preprocessor define |
| Apr 25 21:17:23 2023 | Fixed error that mistakenly "corrected" capillary pore volume for its porosity |
| Apr 24 21:45:58 2023 | Fixed one error in Lattice::writePoreSizeDistribution. Phase porosity still not being calculated |
| Apr 22 21:45:58 2023 | Fixed more runtime errors. Almost working but won't empty pores |
| Apr 21 17:46:33 2023 | Fixed multiple runtime errors, but still has some |
| Apr 20 14:19:42 2023 | Composition-dependent subvoxel porosity compiles; need to test runtime behavior |
| Apr 19 18:21:20 2023 | Created a branch for composition-dependent subvoxel porosity |
| Apr 14 23:54:13 2023 | Accounting for subvoxel porosity of microstructure phases |
| Mar 30 16:21:47 2023 | Adjust microstructure molar volumes based on interhydrate porosity |
| Jan 14 17:12:41 2023 | Corrected local path for GEMS-3K node.h include file |
| Jan  3 13:52:07 2023 | Tagged as version 2.6 |
| Jan  3 13:52:07 2023 | Assume XML schema files are in local working directory |
| Dec 29 23:28:35 2022 | Made silica fume a kinetic phase with a custom dissolution rate equation |
| Dec 14 13:52:07 2022 | Minor formatting changes |
| Jul 28 15:16:48 2022 | Reduced code's dependence on hard-wired phase names |
| Mar 25 19:35:46 2021 | Improved readability of KineticModel and ChemicalSystem classes |
| Mar 22 22:01:18 2021 | Added GEMS3K standalone library as a subtree |
| Mar 11 16:28:57 2021 | Corrected the chemistry.xml file in two of the sample-input examples |
| Mar  9 22:21:43 2021 | Fixed a bug that causes runtime errors when the chemistry.xml file has kinetically controlled phases separated by other types of phases |
| Mar  5 15:17:31 2021 | Cleaned up some comments and added new sample-input examples |
| Mar  3 15:17:31 2021 | Added ability to "tweak" the IC compositions slightly when a GEM run error is encountered, in the hopes of re-running the GEM run method successfully |
| Mar  1 12:02:18 2021 | Water-cement ratio and scaled masses now inferred directly from the microstructure |
| Feb  6 16:56:47 2021 | Better exception handling; better interface to new GEMS3K standalone library |
| Jan 14 19:33:52 2021 | Can now specifiy initial solution composition in input file |
| Jan 11 19:03:08 2021 | Removed growthtemplate functionality from input |
| Jan  8 20:54:08 2021 | Improved vcctl2thames program and integrated within cmake |
| Jan  5 20:14:32 2021 | Added a guide for preparing input files |
| Jan  4 16:40:34 2021 | Added verbose option, plus better control over output images |
| Jan  1 20:28:16 2021 | Better handling of void production under sealed conditions |
| Dec 28 15:38:15 2020 | Simulation should stop if capillary pore water is exhausted |
| Dec 22 11:36:29 2020 | Under sealed conditions, now emptying water from the largest capillary pores |
| Dec 18 13:22:55 2020 | Removed hard-wired code for ettringite crystallization pressure |
| Dec 16 13:15:34 2020 | Minor change to Controller |
| Dec  7 13:25:09 2020 | Added new input and output examples |
| Dec  6 20:16:30 2020 | Added ability to do saturated or sealed curing, plus adjustable wc ratio in the chemistry.xml file alone |
| Dec  3 17:44:43 2020 | Enabling w/c ratio to determine initial water content (still needs work).  Also moved GEMS3K library outside of project |
| Nov 20 17:31:57 2020 | New input files and modified simulation type flags behavior |
| Nov 20 16:00:33 2020 | Damage file only created when sulfate attack simulation |
| Nov 20 13:19:55 2020 | Slight tweak to gitignore file |
| Dec 17 15:23:49 2019 | Minor changes to DISCLAIMER.md |
| Dec 17 15:22:24 2019 | Modified gitignore file |
| Dec 17 15:18:08 2019 | Modified gitignore file |
| Dec 17 15:15:46 2019 | Added a gitignore file |
| Jun 26 15:32:47 2019 | Better time logging in output files and more accurate image generation |
| Apr 23 16:52:47 2019 | Modified lattice dimensions when adding sites (Lattice::addSite) and modified myconfig.h so that it gets customized when cmake is executed |
| Mar 08 17:46:56 2019 | Added build directory |
| Apr 13 14:55:26 2017 | Added table markdown for README.md |
| Apr 13 14:48:12 2017 | Updated formatting of README.md |
| Apr 13 14:45:35 2017 | Testing markdown formatting of bulleted lists |
| Apr 13 14:37:01 2017 | Initial committing of THAMES for GitHub remote repository |
| Apr 13 13:12:36 2017 | Update README.md and corrected name spelling |
| Apr 27 11:22:22 2015 | Initial commit |
