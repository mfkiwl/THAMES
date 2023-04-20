# THAMES

## CHANGES

| Date                     | Description                         |
|:-------------------------|:-------------------------------------|
| May 15 10:52:34 2023 | Updated tag to 2.6.1; updated install instructions |
| May 13 14:05:20 2023 | Fixed minor calculation error in w/c effect of the PK model; cleaned up test cases |
| May 6 16:37:11 2023 | Fix behavior when GEM_run fail happens on first try |
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

-----------------------------------------------------------------------------
