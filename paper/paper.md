---
title: 'GRDzhadzha: A code for evolving relativistic matter on analytic metric backgrounds'
tags:
  - c++
  - MPI
  - Open MP
  - vector intrinsics
  - gravity
  - general relativity
  - numerical relativity
authors:
- name: Josu C. Aurrekoetxea
  orcid: 0000-0001-9584-5791
  affiliation: 1
- name: Jamie Bamber
  orcid: 0000-0001-7181-3365
  affiliation: 1
- name: Sam E. Brady
  orcid: 0009-0000-5568-839X
  affiliation: 2
- name: Katy Clough
  orcid: 0000-0001-8841-1522
  affiliation: 2
- name: Thomas Helfer
  orcid: 0000-0001-6880-1005
  affiliation: 3
- name: James Marsden 
  orcid: 
  affiliation: 1
- name: Miren Radia
  orcid: 0000-0001-8861-2025
  affiliation: 4
- name: Dina Traykova
  orcid: 0000-0002-3451-0987
  affiliation: 5
- name: Zipeng Wang
  orcid: 0000-0002-4745-8209
  affiliation: 3
affiliations:
- name: Astrophysics, University of Oxford, Denys Wilkinson Building, Keble Road, Oxford OX1 3RH, United Kingdom
  index: 1
- name: School of Mathematical Sciences, Queen Mary University of London, Mile End Road, London E1 4NS, United Kingdom
  index: 2
- name: Henry A. Rowland Department of Physics & Astronomy, Johns Hopkins University, 3701 San Martin Drive, Baltimore, Maryland (MD) 21218, United States
  index: 3
- name: Department of Applied Mathematics and Theoretical Physics (DAMTP), University of Cambridge, Centre for Mathematical Sciences, Wilberforce Road, Cambridge CB3 0WA, United Kingdom
  index: 4
- name: Max Planck Institute for Gravitational Physics (Albert Einstein Institute), Am Mühlenberg 1, Potsdam-Golm, 14476, Germany
  index: 5
date: 14 Aug 2023
bibliography: paper.bib

---

# Summary

Strong gravity environments such as those around black holes provide us with unique opportunities to study questions in fundamental physics (see e.g [@Macedo:2013qea;@Barausse:2014pra;@Barack:2018yly;@Bertone:2019irm]), such as the existence and properties of dark matter and dark energy. Characterising the behaviour of new fields and other types of matter in highly relativistic environments usually necessitates numerical simulations unless one imposes significant symmetries.
Therefore we need to turn to numerical methods to study the dynamics and evolution of the complex systems of black holes and other compact objects in different environments, using numerical relativity (NR).
These methods allow us to split the four-dimensional Einstein equations into three-dimensional spatial hypersurfaces and a time-like direction.
Then if a solution is known at the initial spatial hypersurface, it can be numerically evolved in time, where an analytic solution no longer exists.
Whilst the tools of NR provide the most complete (i.e., approximation free) method for evolving matter in such environments, in many cases of interest, the density of the matter components is negligible in comparison to the curvature scales of the background spacetime metric [@Clough:2021qlv]. In such cases it is a reasonable approximation to neglect the backreaction of the matter environment onto the metric and treat it as fixed (assuming the background itself is stationary or otherwise has an analytic form).

In such cases, one does not need to evolve all the metric degrees of freedom as in NR, but only the additional matter ones. 
It is possible to do this using any NR code in a trivial way by setting the evolution of the metric variables to zero, but this is clearly rather inefficient. This code, GRDzhadzha, directly evolves the matter variables on an analytically specified background.
This significantly speeds up the computation time and reduces the resources needed (both in terms of CPU hours and storage) to perform a given simulation. 
The code is based on the publicly available NR code GRChombo [@Clough:2015sqa;@Andrade:2021rbd], which itself uses the open source Chombo framework [@Adams:2015kgr] for solving PDEs. 

In the following sections we discuss the key features and applications of the code, and give an indication of the efficiencies that can be achieved compared to a standard NR code.

# Key features of GRDzhadzha

GRDzhadzha inherits many of the features of GRChombo and Chombo, but avoids the complications introduced when evolving the metric. The key features are:

- Background metrics: The currently available backgrounds in the code are a static Kerr black hole in horizon penetrating Kerr-Schild coordinates and a boosted black hole in isotropic Schwarzschild coordinates. As the code is templated over the background, it can easily be changed or adapted to other coordinate systems for different problems without major code modification.
- Matter evolution: The code calculates the evolution for the matter variables on the metric background using an ADM decomposition [@Arnowitt:1962hi;@York:1978gql] in space and time - currently we have implemented a real and a complex scalar field as examples of matter types. Again the code is templated over the matter class so that the matter types can be exchanged with minimal modification.
- Accuracy: The metric values and their derivatives are calculated exactly at each point, whereas the matter fields are evolved with a 4th order Runge-Kutta time integration and their derivatives calculated with the same finite difference stencils used in GRChombo (4th and 6th order are currently available).
- Boundary Conditions: GRDzhadzha inherits all the available boundary conditions in GRChombo, namely, extrapolating (extrapolating the field value radially from values within the numerical grid), Sommerfeld (radiative) [@Sommerfeld1912], reflective and periodic. 
- Initial Conditions: The current examples provide initial data for real and complex scalar field matter. Since backreaction is ignored, there are no constraint equations to satisfy in the case of a scalar field, and the initial data can be freely specified.
- Diagnostics:  GRDzhadzha has routines for verifying the conservation of matter energy densities, angular and linear momentum densities, and their fluxes, as discussed in [@Clough:2021qlv;@Croft:2022gks].
- C++ class structure: Following the structure of GRChombo, GRDzhadzha is also written in C++ and uses object oriented programming (OOP) and templating.
- Parallelism: GRChombo uses hybrid OpenMP/MPI parallelism with explicit vectorisation of the evolution equations via intrinsics, and is AVX-512 compliant.
- Adaptive Mesh Refinement: The code inherits the flexible AMR grid structure of Chombo, which provides Berger-Oliger style [@Berger:1984zza] AMR with block-structured Berger-Rigoutsos grid generation [@Berger:1991]. Depending on the problem, the user may specify the refinement to be triggered by the matter or the background spacetime [@Radia:2021smk]. One nice feature is that one does not need to resolve the horizon of the black hole unless matter is present at that location, so for an incoming wave a lot of storage and processing time can be saved by only resolving the wave, and not the spacetime background.

# Statement of Need

As mentioned in the introduction, any numerical relativity code like GRChombo can undertake these simulations. Examples of these include the Einstein Toolkit (http://einsteintoolkit.org/), with its related Cactus (http://cactuscode.org) [@Loffler:2011ay;@Schnetter:2003rb], and Kranc (http://kranccode.org) [@Husa:2004ip] infrastructure used by LEAN [@Sperhake:2006cy;@Zilhao:2010sr] and Canuda (https://bitbucket.org/canuda) [@Witek:2018dmd]. Other notable but non-public codes include BAM [@Bruegmann:2006ulg;@Marronetti:2007ya], AMSS-NCKU [@Galaviz:2010mx], PAMR/AMRD and HAD [@East:2011aa;@Neilsen:2007ua]. Codes such as SPeC [@Pfeiffer:2002wt] and bamps [@Hilditch:2015aba] implement the generalised harmonic formulation of the Einstein equations using a pseudospectral method, and discontinuous Galerkin methods are used in SpECTRE (https://spectre-code.org) [@deppe_nils_2021_4734670;@Kidder:2016hev;@Cao:2018vhw]. NRPy (http://astro.phys.wvu.edu/bhathome) [@Ruchlin:2017com] is a code aimed for use on non-HPC systems, which generate C code from Python, and uses adapted coordinate systems to minimise computational costs. CosmoGRaPH (https://cwru-pat.github.io/cosmograph) [@Mertens:2015ttp] and GRAMSES [@Barrera-Hinojosa:2019mzo] are among several NR codes targeted at cosmological applications (see [@Adamek:2020jmr] for a comparison) and which also employ particle methods. Simflowny (https://bitbucket.org/iac3/simflowny/wiki/Home) [@Palenzuela:2018sly], like CosmoGRaPH, is based on the SAMRAI infrastructure, and has targeted fluid and MHD applications. GRAthena++ [@Daszuta:2021ecf] makes use of oct-tree AMR to maximise scaling.

Whilst there exist many NR codes (both public and private), which can in principle be used to perform simulations of fundamental fields on a fixed BH background, most do not have the efficiency advantages of GRDzhadzha^[As far as we are aware, only NRPy and Canuda offer the same functionality. Some private codes also have such capabilities (see e.g. [@Traykova:2017zrn], based on [@Braden:2014cra]). Other codes may have similar features that are not explicitly separated out, so this makes it difficult to identify them.].
In particular, the fact that the ADM variables and their derivatives are not evolved or stored on the grid saves both a lot of simulation run time, as well as output file storage space.
To get a rough idea of the improvement in storage and CPU hours one can achieve, we performed a short test simulation using GRDzhadzha and compared it simulation performed using the full NR capabilities of GRChombo.
We find that on average GRDzhadzha is 15-20 times faster than GRChombo and requires about 3 times less file storage.
An additional advantage of this code versus using a full NR code, for problems with negligible backreaction, is that here the metric variables are calculated analytically at every point on the grid, which significantly decreases the margin for numerical error, and means that resolution can be focussed on the matter location, and not the spacetime curvature.

It is important to note that whilst backreaction is neglected in the metric calculation, this does not mean that the backreaction effects cannot be calculated. Fixed background simulations provide a first order (in the density) estimate of the gravitational effects caused by the matter, taking into account their relativistic behaviour. This is discussed further in [@Clough:2021qlv] and some examples using the approach are [@Bamber:2020bpu;@Traykova:2021dua;@Traykova:2023qyv]. 

Since the interface and structure of the code is very close to the GRChombo numerical relativity code, it is possible for the results of these fixed background simulations to be used as initial data in full numerical relativity simulations (and vice versa), as was done in [@Bamber:2022pbs]. Therefore if the backreaction is found to be significant due to some growth mechanism, the simulation can be continued in full NR.

# Key research projects using GRDzhadzha

So far the code has been used to study a range of fundamental physics problems:

- Studying the interference patterns in neutrino flavour oscillations around a static black hole [@Alexandre:2018crg]

- Growth of scalar hair around a Schwarzschild [@Clough:2019jpm] and a Kerr [@Bamber:2020bpu] black hole
  
- Determining the relativistic drag forces on a Schwarzschild black hole moving through a cloud of scalar field dark matter [@Traykova:2021dua;@Traykova:2023qyv]

![Formation of overdense tail of scalar dark matter behind a moving BH, due to dynamical friction](Figures/friction.png){width=45%} ![The evolution of scalar clouds around black hole binaries in, the initial conditions were generated with a modified version of GRDzhadzha](Figures/accretion.png){width=45%}
*Some examples of the physics studied with GRDzhadzha. The left image is taken from [@Traykova:2021dua] in which the dynamical friction of light dark matter was studied in the relativistic regime. The right image is from a study of the scalar clouds around black hole binaries in [@Bamber:2022pbs], in which the initial conditions were generated with a modified version of GRDzhadzha.*

- Studying the dynamical friction effects on a Kerr black hole (Magnus effect) [@Wang:inprep]
  
- Superradiance with self-interacting vector field [@Clough:2022ygm] and with spatially varying mass [@Wang:2022hra]
  
- BH mergers in wave dark matter environments [@Bamber:2022pbs]


# Acknowledgements

We thank the GRChombo collaboration (http://www.grchombo.org) for their support and code development work. JB and JM acknowledge funding from a UK Science and Technology Facilities Council (STFC) studentship. JCA acknowledges funding from the Beecroft Trust and The Queen’s College via an extraordinary Junior Research Fellowship (eJRF). KC acknowledges funding from the UKRI Ernest Rutherford Fellowship (grant number ST/V003240/1). SEB is supported by a QMUL Principal studentship. ZW is supported by NSF Grants No.~AST-2006538, PHY-2207502, PHY-090003 and PHY-20043, and by NASA Grants No. 20-LPS20-0011 and 21-ATP21-0010. A few of the projects using this code have received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant Agreement No 693024).

GRDzhadzha users have benefited from the provision of HPC resources from:

 * DiRAC (Distributed Research utilising Advanced Computing) resources under the projects ACSP218, ACSP191, ACTP183, ACTP186 and ACTP316. Systems used include:
   
    - Cambridge Service for Data Driven Discovery (CSD3), part of which is operated by the University of Cambridge Research Computing on behalf of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The DiRAC component of CSD3 was funded by BEIS capital funding via STFC capital grants ST/P002307/1 and ST/R002452/1 and STFC operations grant ST/R00689X/1. DiRAC is part of the National e-Infrastructure.

    - DiRAC Data Intensive service at Leicester, operated by the University of Leicester IT Services, which forms part of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The equipment was funded by BEIS capital funding via STFC capital grants ST/K000373/1 and ST/R002363/1 and STFC DiRAC Operations grant ST/R001014/1. DiRAC is part of the National e-Infrastructure.

    - DiRAC at Durham facility managed by the Institute for Computational Cosmology on behalf of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The equipment was funded by BEIS capital funding via STFC capital grants ST/P002293/1 and ST/R002371/1, Durham University and STFC operations grant ST/R000832/1. DiRAC is part of the National e-Infrastructure.

    - DiRAC Complexity system, operated by the University of Leicester IT Services, which forms part of the STFC DiRAC HPC Facility (www.dirac.ac.uk). This equipment is funded by BIS National E-Infrastructure capital grant ST/K000373/1 and STFC DiRAC Operations grant ST/K0003259/1. DiRAC is part of the National e-Infrastructure.
    
  * Sakura, Cobra and Raven clusters at the Max Planck Computing and Data Facility (MPCDF) in Garching, Germany,

  * PRACE (Partnership for Advanced Computing in Europe) resources under grant numbers 2018194669, 2020225359. Systems used include:

    - SuperMUCNG, Leibniz Supercomputing Center (LRZ), Germany

    - JUWELS, Juelich Supercomputing Centre (JSC), Germany
    
  * the Argonne Leadership Computing Facility, including the Joint Laboratory for System Evaluation (JLSE), which is a U.S. Department of Energy (DOE) Office of Science User Facility supported under Contract DE-AC02-06CH11357.

  * the Texas Advanced Computing Center (TACC) at the University of Austin HPC and visualization resources URL: http://www.tacc.utexas.edu, and the San Diego Supercomputing Center (SDSC) URL: https://www.sdsc.edu, under project PHY-20043 and XSEDE Grant No. NSF-PHY-090003

  * the Glamdring cluster, Astrophysics, Oxford, UK.

# References
