# MassiveBlack-II

Tools for accessing the MassiveBlack-II simulations and derived data products hosted by CMU's Coma cluster.

### Dependencies

* pymysql
* fitsio
* treecorr (https://github.com/EiffL/TreeCorr is required for ED correlations)
* halotools (https://github.com/duncandc/halotools/tree/alignments_devel)
* yaml
* astropy

### Cosmology:

n_s = 0.968

sigma_8 = 0.816

Omega_de = 0.725

Omega_m = 0.275

Omega_b = 0.046

Omega_k = 0.0

h = 0.701

w0 = -1.0

wa = 0.0

### Overview

This is my attempt to create an end-to-end pipeline from simulation products (particles and group information) to galaxy catalogues and summary statistics.
This is still a work in progress. We're currently still quite dependent on various scattered bits of code and data products built (but not documented in any meaningful way that we know of) by Ananth Tenneti.


### Usage

The data processing is handled by a series of calls to the scripts in mbii/pipeline.
Adjustable features are set by yaml files (coma-specific examples can be found in mbii/config).

Assuming the paths in the config files are set correctly, mbii is in the PYTHONPATH and there are no dependencies missing one should be able to do.

#### 1. Build subhalo catalogues from particle data

This step entails reading in the SubFind particle positions in each halo, working out the inertia tensors and saving the flattened results as columns in a FITS file.

`python -m mbii.pipeline.calculate_shapes --config config/fiducial_cat.yaml`

#### 2. Postprocess the catalogues 

This step involves various matching, calculation of derived quantities and application of user-specified cuts.

`python -m mbii.pipeline.construct_catalogue --config config/fiducial_cat.yaml`

#### 3. Symmetrise the catalogues

This is is only useful for specific applications. Spin satellite galaxies about the halo centres, maintaining the relative orientation to those centres. Should leave the central objects unchanged.

Output/input files specified in the `symmetrisation` section of the config file.

`python -m mbii.pipeline.symmetrise --config config/fiducial_cat.yaml`

#### 4. Calculate two point statistics

Reads in a galaxy catalogue in the format produced by step 2. Outputs the result as text files.

`python -m mbii.pipeline.calculate_2pt --config config/fiducial_cat.yaml`

### Other config options

cuts : column names, lower then upper bounds

snapshot : which redshift snaptshot to use. 85 is the lowest z~0.

output : place to write the new postprocessed catalogue

errors :
   
     nsub : number of jackknife volumes to use to calculate errorbars.
    
2pt:

     ctypes : whitespace separated list of two point correlations to calculate. 
             See Mandelbaum et al 2010 (https://arxiv.org/pdf/0911.5347.pdf)
             and Joachimi et al 2010 (https://arxiv.org/pdf/1008.3491.pdf)
             
             Allowed:
             
             ed : orientation - separation vector 3D correlation
             ee : orientation - orientation 3D correlation
             gg : position - position 3D correlation
             iiplus_3d : eq 10a of MB10
             iiminus_3d : eq 10b of MB10
             giplus_3d : eq 8 of MB10
             giminus_3d :
             iiplus_proj : 
             iiminus_proj : 
             giplus_proj : eq 12 of J10
             giminus_proj : 

    errors : Whether or not to calculate jackknife errorbars, bool
   
    shapes : location of the input shape/position catalogue
   
    savedir : place to save the output text files
   
    split : quantity by which to split the catalogue before calculation (for example, for central/satellite separation).
   
    split_val : value of the above quantity about which to split the catalogue
   
    rmin : minimum 3D separation
   
    rmax : maximum 3D separation
   
    nbin : number of logspaced bins to use
   
    rpmin : minimum perpendicular separation
   
    rpmax : maximum perpendicular separation
   
    nrpbin : number of logspaced rp to use
   
    pi_max : upper limit of line of sight integration for projected statistics
    
    



## Historical Stuff

### The complete galaxy properties of MassiveBlack II simulation

For each posted redshift it is a complete rip-off of the simulation
data.

We have position velocity id of every particle that is in a group;
and also the often used properties of groups, plus predicted band luminosity
for various filter sets.

Binary files are in intel endian.

Recommend to use mmap to access some of the huge files;
for example, the star properties at 085 are ~10 GB in size.

readsubhalo.py implements SnapDir objects which can be used to access
the data in python. 

==================

KNOWN PROBLEMS
+  The old (before Sep 17 2013) band filters were off by a factor of wavelength -- do not use them

================


Data model: 
----------------
           id
particles  0 1   2 3 4    5 6 7 8 9 10 11     12 13 14 ....
subhalo    |0--| |1----|  |2-contamination|   |3-----|  ....
group      [0-----------------------------]   [1------------

Particles are first cut into groups.
Within groups, gravitationally bound objects are subhalos. 
unbound particles are contaminations.

directory structure:
--------------------

The isolated subhalos are stored in 
  {snapshotid}/

Each archive has the following files/dirs:

meta data:

    /header.txt

    the value of flag_double(header) is important. It decides the precision of
    the particle properties
    
    the value of redshift(header) is also handy. 

group tab:
    grouphalotab.raw is a C struct array of the following type:

    use SnapDir.readgroup() to read in the tab.

    groupdtype = numpy.dtype([
        ('mass', 'f4'), 
        ('len', 'i4'), 
        ('pos', ('f4', 3)), 
        ('vel', ('f4', 3)), 
        ('nhalo', 'i4'),   
        ('massbytype', ('f4', 6)),  
        ('lenbytype', ('u4',6)), 
       ])

    One group per entry, ordered by number of particles

subhalo tab:
    subhalotab.raw is a C struct of the following type:

    one subhalo per entry, plus one entry for the contamination
    For example, if the first group has nhalo = 100, then there will
    be 101 subhalo entries. 

    This is indicated in subhalo/type.raw (u1). If the type is 0
    then entry is a subhalo. if the type is 1, the entry is a contamination.

        load('subhalo', 'type') reads in the types.

    ** Be aware of the contaminations for subhalo analysis **

        load('subhalo', 'tab') reads in the tab.

    subdtype = numpy.dtype([
        ('mass', 'f4'),  # mass, nan for contamination
        ('len', 'i4'), 
        ('pos', ('f4', 3)), 
        ('vel', ('f4', 3)),   
        ('vdisp', 'f4'),  
        ('vcirc', 'f4'),  
        ('rcirc', 'f4'),  
        ('parent', 'i4'),  
        ('massbytype', ('f4', 6)),  
        ('lenbytype', ('u4',6)), 
        ('unused', 'f4'),    
        ('groupid', 'u4'), 
       ])
     
    the luminosity of subhalos at various band are also available (described in
        next section).

    For easying the access, we have performed some joint queries and provide
    the following extra properties of subhalos:

     ( use with SnapDir.load('subhalo', propertyname) )

    subhalo/type.raw  : int8 array of type of the entry: 
                        0 for subhalo 
                        1 for contamination
    subhalo/sfr.raw  : float32 array of total star formation of each entry
    subhalo/bhmass.raw : float32 array of total bhmass of each entry
    subhalo/bhmdot.raw : float32 array of total bhmdot of each entry
    subhalo/bhmassive.raw : property of most massive bh of each entry
    subhalo/bhluminous.raw : property of most luminous bh of each entry

    the type of entries in bhmassive and bhluminous are:
        bhdtype = numpy.dtype([
            ('pos', ('f4', 3)),
            ('vel', ('f4', 3)),
            ('id', 'u8'),
            ('bhmass', 'f8'),
            ('bhmdot', 'f8')])
        bhmass is NaN if there is no bh in that subhalo or contamination.
    
Particle properties:
    /{particletype}/{blockname}.raw
    
    use SnapDir.load(particletype, blockname)
       [optionally give g=subhalotab or g=grouptab to access groups/subhalos], see doc.

    Particle types are 
       0: Gas
       1: Dark mater    
       2&3: unused
       4: Star
       5: BH

    The following blocks have fixed data type:
        ('pos', ('f4', 3)),
        ('vel', ('f4', 3)),
        ('mass', 'f4'),
        ('id', 'u8'),
        ('type', 'u1')

    The other blocks, currently including 
        0/sfr.raw 
        4/met.raw 
        5/bhmass.raw
        5/bhmdot.raw

    are floating numbers, width of witch decided by 
    the value of flag_double in header.txt.
    (simply match the flag_double line)
    float64 : if flag_double is 1
    float32 : if flag_double is 0

    One particle per entry in any of the blocks.
    The particles are ordered in chunks; one chunk per subhalo/contamination.
    Look at the lenbytype value in subhalo tab or group tab.

    For example, 
    if the first group has lenbytype == [50000, 100000, 0, 0, 50000, 100],
    and the second group has lenbytype == [10000, 20000, 0, 0, 10000, 10],

    0:50000 in 0/*.raw are gas in the first group,
    0:100000 in 1/*.raw are dark matter in the first group,
    0:50000 in 4/*.raw are star in the first group,
    0:100 in 5/*.raw are black hole in the first group.

    50000:60000 in 0/*.raw are gas in the second group,
    100000:12000 in 1/*.raw are dark matter in the second group,
    50000:60000 entries in 4/*.raw are star in the second group,
    100:110 entries in 5/*.raw are black hole in the second group.

Subhalo spectra SED and band luminosity :

    The full SED spectra of subhalos is mocked according to Wilkins et al. 2013.
    use SnapDir.readSED to read the SED for a given subhalo. Note that
    if the subhalo has no stars there is no SED and the function returns 0
    the returned SED has nebula lines.

    We also have the band luminosity of the subhalos, integrated from the SED.

    For subhalos with no stars, and for the subhalo entry corresponding to
    contamination, all bands are zero.

    Files are named by for example sbuhaloHST.ACS.raw.
    Refer to readsubhalo.py and bandfilterdtype for all the available bands.

    Two categories of filters `RfFilter' for restframe, and 
    `ObsFilter' for observed frame.

    Subhalos:
    --------
    The band luminosity can be loaded like other 'subhalo' variables, eg, 

        SnapDir.load('subhalo', 'RfFilter/SDSS.SDSS/i')

        (saved  under RfFilter/SDSS.SDSS/i , in 'f4')

        They are in f4 of and of unit 1e28 Erg/s/Hz/h.

        L_band = nu_mean * L_nu_mean / Delta mu

        The magnitude of the band is given by:
        
        M = -2.5 * numpy.log10(L) + 34.1  # L in SI, Joule

    Star particles:
    ----
    The band luminosity of individual star particles can be loaded via

        (eg)
        SnapDir.loadstarLband('RfFilter', 'SDSS.SDSS/r')
    
    IGM absoprtion and nebula lines are included.

    The bands are saved in a compressed format:
     1. The filtered band luminosity per unit mass of 260 time/metallicity 
        bins are calculated and saved in SED/SEDtableXXXX.
     2. 4/SEDindex.raw saves the index of the table entry that saves the
        band luminosity per unit mass.   (Erg/s/Hz/[1e10Msun])
     3. Multiply the entry by the mass of the star particle
        gives its true band luminosity.

     See Snapdir.loadstarLband for more details.
