# LFSM-0.1
## Long-frequency Sky Model (LFSM)
--------
The LFSM implementation of the NE2001 electron model using in 3D emissivity.

The output can be stored either in memory by means of a numpy array, or in a HDF5 format file.

## Installation
```
as root user:
python setup.py install or python setup.py develop
for user only:
python setup.py install --user
```



To use this code, one can simply do::

    >>> from LFSM.absorp_sky_map.global_spectrum_new import absorption_JRZ
    >>> f = absorption_JRZ(v = v, nside = nside, clumping_factor = 1., index_type = 'constant_index_minus_I_E', distance = dist, test = False, emi_form  = 'exp',I_E_form = 'seiffert',R0_R1_equal=True,using_raw_diffuse = False,only_fit_Anu = False)
    >>> f.mpi()
    v int: frequency in MHz of the input skymap  \\
    nside int: the NSIDE of the skymap in healpix mode \\
    clumping_factor float: the clumping factor influnce the absorption, the value set to one in our model. \\
    index_type str: ['constant_index_minus_I_E','freq_dependence_index_minus_I_E','pixel_dependence_index_minus_I_E'], one of them can be choose with different type of spectral index. \\
    distance kpc: the integration max distance for galaxy absorption. \\
    test: normally False, if True will do some test programm. \\
    emi_form str: ['exp','sech'] the distribution form of emissivity, normally choosen 'exponantial'. \\
    I_E_form str: ['seiffert'], the form of extragalactic component except for CMB. \\
    R0_R1_equal bool: in this paper we fixed R0 equals to R1 in emissivity. \\
    using_raw_diffue bool: the input data for fitting parameter of emissivity, if True the data will be smoothed by hp.smooth() method. \\
    only_fit_Anu: fixed False. \\
There are an example in LFSM/test_ground, one can simpyly go to the dir and run test.py
## Documentation
Documentation can be found at `<https://lfsm-01.readthedocs.io/en/latest/>`_.
