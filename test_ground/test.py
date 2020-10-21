from LFSM.absorp_sky_map.global_spectrum_new import absorption_JRZ
v =1.
nside = 2**6
dist = 50.

cla = absorption_JRZ(v = v, nside = nside, clumping_factor = 1., index_type = 'constant_index_minus_I_E', distance = dist, test = False, emi_form  = 'exp',I_E_form = 'seiffert',R0_R1_equal=True,using_raw_diffuse = False,only_fit_Anu = False)
cla.mpi()
