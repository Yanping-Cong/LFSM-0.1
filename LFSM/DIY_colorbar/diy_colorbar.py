import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np; np.random.seed(1)
import healpy as hp

def diy_colorbar(_min,_max,unit=r'log$\tau$',fontsize=12,tick_numbers=8):
    #tick_numbers change how many ticks in colorbar
    #size decide the distance between map and colorbar;shrink change the length of colorbar;aspect change the height of colorbar
    hp.graticule(30,coord="G",local=True)
    #fontsize = 'medium'
    hp.projtext(0., 0, r'$0\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(30., 0., r'$30\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(60., 0., r'$60\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(90., 0., r'$90\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(120., 0., r'$120\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(150, 0., r'$150\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(180., 0., r'$180\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(210., 0, r'$210\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(240., 0, r'$240\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(270., 0, r'$270\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(300., 0, r'$300\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(330., 0, r'$330\degree$', lonlat=True, coord='G',color = 'w',fontsize = fontsize)
    hp.projtext(0., 30., r'$30\degree$', lonlat=True, coord='G',color = 'black',fontsize = fontsize)
    hp.projtext(0., 60., r'$60\degree$', lonlat=True, coord='G',color = 'black',fontsize = fontsize)
    hp.projtext(0., -30, r'$-30\degree$', lonlat=True, coord='G',color = 'black',fontsize = fontsize)
    hp.projtext(0., -60, r'$-60\degree$', lonlat=True, coord='G',color = 'black',fontsize = fontsize)



    fig = plt.gcf()
    ax = plt.gca()
    image = ax.get_images()[0]


    #size decide the distance between map and colorbar;shrink change the length of colorbar;aspect change the height of colorbar
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size="20%",pad=0.1, pack_start=True)
    fig.add_axes(cax)
    if True:
        if 10 * round(_max - _min,1) % 2 == 0 and _max < 10:
            cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks=np.round(np.arange(_min,_max+0.1,0.2),1),orientation="horizontal")
        elif 10 * round(_max - _min,1) % 3 == 0 and _max <10:
            cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks=np.round(np.arange(_min,_max+0.1,0.3),1),orientation="horizontal")
        elif int(_max-_min) > 20:
            cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks=np.round(np.arange(_min,_max+0.1,20),1),orientation="horizontal")
        else:
            cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks=np.round(np.arange(_min,_max+0.1,0.1),1),orientation="horizontal")
    if False:
        #cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks = np.int(np.linspace(_min,_max,3)),orientation="horizontal")
        cbar = fig.colorbar(image,shrink=0.5,aspect=25,orientation="horizontal")
    
    #cbar.ax.locator_params(nbins=tick_numbers)
    cbar.ax.locator_params()
    cbar.set_label(unit,size=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

def diy_colorbar_C(_min,_max,unit=r'log$\tau$',fontsize=12,tick_numbers=8):
    #tick_numbers change how many ticks in colorbar
    #size decide the distance between map and colorbar;shrink change the length of colorbar;aspect change the height of colorbar
    hp.graticule(30,coord="C",local=True)
    #fontsize = 'medium'
    hp.projtext(0., 0, r'$0\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(30., 0., r'$30\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(60., 0., r'$60\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(90., 0., r'$90\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(120., 0., r'$120\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(150, 0., r'$150\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(180., 0., r'$180\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(210., 0, r'$210\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(240., 0, r'$240\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(270., 0, r'$270\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(300., 0, r'$300\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(330., 0, r'$330\degree$', lonlat=True, coord='C',color = 'w',fontsize = fontsize)
    hp.projtext(0., 30., r'$30\degree$', lonlat=True, coord='C',color = 'black',fontsize = fontsize)
    hp.projtext(0., 60., r'$60\degree$', lonlat=True, coord='C',color = 'black',fontsize = fontsize)
    hp.projtext(0., -30, r'$-30\degree$', lonlat=True, coord='C',color = 'black',fontsize = fontsize)
    hp.projtext(0., -60, r'$-60\degree$', lonlat=True, coord='C',color = 'black',fontsize = fontsize)



    fig = plt.gcf()
    ax = plt.gca()
    image = ax.get_images()[0]


    #size decide the distance between map and colorbar;shrink change the length of colorbar;aspect change the height of colorbar
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size="20%",pad=0.1, pack_start=True)
    fig.add_axes(cax)
    if True:
        if 10 * round(_max - _min,1) % 2 == 0 and _max < 10:
            cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks=np.round(np.arange(_min,_max+0.1,0.2),1),orientation="horizontal")
        elif 10 * round(_max - _min,1) % 3 == 0 and _max <10:
            cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks=np.round(np.arange(_min,_max+0.1,0.3),1),orientation="horizontal")
        elif int(_max-_min) > 20:
            cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks=np.round(np.arange(_min,_max+0.1,20),1),orientation="horizontal")
        else:
            cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks=np.round(np.arange(_min,_max+0.1,0.1),1),orientation="horizontal")
    if False:
        #cbar = fig.colorbar(image,shrink=0.5,aspect=25, ticks = np.int(np.linspace(_min,_max,3)),orientation="horizontal")
        cbar = fig.colorbar(image,shrink=0.5,aspect=25,orientation="horizontal")
    
    #cbar.ax.locator_params(nbins=tick_numbers)
    cbar.ax.locator_params()
    cbar.set_label(unit,size=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
