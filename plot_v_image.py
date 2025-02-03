#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.cm as cm
from matplotlib import colors as mpl_colors
from vel_col_cmap import vel_col
from pwr_col_cmap import pwr_col
from median_filter import medianFilter
from date_strings import cnv_datetimestr_dtlist
#from mag_continents import magContinents

from subsolar import subSolarPoint

import argparse

import aacgm

import cartopy.crs as ccrs
import cartopy.feature as cfeature


parser = argparse.ArgumentParser(description='velocity imaging argument parser')
# parser.add_argument('site_code', type=str, nargs=1,
#                     help='station code')
parser.add_argument('fname',type=str,nargs=1)

parser.add_argument("--param", required=False, type=str, default='velocity',
                    help="No parameter specified. Plotting velocity")

parser.add_argument("--pthresh", required=False, type=float, default=0.001)
parser.add_argument('--med_filt', required=False, dest='med_filt', action='store_true', default=False)
parser.add_argument('--neighbors', required=False, type=int, default=1)
parser.add_argument('--level', required=False, type=int, default=1)

parser.add_argument('--min_lat', dest='in_min_lat', required=False, type=float, default=None)
parser.add_argument('--min_lon', dest='in_min_lon', required=False, type=float, default=None)
parser.add_argument('--max_lat', dest='in_max_lat', required=False, type=float, default=None)
parser.add_argument('--max_lon', dest='in_max_lon', required=False, type=float, default=None)

parser.add_argument("--max_range",required=False, type=int, default=np.nan)
parser.add_argument('--scale', required=False, type=int, nargs=2, default=[])

parser.add_argument('--png', required=False, dest='png', action='store_true', default=False)
parser.add_argument('--mag', required=False, dest='mag', action='store_true', default=False)

args = parser.parse_args()
# site_code=args.site_code[0]
fname=str(args.fname[0])
param = args.param
scale = args.scale
png = args.png
pthresh=args.pthresh
med_filt=args.med_filt
neighbors=args.neighbors
level=args.level
in_min_lat=args.in_min_lat
in_min_lon=args.in_min_lon
in_max_lat=args.in_max_lat
in_max_lon=args.in_max_lon
max_range=args.max_range
mag=args.mag

badvalue=-99999.0
nrang=int(75)
nang=int(60)


from read_v_image import readVimage
data=readVimage(fname)

     # if keyword_set(no_gs) then begin
     #    q=where(pwr_ar lt 0.5 or abs(vel_ar) lt 30.)
     #    pwr_ar[q]=badvalue
     #    vel_ar[q]=badvalue
     # endif

date=data["date"]
nrang=data["nrang"]
nang=data["nang"]
lat_ar=data["lat_ar"]
lon_ar=data["lon_ar"]
azm_ar=data["azm_ar"]
vel_ar=data["vel_ar"]
pwr_ar=10*np.log10(data["pwr_ar"])
wid_ar=data["wid_ar"]
print(date)

neg_lons=lon_ar<0
lon_ar[neg_lons]+=360


min_lat=np.min(lat_ar)
min_lon=np.min(lon_ar)
max_lat=np.max(lat_ar)
max_lon=np.max(lon_ar)

if in_min_lat != None: min_lat=in_min_lat
if in_min_lon != None: min_lon=in_min_lon
if in_max_lat != None: max_lat=in_max_lat
if in_max_lon != None: max_lon=in_max_lon

if mag:
    aacgm.set_datetime(date.year,date.month,date.day,date.hour,date.minute,date.second)
    noonlon=aacgm.inv_mlt_convert(date.year,date.month,date.day,date.hour,date.minute,date.second,12)
    coords='mag'
else:
    slat,noonlon=subSolarPoint(date)
    coords='geo'



print(min_lat,max_lat,min_lon,max_lon)

if max_lat < 0:
    south=True
else:
    south=False



if scale == []:
    if(param == 'velocity'):
        scale = [-1000, 1000]
        label='V (m/s)'
    elif(param == 'power'):
        scale = [0, 100]
        label='p_l (dB)'
        gs_nogray=True
    elif(param == 'width'):
        scale = [0, 250]
        label='w_l (m/s)'
    elif(param == 'elevation'):
        scale = [0, 50]
        label=r'\alpha (^\circ)'
        gs_nogray=True
    elif(param == 'phi0'):
        scale = [-numpy.pi, numpy.pi]
        label=r'\Phi (^\circ)'
else:
    if(param == 'velocity'):
        label='V (m/s)'
    elif(param == 'power'):
        label='p_l (dB)'
        gs_nogray=True
    elif(param == 'width'):
        label='w_l (m/s)'
    elif(param == 'elevation'):
        label=r'\alpha (^\circ)'
        gs_nogray=True
    elif(param == 'phi0'):
        label=r'\Phi (^\circ)'

    
if(param == 'velocity'):
    from vel_col_cmap import vel_col
    # if plot_gscat:
    #     cmap=vel_col(groundscat=True)
    # else:
    cmap=vel_col()
else:
    from pwr_col_cmap import pwr_col
    cmap=pwr_col()
    
norm = mpl_colors.Normalize(vmin=scale[0], vmax=scale[1])

inds=np.where(pwr_ar < pthresh)

vel_ar[inds]=np.nan
pwr_ar[inds]=np.nan
wid_ar[inds]=np.nan

year=date.year
month=date.month
day=date.day
hour=date.hour
minute=date.minute

c_lon=np.average(lon_ar)
c_lat=np.average(lat_ar)

fig=plt.figure(figsize=(8,8))


if np.abs(max_lon-min_lon)>200:
    min_lon=0
    max_lon=359
    if south:
        central_lat=-90
        central_lon=noonlon
    else:
        central_lat=90
        central_lon=noonlon-180
        

else:
    central_lat=(min_lat+max_lat)/2
    central_lon=(np.max(c_lon)+np.min(c_lon))/2


print("MIN_LAT: ",min_lat)
print("MAX_LAT: ",max_lat)

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
    
fig.clf()    
crs=ccrs.AzimuthalEquidistant(central_longitude=central_lon,
                              central_latitude=central_lat,
                              false_easting=0.0,
                              false_northing=0.0,
                              globe=None)

ax = fig.add_subplot(1,1,1,projection=crs)
ax.set_extent((min_lon,max_lon,min_lat,max_lat),ccrs.PlateCarree())
ax.gridlines()
if np.abs(max_lon-min_lon) > 270:
    ax.set_boundary(circle, transform=ax.transAxes)



#if mag:
#    magContinents(ax,boundary=True,alpha=.6)
#else:
#    ax.add_feature(cfeature.OCEAN,alpha=.6)

ax.add_feature(cfeature.LAND.with_scale('50m'), edgecolor='black',fc='none')
ax.gridlines()


from date_strings import make_date_time_str
d={"yr":year,"mo":month,"dy":day,"hr":hour,"mt":minute,"sc":0}
dt_tm_str=make_date_time_str(d,colons=False,seconds=False)


if param == 'velocity':
    if (med_filt):
        vel_ar=medianFilter(vel_ar,level=level,min_neighbors=neighbors)
    psm = ax.pcolormesh(lon_ar,lat_ar,vel_ar, cmap=cmap, rasterized=True,
                        vmin=scale[0], vmax=scale[1], alpha=.6,
                        transform=ccrs.PlateCarree())
    figname=dt_tm_str+'_vel'+'.png'


if param == 'power':
    if (med_filt):
        pwr_ar=medianFilter(pwr_ar,level=level,min_neighbors=neighbors)

    scale[1]=np.nanmax(pwr_ar)

    print(np.nanmax(pwr_ar))
    
    psm = ax.pcolormesh(lon_ar,lat_ar,pwr_ar, cmap=cmap, rasterized=True,
                        vmin=scale[0], vmax=scale[1], alpha=.8,
                        transform=ccrs.PlateCarree())
    figname=dt_tm_str+'_pwr'+'.png'


if param == 'width':
    if (med_filt):
        pwr_ar=medianFilter(wid_ar,level=level,min_neighbors=neighbors)
        
    psm = ax.pcolormesh(lon_ar,lat_ar,wid_ar, cmap=cmap, rasterized=True,
                        vmin=scale[0], vmax=scale[1], alpha=.8,
                        transform=ccrs.PlateCarree())
    figname=dt_tm_str+'_wid'+'.png'

    
ax.set_title(dt_tm_str)


transform = ccrs.PlateCarree()._as_mpl_transform(ax)
col='k'
lbfsize='large'
ax.gridlines()    
if np.abs(max_lon-min_lon)>200:
    if south:
        ax.annotate("-60",[0,-60],xycoords=transform,xytext=(0,-60),c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("-70",[0,-70],xycoords=transform,xytext=(0,-70),c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("-80",[0,-80],xycoords=transform,xytext=(0,-80),c=col,fontsize=lbfsize,zorder=10,clip_on=True)

        ax.annotate("Noon",[.5,.85],xycoords='axes fraction',xytext=(.5,.97),horizontalalignment='center',c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("Dusk",[.5,.85],xycoords='axes fraction',xytext=(.95,.5),horizontalalignment='center',c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("Midnight",[.5,.85],xycoords='axes fraction',xytext=(.5,.03),horizontalalignment='center',c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("Dawn",[.5,.85],xycoords='axes fraction',xytext=(.05,.5),horizontalalignment='center',c=col,fontsize=lbfsize,zorder=10,clip_on=True)

    else:
        ax.annotate("60",[0,60],xycoords=transform,xytext=(0,60),c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("70",[0,70],xycoords=transform,xytext=(0,70),c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("80",[0,80],xycoords=transform,xytext=(0,80),c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        
        ax.annotate("Noon",[.5,.85],xycoords='axes fraction',xytext=(.5,.97),horizontalalignment='center',c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("Dusk",[.5,.85],xycoords='axes fraction',xytext=(.05,.5),horizontalalignment='center',c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("Midnight",[.5,.85],xycoords='axes fraction',xytext=(.5,.03),horizontalalignment='center',c=col,fontsize=lbfsize,zorder=10,clip_on=True)
        ax.annotate("Dawn",[.5,.85],xycoords='axes fraction',xytext=(.95,.5),horizontalalignment='center',c=col,fontsize=lbfsize,zorder=10,clip_on=True)

cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
posn=ax.get_position()

barx0=posn.x0+posn.width+0.01
barx1=barx0+.015
y0=posn.y0+.01*posn.height
y1=y0+.5*posn.height
x0_label=(barx0+barx1)/2
# print(posn)

norm = mpl_colors.Normalize(vmin=scale[0], vmax=scale[1])

cbar_ax.set_position([barx0,y0, 0.015*posn.width, .5*posn.height])
cbar=fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbar_ax)
cbar_ax.annotate(label,(0,1),xytext=(0,1.1),xycoords="axes fraction")

if png:
    print(figname)
    fig.savefig(figname,bbox_inches='tight',dpi=200)
else:
    plt.show()
