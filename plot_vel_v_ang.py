#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.cm as cm
from matplotlib import colors as mpl_colors
import os
import sys
from vel_col_cmap import vel_col
from pwr_col_cmap import pwr_col
import PyDarn
import glob

import argparse
parser = argparse.ArgumentParser(description='velocity imaging argument parser')
parser.add_argument('filename',type=str,nargs=1,
                    help='date and time (yyyymmddhhmm)')

parser.add_argument('site_code', type=str, nargs=1,
                    help='station code')

parser.add_argument("--param", required=False, type=str, default='velocity',
                    help="No parameter specified. Plotting velocity")

parser.add_argument("--pthresh", required=False, type=float, default=2)
parser.add_argument("--med_filt",required=False, type=int, default=np.nan)
parser.add_argument("--max_range",required=False, type=int, default=np.nan)
parser.add_argument("--scale", required=False,type=float, default=None)
parser.add_argument('--png', required=False, dest='png', action='store_true', default=False)
parser.add_argument('--neighbors', required=False, type=int, default=1)

parser.add_argument('--min_lat', dest='in_min_lat', required=False, type=float, default=None)
parser.add_argument('--min_lon', dest='in_min_lon', required=False, type=float, default=None)
parser.add_argument('--max_lat', dest='in_max_lat', required=False, type=float, default=None)
parser.add_argument('--max_lon', dest='in_max_lon', required=False, type=float, default=None)


args = parser.parse_args()
site_code=args.site_code[0]
filename=str(args.filename[0])
param = args.param
scale = args.scale
png = args.png
pthresh=args.pthresh
med_filt=args.med_filt
max_range=args.max_range
neighbors=args.neighbors
in_min_lat=args.in_min_lat
in_min_lon=args.in_min_lon
in_max_lat=args.in_max_lat
in_max_lon=args.in_max_lon

print('param: ',param)
print('site: ',site_code)
print('filename: ',filename)

vel_ar=None

print('filename: ',filename)
    
basename = os.path.basename(filename)

strs=basename.split(".")
date=strs[0]
beam=strs[1]

f=open(filename)
date=f.readline().split()
[nrang, nang, smsep]=f.readline().split()
[bm, tfreq, mpinc]=f.readline().split()
    
if np.isnan(max_range):        
    nrang=int(nrang)
else:
    nrang=int(max_range)
        
nang=int(nang)
if not isinstance(vel_ar,np.ndarray):
    vel_ar=np.zeros((nrang,nang))
    pwr_ar=np.zeros((nrang,nang))
    wid_ar=np.zeros((nrang,nang))
    azm_ar=np.zeros((nrang,nang))
    count=np.zeros((nrang,nang))
        
for line in f:
        
    strs=line.split()
    jr=int(strs[0])
    ja=int(strs[1])
    if jr < nrang:
        azm_ar[jr,ja]+=float(strs[2])
        pwr_ar[jr,ja]+=float(strs[3])
        wid_ar[jr,ja]+=float(strs[4])
        vel_ar[jr,ja]+=float(strs[5])
        count[jr,ja]+=1
        
        
rsep=int(float(smsep)*0.15)
nrang=int(nrang)
bmsep=float(1)
nang=int(nang)
coords='geo'
lon_grid=np.zeros((nrang,nang))
lat_grid=np.zeros((nrang,nang))

    
for jr in range(int(nrang)):
    for ja in range(int(nang)):
        if count[jr,ja] > 0:
            azm_ar[jr,ja]/=count[jr,ja]
            pwr_ar[jr,ja]/=count[jr,ja]
            wid_ar[jr,ja]/=count[jr,ja]
            vel_ar[jr,ja]/=count[jr,ja]


print("pthresh: ",pthresh)
for jr in range(1,nrang-1):
    for ja in range(1,nang-1):
        # if (pwr_ar[jr-1,ja] < pthresh and pwr_ar[jr+1,ja] < pthresh
        # and pwr_ar[jr,ja+1] < pthresh and pwr_ar[jr,ja+1] < pthresh):
        if( pwr_ar[jr,ja]>pthresh):
            # print(jr,ja,pwr_ar[jr,ja],count[jr,ja])
            pwr_ar[jr,ja]=pwr_ar[jr,ja]
            wid_ar[jr,ja]=wid_ar[jr,ja]
            vel_ar[jr,ja]=vel_ar[jr,ja]
        else:
            pwr_ar[jr,ja]=0
            wid_ar[jr,ja]=0
            vel_ar[jr,ja]=0
            
for jr in range(nrang):
    pwr_ar[jr,0]=0
    wid_ar[jr,0]=0
    vel_ar[jr,0]=0
    pwr_ar[jr,1]=0
    wid_ar[jr,1]=0
    vel_ar[jr,1]=0
    pwr_ar[jr,2]=0
    wid_ar[jr,2]=0
    vel_ar[jr,2]=0    
    pwr_ar[jr,nang-1]=0
    wid_ar[jr,nang-1]=0
    vel_ar[jr,nang-1]=0


vel_ar=np.where(vel_ar==0,np.nan,vel_ar)
if not np.isnan(med_filt):
    from median_filt import medianFilter
    vel_ar=medianFilter(vel_ar,level=med_filt,min_neighbors=neighbors)


vel_plt = np.ma.masked_where(((np.isfinite(vel_ar) != 1) | (np.isfinite(pwr_ar) != 1)), vel_ar)
pwr_plt = np.ma.masked_where(((np.isfinite(vel_ar) != 1) | (np.isfinite(pwr_ar) != 1)), pwr_ar)
wid_plt = np.ma.masked_where(((np.isfinite(vel_ar) != 1) | (np.isfinite(pwr_ar) != 1)), wid_ar)

print(np.max(pwr_plt))

if param == 'velocity':
    cmap=vel_col()
    if scale == None:
        mxval=1000
    else:
        mxval=scale
    mnval=-mxval

if param == 'power':
    cmap=pwr_col()
    inds=np.where(pwr_plt != 0)
    pwr_plt[inds]=20*np.log10(pwr_plt[inds])
    if scale == None:
        mxval=np.max(pwr_plt)
    else:
        mxval=scale
        
    mnval=mxval-30

if param == 'width':
    cmap=pwr_col()
    if scale == None:
        mxval=500
    else:
        mxval=scale
        
    mnval=0
    

import datetime as dt
from PyDarn.radar import site, radFov
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from date_strings import make_date_str, make_date_time_str

year=int(date[0])
month=int(date[1])
day=int(date[2])
hour=int(date[3])
minute=int(date[4])

date=dt.datetime(*map(int,date))

siteS = site(code=site_code, dt=date)
fov = radFov.fov(site=siteS, rsep=rsep, ngates=nrang+1,
                              nbeams=nang, bmsep=bmsep, coords=coords,
                              date_time=date)

neglons=fov.lonFull<0
fov.lonFull[neglons]+=360.

c_lon=np.average(fov.lonFull)
c_lat=np.average(fov.latFull)
min_lat=np.min(fov.latFull)
max_lat=np.max(fov.latFull)
min_lon=np.min(fov.lonFull)
max_lon=np.max(fov.lonFull)

if in_min_lat != None: min_lat=in_min_lat
if in_min_lon != None: min_lon=in_min_lon
if in_max_lat != None: max_lat=in_max_lat
if in_max_lon != None: max_lon=in_max_lon


print(min_lon,max_lon)
    
for jr in range(nrang):
    for ja in range(nang):
        lon_grid[jr,ja]=fov.lonFull[ja+1,jr+1]
        lat_grid[jr,ja]=fov.latFull[ja+1,jr+1]

fig=plt.figure(figsize=(8,8))


theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

fig.clf()    
crs=ccrs.AzimuthalEquidistant(central_longitude=c_lon,
                              central_latitude=c_lat,
                              false_easting=0.0,
                              false_northing=0.0,
                              globe=None)

ax = fig.add_subplot(1,1,1,projection=crs)
ax.set_extent((min_lon,max_lon,min_lat,max_lat),ccrs.PlateCarree())
ax.gridlines()
if np.abs(max_lon-min_lon) > 270:
    ax.set_boundary(circle, transform=ax.transAxes)


mag=False
if mag:
    magContinents(ax,boundary=True,alpha=.6)
else:
    ax.add_feature(cfeature.LAND.with_scale('50m'), edgecolor='black',fc='none')
    ax.gridlines()

ax.set_title(date)

from date_strings import make_date_time_str
d={"yr":year,"mo":month,"dy":day,"hr":hour,"mt":minute,"sc":0}
dt_tm_str=make_date_time_str(d,colons=False,seconds=False)

if param == 'velocity':
    psm = ax.pcolormesh(lon_grid,lat_grid,vel_plt, cmap=cmap, rasterized=True,
                        vmin=mnval, vmax=mxval, alpha=.6,
                        transform=ccrs.PlateCarree())
    figname=dt_tm_str+'_'+beam+'_vel.png'
    label='V (m/s)'

if param == 'power':
    psm = ax.pcolormesh(lon_grid,lat_grid,pwr_plt, cmap=cmap, rasterized=True,
                        vmin=mnval, vmax=mxval, alpha=.8,
                        transform=ccrs.PlateCarree())
    figname=dt_tm_str+'_'+beam+'_pwr.png'
    label='p_l (dB)'

if param == 'width':
    psm = ax.pcolormesh(lon_grid,lat_grid,wid_plt, cmap=cmap, rasterized=True,
                        vmin=mnval, vmax=mxval, alpha=.8,
                        transform=ccrs.PlateCarree())
    figname=dt_tm_str+'_'+beam+'_wid.png'
    label='w_l (m/s)'




cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
posn=ax.get_position()

barx0=posn.x0+posn.width+0.01
barx1=barx0+.015
y0=posn.y0+.01*posn.height
y1=y0+.5*posn.height
x0_label=(barx0+barx1)/2
# print(posn)

norm = mpl_colors.Normalize(vmin=mnval, vmax=mxval)

cbar_ax.set_position([barx0,y0, 0.015*posn.width, .5*posn.height])
cbar=fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbar_ax)
cbar_ax.annotate(label,(0,1),xytext=(0,1.1),xycoords="axes fraction")

ax.annotate( "Transmit Beam: "+str(bm),(0.1,0.95),xytext=(0.8,1.025),xycoords="axes fraction")

if png:
    print(figname)
    plt.savefig(figname,bbox_inches='tight',dpi=200)
else:
    plt.show()


