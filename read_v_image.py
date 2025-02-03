import numpy as np
import datetime as dt
from date_strings import cnv_datetimestr_dtlist

def readVimage(fname):

    badvalue=-99999.0

    strs=fname.split()
    date=strs[0]

    f=open(fname)
    date=f.readline().split()
    [nrang,nang]=f.readline().split()
    site_code=f.readline().split()

    nang=int(nang)
    nrang=int(nrang)
    
    rdar=[]
    for line in f:
        rdar.append(line)
        
    vel_ar=np.zeros((nrang,nang))
    pwr_ar=np.zeros((nrang,nang))
    wid_ar=np.zeros((nrang,nang))
    azm_ar=np.zeros((nrang,nang))
    count=np.zeros((nrang,nang))
    lat_ar=np.zeros((nrang,nang))
    lon_ar=np.zeros((nrang,nang))
    
    for j in range(len(rdar)):
        line=rdar[j].split()
        jr=int(line[0])
        ja=int(line[1])
        if jr < nrang:
            if float(line[2]) == badvalue:
                vel_ar[jr,ja]=np.nan
                pwr_ar[jr,ja]=np.nan
                wid_ar[jr,ja]=np.nan
            else:
                vel_ar[jr,ja]=float(line[2])
                pwr_ar[jr,ja]=float(line[3])
                wid_ar[jr,ja]=float(line[4])
            
            lat_ar[jr,ja]=float(line[5])
            lon_ar[jr,ja]=float(line[6])
            azm_ar[jr,ja]=float(line[7])

    badvalue=np.nan

    vel_ar[:,0]=badvalue
    vel_ar[:,1]=badvalue
    vel_ar[:,nang-1]=badvalue
    vel_ar[:,nang-2]=badvalue
    pwr_ar[:,0]=badvalue
    pwr_ar[:,1]=badvalue
    pwr_ar[:,nang-1]=badvalue
    pwr_ar[:,nang-2]=badvalue
    wid_ar[:,0]=badvalue
    wid_ar[:,1]=badvalue
    wid_ar[:,nang-1]=badvalue
    wid_ar[:,nang-2]=badvalue
    
    [year, month, day, hour, minute]=cnv_datetimestr_dtlist(date[0])
    date=dt.datetime(year,month,day,hour,minute)

    return({"date":date, "nrang":nrang, "nang":nang,
            "vel_ar":vel_ar, "pwr_ar":pwr_ar, "wid_ar":wid_ar,
            "lat_ar":lat_ar, "lon_ar":lon_ar, "azm_ar":azm_ar})
