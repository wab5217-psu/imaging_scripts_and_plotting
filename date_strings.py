import datetime as dt

def make_date_str(d):
    istrs=["00","01","02","03","04","05","06","07","08","09"]

    if d.mo < 10: mostr=istrs[d.mo]
    else: mostr=str(d.mo)

    if d.dy < 10: dystr=istrs[d.dy]
    else: dystr=str(d.dy)

    dtstr=str(d.yr)+mostr+dystr
    return dtstr

def make_date_time_str(d,colons=True,seconds=True):
    istrs=["00","01","02","03","04","05","06","07","08","09"]

    if isinstance(d,dict):
        if d["mo"] < 10: mostr=istrs[d["mo"]]
        else: mostr=str(d["mo"])
        
        if d["dy"] < 10: dystr=istrs[d["dy"]]
        else: dystr=str(d["dy"])
        
        if d["hr"] < 10: hrstr=istrs[d["hr"]]
        else: hrstr=str(d["hr"])
        
        if d["mt"] < 10: mtstr=istrs[d["mt"]]
        else: mtstr=str(d["mt"])
        
        if d["sc"] < 10: scstr=istrs[d["sc"]]
        else: scstr=str(d["sc"])

        if colons:
            dtstr=str(d["yr"])+mostr+dystr+"-"+hrstr+":"+mtstr
        else:
            dtstr=str(d["yr"])+mostr+dystr+"-"+hrstr+mtstr

        if seconds:
            if colons:
                dtstr+=":"+scstr
            else:
                dtstr+=scstr        
    else:
        if d.mo < 10: mostr=istrs[d.mo]
        else: mostr=str(d.mo)
        
        if d.dy < 10: dystr=istrs[d.dy]
        else: dystr=str(d.dy)
        
        if d.hr < 10: hrstr=istrs[d.hr]
        else: hrstr=str(d.hr)
        
        if d.mt < 10: mtstr=istrs[d.mt]
        else: mtstr=str(d.mt)
        
        if d.sc < 10: scstr=istrs[d.sc]
        else: scstr=str(d.sc)
        
        dtstr=str(d.yr)+mostr+dystr+"-"+hrstr+":"+mtstr+":"+scstr
        
    return dtstr

def cnv_datetimestr_dtlist(d):
    d1=int(d)
    yr=int(d1/100000000)

    d1=d1-100000000*yr
    mo=int(d1/1000000)

    d1=d1-mo*1000000
    dy=int(d1/10000)

    d1=d1-10000*dy
    hr=int(d1/100)

    d1=d1-100*hr
    mt=d1

    return([yr,mo,dy,hr,mt])

def cnv_datetimestr_datetime(d):
    d1=int(d)
    yr=int(d1/100000000)

    d1=d1-100000000*yr
    mo=int(d1/1000000)

    d1=d1-mo*1000000
    dy=int(d1/10000)

    d1=d1-10000*dy
    hr=int(d1/100)

    d1=d1-100*hr
    mt=d1

    
    return(dt.datetime(yr,mo,dy,hr,mt,0))
