import numpy as np

def medianFilter(in_ar,level=0,min_neighbors=1,filter_value=0,param=None):
    """
    returns an array of same shape as input
    for level=0 removes any point that have less then min_neighbors adjacent points
    for level=1 replaces points by average of surrounding points
    """

    if isinstance(in_ar,np.ndarray):
        dim0=len(in_ar)
        dim1=len(in_ar[0])
    else:

        [dim0,dim1]=in_ar.shape
        
    hold_ar=np.zeros([dim0,dim1])
    hold_ar.fill(np.nan)

    if np.isnan(filter_value):
        neighbor_ar=np.where(np.isnan(in_ar)==1,0,1)
    else:
        neighbor_ar=np.where(in_ar==filter_value,0,1)        


    mask=np.zeros([3,3])
    mask[0,:]=1
    mask[2,:]=1
    mask[1,0]=1
    mask[1,2]=1

    for i0 in range(1,dim0-1):
        for i1 in range(1,dim1-1):

            if level == 0:
                ncount=np.sum(mask*neighbor_ar[i0-1:i0+2,i1-1:i1+2])
                if ncount >= min_neighbors:
                    hold_ar[i0,i1]=in_ar[i0,i1]
                # else:
                #     hold_ar[i0,i1]=np.nan

            if level == 1:
                count=0
                value=0
                for i2 in range(3):
                    for i3 in range(3):
                        if (np.isfinite(in_ar[i0-1+i2,i1-1+i3])):
                            count+=1
                            value+=in_ar[i0-1+i2,i1-1+i3]
                            
                if ( np.isfinite(in_ar[i0,i1]) )and( count > min_neighbors+1 ):
                    hold_ar[i0,i1]=value/count
                # else:
                #     hold_ar[i0,i1]=np.nan
                    
    return(hold_ar)
