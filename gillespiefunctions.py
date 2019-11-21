import random
import numpy as np

def gillespie (final_time):
    ko = 0.2
    k1 = 0.01
    t=0
    m=1
    time_list = [0]
    mRNA_list = [0]
    while t < final_time:
        a0=ko
        a1=k1*m
        t1=np.random.exponential(np.divide(1,a0))
        t2=np.random.exponential(np.divide(1,a1))
        if t1<t2:  
            m+=1
            t+=t1
        if t1>t2:
            if m>0:   
                m+=-1
                t+= t2
        mRNA_list.append(m)
        time_list.append(t)
    return (time_list, mRNA_list)


def celldivision_gillespie (final_time):
    ko = 0.2
    k1 = 0.01
    t=0
    m=1
    time_list = [0]
    mRNA_list = [0]
    mRNA_division = [0]
    division = 1
    
    while t < final_time:
        a0=ko
        a1=k1*m
        t1=np.random.exponential(np.divide(1,a0))
        t2=np.random.exponential(np.divide(1,a1))
        if t1<t2:  
            m+=1
            t+=t1
        if t1>t2:
            if m>0:   
                m+=-1
                t+= t2
        if t>=1200*division:
            m/=2
            division+=1
            mRNA_division.append(m)
        mRNA_list.append(m)
        time_list.append(t)
    return (time_list, mRNA_list, mRNA_division)


def celldivision_gillespie_total (final_time):
    ko = 0.2
    k1 = 0.01
    t=0
    m=1
    time_list = [0]
    mRNA_list = [0]
    mRNA_division = [0]
    generation = 1
    m_tot=0
    m_tot_l=[0]
    while t < final_time:
        a0=ko
        a1=k1*m
        t1=np.random.exponential(np.divide(1,a0))
        t2=np.random.exponential(np.divide(1,a1))
        if t1<t2:  
            m+=1
            t+=t1
        if t1>t2:
            if m>0:   
                m+=-1
                t+= t2
        m_tot=m*2**(generation-1)
        if t>=1200*generation:
            m/=2
            generation+=1
            mRNA_division.append(m)
        m_tot_l.append(m_tot)
        mRNA_list.append(m)
        time_list.append(t)
    return (time_list, mRNA_list, m_tot_l, )


