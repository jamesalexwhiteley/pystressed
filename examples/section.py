import matplotlib.pyplot as plt
from pystressed.models import SectionForce
from pystressed.SLS import plot_magnel, optimize_magnel, optimize_and_plot_magnel

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams.update({'font.size': 12})

# Author: James Whiteley (github.com/jamesalexwhiteley) 

# transfer 
A = 0.16e6 # mm2 
fc = 20 # N/mm2 
ft = 0 # N/mm2
Ztop = 34.4e6 # mm3
Zbot = 29.3e6 # mm3
Mmax = 430e6 # Nmm
Mmin = 108e6 # Nmm
losses = 1.0
ebounds = [0, 500] # mm
transfer = SectionForce(A=A, fc=fc, ft=ft, Ztop=Ztop, Zbot=Zbot, Mmax=Mmax, Mmin=Mmin, losses=losses)

# service 
A = 0.16e6 # mm2
fc = 20 # N/mm2 
ft = 0 # N/mm2
Ztop = 34.4e6 # mm3 
Zbot = 29.3e6 # mm3     
Mmax = 530e6 # Nmm 
Mmin = 158e6 # Nmm 
losses = 1.0
ebounds = [0, 500] # mm 
service = SectionForce(A=A, fc=fc, ft=ft, Ztop=Ztop, Zbot=Zbot, Mmax=Mmax, Mmin=Mmin, losses=losses) 

# design at SLS 
# optimize_magnel(transfer=transfer, service=service, ebounds=ebounds, mode='min', output=True) 
# plot_magnel(transfer=transfer, service=service, line_ext=2.2) 
optimize_and_plot_magnel(transfer, service, ebounds, output=False, line_ext=2.0)
