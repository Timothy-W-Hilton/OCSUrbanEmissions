"""make a dummy PBL file according to Medeiros et al (2005), Fig. 8

REFERENCES
Medeiros, B., A. Hall, and B. Stevens, 2005: What Controls the Mean
Depth of the PBL?. J. Climate, 18, 3157–3172,
https://doi.org/10.1175/JCLI3417.1
"""

import numpy as np

Nt = 81
ONELAYER = 1
Nx = 50
Ny = 50
WRF_TSTEP = 6  # WRF time step, hours
UTC_PST_OFFSET = 8  # offset from UTC to Pacific Standard Time, hours

# PBL height follows Medeiros et al (2005), Fig. 8.  PBL height is 100
# m at 04:00 PST, 1100 m at 10:00 PST, 1600 m at 16:00 PST, and 100 m
# at 22:00 PST.
PBL_HEIGHT = {4: 100, 10: 1100, 16: 1600, 22: 100}

if __name__ == "__main__":
    pbl = np.zeros((Nt, ONELAYER, Nx, Ny))
    for t in range(Nt):
        t_utc = (t * WRF_TSTEP) % 24  # simulation hour of day, UTC
        t_pst = (t_utc - UTC_PST_OFFSET) % 24
        print('t: {}, tUTC: {}, tPST: {}, PBL: {} m'.format(
            t, t_utc, t_pst, PBL_HEIGHT[t_pst]))
        pbl[t, 0, :, :] = PBL_HEIGHT[t_pst]


    # write pbl to ASCII file for M3FAKE
    # (https://www.cmascenter.org/ioapi/documentation/all_versions/html/M3FAKE.html)
    np.savetxt('dummy_pbl.csv',
               pbl.reshape((-1, Nx)),
               delimiter=',', fmt='%1.2f')
