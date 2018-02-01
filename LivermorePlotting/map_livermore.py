import os.path
import configparser
import numpy as np

from map_tools_twh.map_tools_twh import Fig
from map_tools_twh.map_tools_twh import Livermore_prj
from map_tools_twh.map_tools_twh import Livermore_Mapper

from stem_pytools import STEM_parsers as sp
from stem_pytools import domain

from timutils import colormap_nlevs

fname_conf = os.path.join('/', 'nfs', 'pic.es', 'user', 't',
                          'thilton', 'Code', 'UAB_pytools',
                          'STEM_livermore.conf')


def umol_m2_s_to_ppt(umol):
    """convert [OCS] from umol m-2 s-1 to ppt
    """
    return(umol * 1e12)


def get_tower_coords(conf, tower_string):
    """split tower coordinate string to list of floats
    """
    return(list(map(float, conf['towers'][tower_string].split(', '))))


def get_PBL_zlev(conf):
    """calculate planetary boundry layer (PBL) Z level from WRF meteo3d

    Here PBL is defined by a temperature inversion: going up from the
    surface, the first WRF vertical level whose T is higher than the
    next-lower level determines the PBL top.

    ARGS: conf (configparser.ConfigParser instance): parsed
       configuration file, must include section paths, item METEO3D

    """
    T = sp.parse_STEM_var(nc_fname=conf['paths']['METEO3D'],
                          varname='T')
    dT = np.diff(T['data'], axis=1)
    idx = np.argmin(dT < 0, axis=1)
    dom = domain.STEM_Domain(fname_topo=conf['paths']['TOPO'])
    dom.get_STEMZ_height(wrfheight_fname=conf['paths']['HEIGHT'])
    pbl_hgt = np.zeros(idx.shape)
    for t in range(pbl_hgt.shape[0]):
        for x in range(pbl_hgt.shape[1]):
            for y in range(pbl_hgt.shape[2]):
                pbl_hgt[t, x, y] = dom.asl[idx[t, x, y], x, y]
    return(pbl_hgt, idx, dom)


def plot_OCS(conf):
    """plot STEM [OCS] concentration on a map
    """
    (lon, lat, topo) = sp.parse_STEM_coordinates(conf['paths']['TOPO'])
    ocs = sp.parse_STEM_var(nc_fname=conf['paths']['AQOUT'],
                            varname=conf['variables']['AQOUT_OCS'])
    ocs['data'] = umol_m2_s_to_ppt(ocs['data'])
    cmap, norm = colormap_nlevs.setup_colormap(
        vmin=float(conf['plotting']['cbar_min']),
        vmax=float(conf['plotting']['cbar_max']),
        extend='neither',
        nlevs=10)

    for this_t in range(len(ocs['t'])):
        print('plotting time step {:02d}'.format(this_t))
        fig = Fig(figsize=(8, 8))
        ax = fig.add_subplot(111,
                             projection=Livermore_prj())

        this_map = Livermore_Mapper(ax=ax)
        cm = this_map.pcolormesh(lon, lat, ocs['data'][this_t, 0, ...],
                                 cmap=cmap, norm=norm)
        this_map.colorbar(cm, label_str='OCS (ppt)')
        for this_tower in conf['towers']:
            tower_lat, tower_lon = get_tower_coords(conf, this_tower)
            this_map.scatter(tower_lon, tower_lat)
        ax.set_title(ocs['t'][this_t].strftime('%d %b %Y %H:%M'))
        fig.savefig(fname='livermore_map{:03d}.png'.format(this_t))


if __name__ == "__main__":

    conf = configparser.ConfigParser()
    conf.read(fname_conf)

    pbl_hgt, idx, dom = get_PBL_zlev(conf)
