"""regrid Ian Baker's SiB COS plant fluxes to 1-km Barcelona domain

The module's top-level function is bcn_main.

Examples:
     %python sib_regrid.py --help
     %python sib_regrid.py 2010 10

Author: Timothy W. Hilton, UC Santa Cruz, <twhilton@ucsc.edu>

"""

import os
import os.path
import subprocess
import argparse
from argparse import RawTextHelpFormatter

from timutils.io import delete_if_exists
from IOAPIpytools.ioapi_pytools import calculate_regrid_matrix, run_regrid


def create_IOAPI_file(fname_in, fname_out, year, month, fname_griddesc):
    """run fortran program SiB_2_IOAPI.F90 to create Models-3 I/O API file from
    SiB COS plant fluxes.

    compiles and runs SiB_2_IOAPI.F90 using the makefile SiB_2_IOAPI.mk.

    ARGS:
        fname_in (string): full path for the SiB plant flux netCDF file
        fname_out (string): full path for the SiB plant flux Models-3
            I/O API file to create.  fname_out is overwritten if it
            exists.

    SEE ALSO:
        `Models-3 I/O API mtxcalc documentation
        <https://www.cmascenter.org/ioapi/documentation/3.1/html/>`_
    """
    os.environ['RAW_SIB_FILE'] = fname_in
    os.environ['OUTPUT'] = fname_out
    os.environ['GRIDDESC'] = fname_griddesc
    subprocess.call('make -f SiB_2_IOAPI.mk', shell=True)
    delete_if_exists(fname_out)
    subprocess.call('./SiB_2_IOAPI.x {} {}'.format(year, month), shell=True)



def create_fCOS_IOAPI(fname_regridded,
                      fname_fCOS_out,
                      in_varname,
                      units_factor,
                      new_units,
                      new_desc):
    """perform units conversions necessary to use regridded SiB COS
    fplant data in STEM.

    The regridded SiB fluxes contain five fluxes: GPP, RE, fCOS
    (mechanistic), fCOS (soil) and fCOS (total).  All units are umol
    m-2 gridcell-1.  This function applies a specified units
    conversion factor to fCOS (mechanistic) and writes a new I/O API
    file.  The variable name is changed to cos to meet current STEM
    code's expectations.

    create_fCOS_IOAPI is a python wrapper for the `m3combo
    <https://www.cmascenter.org/ioapi/documentation/3.1/html/M3COMBO.html>`_
    and `m3edhdr
    <https://www.cmascenter.org/ioapi/documentation/3.1/html/M3EDHDR.html>`_
    tools from the Models-3 I/O API.

    ARGS:
        fname_regridded (string): full path to the regridded fCOS
            Models-3 I/O API file containing GPP, RE, and three fCOS products.
        fname_fCOS_out (string): full path to the fCOS
            Models-3 I/O API file to create.  Will be overwritten if
            it already exists.
        in_varname (string): name of the I/O API variable to extract
            and apply units conversion to.
        units_factor (real): units conversion factor to apply to all
            values of in_varname
        new_units (string): string to place in the 'units' metadata
            field of the variable 'cos' in fname_COS_out.  80
            character maximum.
        new_desc (string): string to place in the vdesc metadata field
            of the variable 'cos' in the file fname_fCOS_out.

    SEE ALSO:
        Models-3 I/O API `documentation
        <https://www.cmascenter.org/ioapi/documentation/3.1/html/>`_
    """

    os.environ['COMBO_FILES'] = 'infile,infile'
    os.environ['COMBO_VBLES'] = 'cos_no_adj,cos'
    os.environ['COMBO_UNITS'] = '{},{}'.format(new_units, new_units)
    os.environ['cos_VBLES'] = '{}'.format(in_varname)
    os.environ['cos_COEFS'] = '{}'.format(units_factor)
    os.environ['cos_OFFSET'] = '0.0'
    os.environ['cos_FILES'] = 'infile'
    os.environ['cos_no_adj_VBLES'] = '{}'.format(in_varname)
    os.environ['cos_no_adj_COEFS'] = '{}'.format(units_factor)
    os.environ['cos_no_adj_OFFSET'] = '0.0'
    os.environ['cos_no_adj_FILES'] = 'infile'
    os.environ['outfile'] = 'OUTPUT_FILE'
    # the logical name of the output file is COMBO_3D -- this differs
    # from m3combo documentation
    os.environ['COMBO_3D'] = fname_fCOS_out
    os.environ['infile'] = fname_regridded

    delete_if_exists(fname_fCOS_out)
    subprocess.call('m3combo << DONE\n'
                    '\n'
                    '\n'
                    '\n'
                    '\n'
                    '\n'
                    '\n'
                    'DONE\n',
                    shell=True)

    #  the m3combo automatically-generated variable description is not
    #  very useful -- replace it with a better one.
    os.environ['INFILE'] = fname_fCOS_out
    subprocess.call('m3edhdr << DONE\n'
                    '\n'   # accept default infile logical name
                    '5\n'  # option 5: edit variable units, description
                    '\n'   # accept current variable name
                    '\n'   # accept current variable units
                    '{new_desc}\n'  # new vdesc
                    '\n'   # accept current variable name
                    '\n'   # accept current variable units
                    '{new_desc}\n'  # new vdesc
                    '10\n'  # option 10: quit
                    'DONE\n'.format(new_desc=new_desc),
                    shell=True)


def bcn_main(year, month, skip_regrid=False, skip_post_process=False):
    """main function for module sib_regrid

    creates two new I/O API files from 'native' 2.5 by 2.0 degree SiB
    COS plant flux netCDF files.  The two files contain:

    SiB 'mechanistic' plant fCOS, calculated within SiB using SiB's
    mechanistic canopy model.

    SiB 'calculated' plant fCOS, calculated from SiB GPP using a
    leaf-scale relative uptake (LRU) of 1.61 and a COS/CO2 ratio of
    1.1.

    The 'native' SiB files are monthly; thus bcn_main takes the
        year and month of the SiB file to process as arguments.

    ARGS:
        year (integer): year of the raw SiB file to process.
        monthy (integer): month of the SiB file to process.
        skip_regrid (True | {False}): if True, skip the regridding
            portion and go directly to postprocessing the regridded fluxes
            into fCOS files for STEM.
        skip_post_process (True | {False}): if True, do not
            postprocess the regridded fluxes into fCOS files for STEM.
    """
    fname_griddesc = '/nfs/pic.es/user/t/thilton/Code/Barcelona/GRIDDESC_BCN'
    fname_mattxt = './SiB_to_1kmBCN_mattxt.txt'
    fname_matrix = './SiB_to_1kmBCN_matrix.txt'
    fname_SiB_raw = os.path.join(
        '/software', 'co2flux', 'SurfaceFluxData', 'SIB',
        'flux_hourly_{}{:02}p001.nc'.format(year, month))
    fname_ioapi = 'SiB_{}{:02}_1.25x1.0_IOAPI.nc'.format(
        year, month)
    fname_regridded = 'SiB_{}{:02}_1kmBCN.nc'.format(
        year, month)
    fname_COS_SiBmech = 'SiB_{}{:02}_1kmBCN_mechanistic.nc'.format(
        year, month)
    fname_COS_SiBcalc = 'SiB_{}{:02}_1kmBCN_calculated.nc'.format(
        year, month)
    fname_GPP = 'SiB_{}{:02}_1kmBCN_GPP.nc'.format(
        year, month)

    print(fname_ioapi)
    print(fname_regridded)
    print(fname_COS_SiBmech)
    print(fname_COS_SiBcalc)
    print(fname_GPP)

    if skip_regrid is False:
        # --------------------------------------------------
        # Calculate grid-to-grid transform matrix.  Use 1000 by 1000
        # subsampling; this results in 1e4 samples per 1.25 degree by 1
        # degree area.

        # calculate_regrid_matrix(fname_griddesc, fname_matrix, fname_mattxt,
        #                         col_refinement=1000, row_refinement=1000)
        # --------------------------------------------------

        create_IOAPI_file(fname_in=fname_SiB_raw,
                          fname_out=fname_ioapi,
                          year=year,
                          month=month,
                          fname_griddesc=fname_griddesc)

    #     run_regrid(fname_ioapi, fname_regridded, fname_matrix, fname_mattxt)

    # if skip_post_process is False:
    #     m2_per_gridcell = 6e4 * 6e4  # STEM grid cells are 60 km by 60 km
    #     mols_per_umol = 1e-6  # moles per micromole
    #     create_fCOS_IOAPI(fname_regridded,
    #                       fname_COS_SiBmech,
    #                       in_varname='OCS_gpp',
    #                       units_factor=(-1.0 * (1 / m2_per_gridcell) *
    #                                     mols_per_umol),
    #                       new_units='mol m-2 s-1',
    #                       new_desc=('SiB mechanistic COS plant flux regridded '
    #                                 'to 124x124 N American grid for STEM'))

    #     LRU = 1.61
    #     COS_CO2_ratio = 1.1
    #     umol_per_pmol = 1e-6
    #     # SiB GPP is umol m-2 s-1.  To calculate fCOS in mol m-2 s-1 I
    #     # need mols_per_umol to take micromoles CO2 GPP to moles CO2
    #     # GPP *and* umol_per_pmol to cancel out the pmol COS per umol
    #     # CO2 in the units of LRU.
    #     create_fCOS_IOAPI(fname_regridded,
    #                       fname_COS_SiBcalc,
    #                       in_varname='GPP',
    #                       units_factor=(-1.0 * (1 / m2_per_gridcell) *
    #                                     mols_per_umol * umol_per_pmol *
    #                                     LRU * COS_CO2_ratio),
    #                       new_units='mol m-2 s-1',
    #                       new_desc=('COS Fplant from SiB CO2 GPP, '
    #                                 '124x124 N American '
    #                                 'grid, LRU=1.61, [COS]/[CO2]=1.1'))

    #     ioapi_pytools.ioapi_const_multiply(
    #         fname_regridded,
    #         fname_GPP,
    #         in_varname='GPP',
    #         out_varname='GPP',
    #         constant_factor=((1 / m2_per_gridcell) *
    #                          mols_per_umol),
    #         new_units='mol m-2 s-1',
    #         new_desc=('CO2 GPP from SiB, '
    #                   '124x124 N American grid'))


        # delete_if_exists(fname_ioapi)
        # delete_if_exists(fname_regridded)


if __name__ == "__main__":
    """ parse command line arguments and pass them to bcn_main.

    Examples:
        %python sib_regrid.py --help
        %python sib_regrid.py 2010 10
    """
    # --------------------------------------------------
    # parse arguments
    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description=(
            ("regrid the SiB CO2 and COS GPP to the 124x124 60km North "
             "American STEM grid")))
    parser.add_argument('year', metavar='year', type=int,
                        help='the year of the SiB data ')
    parser.add_argument('month', metavar='month', type=int,
                        choices=range(1, 13),
                        help='the month of the SiB data ')
    parser.add_argument('--diagnostics',
                        dest='run_diagnostics',
                        action='store_true',
                        help=('if set, some diagnostic plots are created'
                              'after the calculations.'))
    parser.add_argument('--skip_regrid',
                        dest='skip_regrid',
                        action='store_true',
                        help=('if set, skip the regridding portion.'))
    parser.add_argument('--skip_postprocess',
                        dest='skip_post_process',
                        action='store_true',
                        help=('if set, skip the post processing portion.'))
    args = parser.parse_args()
    # --------------------------------------------------

    bcn_main(args.year,
                   args.month,
                   args.skip_regrid,
                   args.skip_post_process)

    if args.run_diagnostics:
        print('diagnostics not available on workaround server')
