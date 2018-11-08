#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
==================================================
PROJECT

Sentinel Synergy Study (S^3)

DESCRIPTION

Retrieval System Preprocessor

Implements/Provides all preprocessing steps required to operate the retrieval
prototype system

EXAMPLES

available options are listed by invoking 'rs_pre.py -h' on the command line


EXIT STATUS

should be 0 in case of success, 1 otherwise

AUTHOR

The Inversion Lab, Michael Vossbeck <Michael.Vossbeck@Inversion-Lab.com>

==================================================
"""

#-- general packages
import sys
import os
import shutil
import numpy as np, numpy.ma as ma
import datetime as dt
import netCDF4 as nc4
from collections import OrderedDict
import logging
import f90nml

#-- signature simulator
import signaturesimulator as ss
import signaturesimulator.satellite_geometry as satgeo
import signaturesimulator.state_vector       as sv


ss_dir_path  = os.path.dirname(os.path.realpath(ss.__file__))
dir_path     = os.path.dirname(os.path.realpath(__file__))
ipt_dir_path = os.path.join(os.path.dirname(dir_path),'input')

#--------------------------------------------------
#     l o g g i n g
#
FileLogger   = logging.getLogger( __file__ )
OUT_HDLR = logging.StreamHandler( sys.stdout )
OUT_HDLR.setFormatter( logging.Formatter( '%(asctime)s %(levelname)s::%(funcName)s:: %(message)s') )
OUT_HDLR.setLevel( logging.INFO )
FileLogger.addHandler( OUT_HDLR )
FileLogger.setLevel( logging.INFO )



#--------------------------------------------------
#     g l o b a l s
#
# MISSION_DICT = OrderedDict( [('S1a', 'Sentinel-1a'), ('S1b', 'Sentinel-1b'),
#                              ('S2a', 'Sentinel-2a'), ('S2b', 'Sentinel-2b')] )
#-- IMPL-NOTE::TLE file does not yet provide 'Sentinel-2b'
MISSION_DICT = OrderedDict( [('S1a', 'Sentinel-1a'),   ('S1b', 'Sentinel-1b'),
                             ('S2a', 'Sentinel-2a'),   ('S2b', 'Sentinel-2b'),
                             ('S1', 'Sentinel-1a/1b'), ('S2', 'Sentinel-2a/2b') ] )


BOUNDS_DICT  = OrderedDict( [('lai_coeff_lobound', 1.e-5),
                             ('lai_coeff_hibound', 1.),
                             ('lai_lobound',       1.e-3),
                             ('lai_hibound',      10.),
                             ('canht_lobound',     1.e-3), # 1mm
                             ('canht_hibound',     4.5),
                             ('sm_lobound',        1.e-3),
                             ('sm_hibound',        0.5)]     )

TIGHT_BOUNDS_DICT  = OrderedDict( [('lai_coeff_lobound', 1.e-5),
                                   ('lai_coeff_hibound', 1.),
                                   ('lai_lobound',       1.e-3),
                                   ('lai_hibound',       8.),
                                   ('canht_lobound',     1.e-3), # 1mm
                                   ('canht_hibound',     3.5),
                                   ('sm_lobound',        0.1),
                                   ('sm_hibound',        0.4)]     )

#-- default time period for the 'pre_synthetic' case
TIME_START_SYN = dt.datetime(2017,5,1)
TIME_END_SYN   = dt.datetime(2017,5,31,23,59,59)


#--
_S1_POL_LST = ['VH','VV']
_S2_BND_LST = ['1', '2', '3', '4', '5', '6', '7', '8', '8a', '9', '10', '11', '12']
NB_S2 = len(_S2_BND_LST)

class ObsTable(object):
    """Class to hold time-series of observations together with satellite overpasses"""
    def __init__(self):
        self.geom    = None  #-- should become instance of SensorGeometry class
        self.data    = None  #-- should after I/O yield a 2D numpy array (nt,ndata)
                             #   ndata=2 for S1, ndata=13 for S2
        self.dataunc = None
        self.sat_flat = None #-- should become list of satellite identifier

class RetrievalSetup(object):

    def __init__(self, **kwargs):
        self.lon        = kwargs.get('lon',None)
        self.lat        = kwargs.get('lat',None)
        self.time_start = kwargs['time_start']
        self.time_end   = kwargs['time_end']
        #-- target schedule
        self.target_select   = kwargs.get('target_select', None)
        self.target_schedule = kwargs.get('target_schedule', None)
        self.date_lst_xtgt   = None

        #-- maps missions to sentinel overpass times/geometries
        self.obs_dct = OrderedDict()
        #-- temporally consecutively ordered observation points
        self.schedule_dct = OrderedDict( [('date_utc', []),
                                          ('sentinel', []),
                                          ('sim_typ',  []),
                                          ('sza',      []),
                                          ('saa',      []),
                                          ('vza',      []),
                                          ('vaa',      [])  ] )
        #-- prior control vector and uncertainties
        self.lai_coeff     = kwargs.get('laicoeff',0.1)
        self.lai_coeff_unc = self.lai_coeff 
        self.prstate       = None #-- becomes numpy array
        self.prstate_unc   = None #-- becomes numpy array
        self.generic_prior     = [1.5, 1.5, 0.3] #LAI[m2/m2], canht[m], SM[m3/m3]
        self.generic_prior_unc = [5.0, 2.0, 0.2] #LAI[m2/m2], canht[m], SM[m3/m3]
        #--
        self.prior_states_file = kwargs.get('prior_states_file', None)
        self.prior_inifile     = kwargs.get('prior_inifile',     None)
        self.use_generic_prior = kwargs.get('use_generic_prior', False)

        #-- observational *output* fill value as expected by retrieval system
        #   (ref. mo_sensimul.f90)
        self.obs_fill_value = -99999.

        #-- observational uncertainty settings
        #-- S1 uncertainty: value of 0.4 (as provided by Bjoern Rommen)
        self.s1_unc_db      = kwargs.get('s1_unc',      0.4)
        self.s1_vv_uncfloor = kwargs.get('s1_vv_uncfloor', None)
        self.s1_vh_uncfloor = kwargs.get('s1_vh_uncfloor', None)
        self.s1_floor_db    = kwargs.get('s1_floor',   -22.)
        self.s1_pol         = kwargs.get('s1_pol',     None)
        self.s2_relunc      = kwargs.get('s2_relunc',  0.05)
        self.s2_uncfloor    = kwargs.get('s2_uncfloor',0.01*self.s2_relunc)
        self.s2_bnds        = kwargs.get('s2_bnds',    None)


        #-- control vector uncertainty settings
        self.lai_coeff_relunc   = kwargs.get('lai_coeff_relunc', 1.0)  # 100%
        self.lai_relunc         = kwargs.get('lai_relunc',       0.5)  #  50%
        self.canht_relunc       = kwargs.get('canht_relunc',     0.5)  #  50%
        self.sm_relunc          = kwargs.get('sm_relunc',        0.5)  #  50%
        self.relunc_by_user     = kwargs.has_key('lai_relunc')
        self.lai_coeff_uncfloor = kwargs.get('lai_coeff_uncfloor', 1.e-3)
        self.lai_uncfloor       = kwargs.get('lai_uncfloor',       1.e-3)
        self.canht_uncfloor     = kwargs.get('canht_uncfloor',     1.e-3)
        self.sm_uncfloor        = kwargs.get('sm_uncfloor',        1.e-3)

        #-- statevector physical bounds
        self.bounds_dct = OrderedDict()
        dct = self.bounds_dct # just abbreviate...
        dct['lai_coeff_lobound'] = kwargs.get('lai_coeff_lobound', BOUNDS_DICT['lai_coeff_lobound'])
        dct['lai_coeff_hibound'] = kwargs.get('lai_coeff_hibound', BOUNDS_DICT['lai_coeff_hibound'])
        dct['lai_lobound']       = kwargs.get('lai_lobound',       BOUNDS_DICT['lai_lobound'])
        dct['lai_hibound']       = kwargs.get('lai_hibound',       BOUNDS_DICT['lai_hibound'])
        dct['canht_lobound']     = kwargs.get('canht_lobound',     BOUNDS_DICT['canht_lobound'])
        dct['canht_hibound']     = kwargs.get('canht_hibound',     BOUNDS_DICT['canht_hibound'])
        dct['sm_lobound']        = kwargs.get('sm_lobound',        BOUNDS_DICT['sm_lobound'])
        dct['sm_hibound']        = kwargs.get('sm_hibound',        BOUNDS_DICT['sm_hibound'])

        #-- dynamical model uncertainties
        self.dynmodunc_inifile = kwargs.get('dynmodunc_inifile',None)
        #   (default values, may be overwritten from external file)
        #
        self.dynmodunc_dct = OrderedDict()
        #-- sigma_LAI=1          dt<5days
        #-- sigma_LAI=2   5days<=dt<10days
        #-- sigma_LAI=5  10days<=dt
        self.dynmodunc_dct['lai']   = {'ndays_lst':[5,10],  'unc_lst':[1.,2.,5.]}
        #-- sigma_canht=1m           dt<10days
        #-- sigma_canht=2m   10days<=dt<30days
        #-- sigma_canht=5m   30days<=dt
        self.dynmodunc_dct['canht'] = {'ndays_lst':[10,30], 'unc_lst':[1,2,5]}
        #-- sigma_sm=0.25           dt<1days
        #-- sigma_sm=0.5     1days<=dt<5days
        #-- sigma_sm=1.0     5days<=dt
        self.dynmodunc_dct['sm']    = {'ndays_lst':[1,5],   'unc_lst':[0.25,0.5,1.]}

        #-- inversion settings
        self.use_prior  = kwargs.get('use_prior', True)
        self.use_model  = kwargs.get('use_model', True)
        self.gtol       = kwargs.get('gtol',     1.e-5)
        self.prior_pert = kwargs.get('prior_pert', 0.)


        #-- NetCDF output generation
        self.zlev     = kwargs.get('zlev',4)
        self.use_zlib = True if self.zlev>0 else False
    # ---__init__---

    def get_npts(self):
        return len(self.schedule_dct['date_utc'])

    def has_obs_s1(self):
        return self.obs_dct.has_key('S1')

    def has_obs_s2(self):
        return self.obs_dct.has_key('S2')

    def load_obs_csv(self, csv_file, date_fmt="%Y/%m/%d %H:%M", mission_lst=None, only_geom=False):
        """Function to read Sentinel observations together with overpass time and
        illumination-view geometries from a csv file which are to be mapped into
        the internal dictionaries of the RetrievalSetup class.

        :param csv_file: name of csv file providing observations
        :type cvs_file:  string
        """

        try:
            obs_data = np.loadtxt(csv_file, delimiter=',', dtype='str')
            msg = "observation data loaded from file ***{}***".format(csv_file)
            FileLogger.info(msg)
        except IOError as exc:
            msg = "could not load observations from csv file ***{}***".format(csv_file)
            msg += " ({})".format(exc)
            FileLogger.fatal(msg)
            raise RuntimeError(msg)

        nt,ncol = obs_data.shape
        date_lst = [ dt.datetime.strptime(obs_data[i,0], date_fmt) for i in xrange(nt) ]
        date_a   = np.array(date_lst)
        time_start_data = date_lst[0]
        time_end_data   = date_lst[-1]
        #-- logging
        msg = "detected ntimepts={} #columns={} in csv file".format(nt, ncol)
        FileLogger.info(msg)

        #-- potential adjustment to specified temporal domain
        if self.time_start!=None:
            time_start = self.time_start
        else:
            time_start = time_start_data
        if self.time_end!=None:
            time_end = self.time_end
        else:
            time_end = time_end_data

        #-- first 8 columns are always:date, vza, vaa, sza, saa, sat_flag, lat, lon

        if ncol==10:
            msg = "start reading S1 observations..."
            FileLogger.info(msg)
            #   date, vza, vaa, sza, saa, sat_flag, lat, lon, vh, vv
            vh_lst = []
            vv_lst = []
            self.obs_dct['S1'] = ObsTable()
            self.obs_dct['S1'].geom = satgeo.SensorGeometry()
            self.obs_dct['S1'].sat_id_lst = []
            #-- abreviate
            sat_geom = self.obs_dct['S1'].geom
            sat_geom.date_utc = []
            sat_geom.vza      = []
            sat_geom.vaa      = []
            sat_geom.sza      = []
            sat_geom.saa      = []
            for i,act_date in enumerate(date_lst):
                if act_date<time_start:
                    continue
                elif act_date>time_end:
                    break
                #-- actual satellite/mission
                act_mission = obs_data[i,7].upper()
                if mission_lst!=None and not act_mission in mission_lst:
                    msg = "observation at date {} is from mission={} and ignored here.".format(
                        act_date.strftime('%Y-%m-%dT%H:%M'), act_mission)
                    FileLogger.info(msg)
                    continue
                #-- read actual geometry
                sat_geom.date_utc.append(act_date)
                sat_geom.vza.append( float(obs_data[i,1]) )
                sat_geom.vaa.append( float(obs_data[i,2]) )
                sat_geom.sza.append( float(obs_data[i,3]) )
                sat_geom.saa.append( float(obs_data[i,4]) )
                #-- lon,lat (columns 5,6) not needed
                #-- satellite flag (column 7)
                self.obs_dct['S1'].sat_id_lst.append(act_mission)
                #-- VH,VV in 0-indexed columns 8,9
                vh_lst.append( float(obs_data[i,8]) )
                vv_lst.append( float(obs_data[i,9]) )

            #-- geometries/satellite flags are done here
            if only_geom:
                return

            #-- turn into arrays
            vh = np.array(vh_lst)
            vv = np.array(vv_lst)
            #-- logging
            msg = "observational backscatter values are assumed to be in linear units!"
            FileLogger.info(msg)
            msg = "VH backscatter values read: VH[linear] min/max={}/{}".format(
                vh.min(), vh.max())
            FileLogger.info(msg)
            msg = "VV backscatter values read: VV[linear] min/max={}/{}".format(
                vv.min(), vv.max())
            FileLogger.info(msg)
            #-- uncertainty computation
            #-- XX_db = XX_db(XX)  = 10*log10(XX)
            #-- XX    = XX(XX_db)  = 10**(XX_db/10)
            #
            # for the uncertainty in linear/raw unit we apply conservative estimation:
            # 2*sXX = [ XX(XX_db+sXX_db) - XX(XX_db-sXX_db) ] (XX=VH,VV)
            #       = [ XX(XX_db)*10**(sXX_db/10.) - XX(XX_db)*10**(-sXX_db/10.)]
            #       = XX(XX_db)*[10**(sXX_db/10.) - 10**(-sXX_db/10.)]
            #       = XX * [10**(sXX_db/10.) - 10**(-sXX_db/10.)]
            ds = 0.5* (10**(self.s1_unc_db/10.) - 10**(-1*self.s1_unc_db/10.))
            #-- S1 uncertainty floor *may* be user-supplied
            if self.s1_vv_uncfloor!=None:
                dsvv_floor = self.s1_vv_uncfloor
            else:
                dsvv_floor = 10**(self.s1_floor_db/10.)*ds
            if self.s1_vh_uncfloor!=None:
                dsvh_floor = self.s1_vh_uncfloor
            else:
                dsvh_floor = 10**(self.s1_floor_db/10.)*ds
            msg = "assuming S1 observational uncertainty of {} [dB] ".format(self.s1_unc_db)
            msg += "yields relative uncertainty of {} [linear unit].".format(ds)
            FileLogger.info(msg)
            msg = "assuming vv={} vh={} S1 observational uncertainty floor [linear unit].".format(
                dsvv_floor, dsvh_floor)
            FileLogger.info(msg)
            svh = np.maximum(vh*ds, dsvh_floor)
            svv = np.maximum(vv*ds, dsvv_floor)
            #-- apply floor value
            nlo_svh = np.count_nonzero(vh*ds<dsvh_floor)
            nlo_svv = np.count_nonzero(vv*ds<dsvv_floor)
            svh = np.maximum(svh, dsvh_floor)
            svv = np.maximum(svv, dsvv_floor)
            msg = "number of applied uncertainty floor values on VH={} VV={}".format(
                nlo_svh, nlo_svv)
            FileLogger.info(msg)
            msg = "determined VH uncertainty in linear units, min/max={}/{}".format(
                svh.min(), svh.max())
            FileLogger.info(msg)
            msg = "determined VV uncertainty in linear units, min/max={}/{}".format(
                svv.min(), svv.max())
            FileLogger.info(msg)
            #-- potential filtering of polarisations
            if not self.s1_pol is None:
                if not 'VH' in self.s1_pol:
                    vh  = self.obs_fill_value
                    svh = self.obs_fill_value
                if not 'VV' in self.s1_pol:
                    vv  = self.obs_fill_value
                    svv = self.obs_fill_value
            #-- 
            nt_use = len(sat_geom.date_utc)
            self.obs_dct['S1'].data = np.empty((nt_use,2), dtype=np.float64) #-- 'VH','VV'
            self.obs_dct['S1'].data[:,0] = vh
            self.obs_dct['S1'].data[:,1] = vv
            self.obs_dct['S1'].dataunc = np.empty((nt_use,2), dtype=np.float64)
            self.obs_dct['S1'].dataunc[:,0] = svh
            self.obs_dct['S1'].dataunc[:,1] = svv
            #-- logging
            msg = "...reading S1 observations DONE"
            FileLogger.info(msg)
        else:
            #-- logging
            msg = "start reading S2 observations..."
            FileLogger.info(msg)
            #   date, vza, vaa, sza, saa, sat_flag, lat, lon, BRF1,...,BRF13
            self.obs_dct['S2'] = ObsTable()
            self.obs_dct['S2'].geom = satgeo.SensorGeometry()
            self.obs_dct['S2'].sat_id_lst = []
            #-- abreviate
            sat_geom = self.obs_dct['S2'].geom
            sat_geom.date_utc = []
            sat_geom.vza      = []
            sat_geom.vaa      = []
            sat_geom.sza      = []
            sat_geom.saa      = []
            brf_lst = [ [] for i in xrange(NB_S2) ] #-- prepare lists for 13 BRF bands
            for i,act_date in enumerate(date_lst):
                if act_date<time_start:
                    continue
                elif act_date>time_end:
                    break
                #-- actual satellite/mission
                act_mission = obs_data[i,7].upper()
                if mission_lst!=None and not act_mission in mission_lst:
                    msg = "observation at date {} is from mission={} and ignored here.".format(
                        act_date.strftime('%Y-%m-%dT%H:%M'), act_mission)
                    FileLogger.info(msg)
                    continue
                #-- read actual geometry
                sat_geom.date_utc.append(act_date)
                sat_geom.vza.append( float(obs_data[i,1]) )
                sat_geom.vaa.append( float(obs_data[i,2]) )
                sat_geom.sza.append( float(obs_data[i,3]) )
                sat_geom.saa.append( float(obs_data[i,4]) )
                #-- lon/lat in columns 5, 6 not used here
                #-- satellite flag
                self.obs_dct['S2'].sat_id_lst.append(obs_data[i,7])
                #-- BRFs start at 0-indexed column 8 in data csv file
                for ib in xrange(NB_S2):
                    icol = ib+8
                    brf_lst[ib].append( float(obs_data[i, icol]) )

            #-- geometries/satellite flags are done here
            if only_geom:
                return
            #--
            nt_use = len(sat_geom.date_utc)
            brf_data = np.empty((nt_use,NB_S2), dtype=np.float64) #-- BRF1-13
            for ib in xrange(NB_S2):
                brf_data[:,ib] = np.array(brf_lst[ib])
            #-- check observational consistency
            nneg = np.count_nonzero( brf_data<0 )
            if nneg>0:
                msg = "detected negative BRF values: nneg={}.".format(nneg)
                msg += " These will be set to fill-value!"
                FileLogger.warn(msg)
                brf_data[ brf_data<0 ] = self.obs_fill_value
            nhi = np.count_nonzero( brf_data>1 )
            if nhi>0:
                msg = "detected high BRF outlier values>1: nout={}.".format(nhi)
                msg += " These will be set to fill-value!"
                FileLogger.warn(msg)
                brf_data[ brf_data>1 ] = self.obs_fill_value

            #-- data uncertainty
            msg = "BRF uncertainty is derived by applying {} relative uncertainty, ".format(
                self.s2_relunc)
            msg += "and an uncertainty floor value of {}".format(self.s2_uncfloor)
            FileLogger.info(msg)
            brf_dataunc = np.maximum(brf_data*self.s2_relunc, self.s2_uncfloor)
            brf_dataunc[ brf_dataunc<0 ] = self.obs_fill_value
            brf_dataunc[ brf_data==self.obs_fill_value ] = self.obs_fill_value
            #-- restriction to seleted bands
            if not self.s2_bnds is None:
                bnd_msk = np.ones((NB_S2,), dtype=np.bool)*True
                bnd_msk[self.s2_bnds] = False
                brf_data[:,bnd_msk]    = self.obs_fill_value
                brf_dataunc[:,bnd_msk] = self.obs_fill_value
            #-- set into structure
            self.obs_dct['S2'].data    = brf_data
            self.obs_dct['S2'].dataunc = brf_dataunc
            #-- logging
            msg = "...reading S2 observations DONE"
            FileLogger.info(msg)
    # ---load_obs_csv---


    def _guess_time_format(self, csv_file):
        """Function that trys to determine the time-format in a states csv file,
        where the time/date is expected to be in the first column.

        Note: this method is to be used internally only and should not be explicitly
        user-invoked.

        :param csv_file: name of csv file providing observations
        :type cvs_file:  string
        """
        import csv

        fmt_lst = ['%Y/%m/%d %H:%M', '%Y-%m-%d %H:%M:%S']

        fmt_found = None

        with open(csv_file,'r') as fp:
            reader = csv.DictReader(fp)
            for i,line in enumerate(reader):
                for k,v in line.iteritems():
                    if k.find('date')>=0: #-- this should be the date column
                        date_str = v
                        break
                if i>0:
                    break

        msg = "found first date in file ---{}---".format(v)
        FileLogger.info(msg)

        for fmt in fmt_lst:
            try:
                dt.datetime.strptime(date_str,fmt)
                fmt_found = fmt
                break
            except ValueError:
                pass

        msg = "detected time-format '{}'".format(fmt_found)
        FileLogger.info(msg)

        return fmt_found
    # ---_guess_time_format---


    def _read_target_schedule(self):
        """Function that determines the target schedule time-point specified by
        the user, either via command-line or a file.

        Note: this method is to be used internally only and should not be explicitly
        user-invoked.
        """
        target_schedule = self.target_schedule
        target_select   = self.target_select

        #-- list of scheduled times to be returned
        self.date_lst_xtgt = []

        if target_schedule!=None:
            if target_select!=None:
                msg = "target_schedule and target_select were both specified, " \
                      "'target_select' will be ignored here!"
                FileLogger.warn(msg)
            msg = "target schedule will be read from file ***{}***".format(target_schedule)
            FileLogger.info(msg)
            with open(target_schedule,'r') as fp:
                for line in fp:
                    if line[0]=='#':
                        continue
                    else:
                        # print "---{}---".format(line.rstrip())
                        date_utc = dt.datetime.strptime(line.rstrip(),'%Y%m%dT%H:%M:%S')
                        if self.time_start!=None and date_utc<self.time_start:
                            continue
                        elif self.time_end!=None and date_utc>self.time_end:
                            continue
                        else:
                            self.date_lst_xtgt.append(date_utc)
                #-- ensure time-increase ordering
                self.date_lst_xtgt = sorted(self.date_lst_xtgt)
                nxtgt = len(self.date_lst_xtgt)
                msg = "...reading target schedule DONE (nxtgt={})".format(nxtgt)
                FileLogger.info(msg)
        elif target_select!=None:
            msg = "target schedule will be determined from specification ---{}---".format(
                target_select)
            FileLogger.info(msg)
            ttgt_min = target_select[0]
            ttgt_max = target_select[1]
            ttgt_delta = target_select[2]
            ttgt_min = dt.datetime.strptime(ttgt_min,'%Y%m%dT%H:%M')
            ttgt_max = dt.datetime.strptime(ttgt_max,'%Y%m%dT%H:%M')
            if ttgt_delta[-1].lower()=='h':
                ttgt_delta = dt.timedelta(hours=float(ttgt_delta[0:-1]))
            elif ttgt_delta[-1].lower()=='d':
                ttgt_delta = dt.timedelta(days=float(ttgt_delta[0:-1]))
            date_utc = ttgt_min
            while date_utc<=ttgt_max:
                if self.time_start!=None and date_utc<self.time_start:
                    pass
                elif self.time_end!=None and date_utc>self.time_end:
                    pass
                else:
                    self.date_lst_xtgt.append(date_utc)
                #-- increment date
                date_utc += ttgt_delta
            msg = "read {} state components at extra target times.".format(len(self.date_lst_xtgt))
            FileLogger.info(msg)


    def setup_common_schedule(self):
        """

        """
        #-- merge dates from different sensors and extra target points
        date_lst = []
        #-- from observations/geometries
        for m,obs in self.obs_dct.iteritems():
            sensor_geom = obs.geom #-- Instance of SensorGeometry
            date_lst += sensor_geom.date_utc

        #-- user target schedule
        self._read_target_schedule()
        date_lst += self.date_lst_xtgt
        #-- make list unique, ensure ascending order
        date_lst = sorted( list(set(date_lst)) )
        ntpts = len(date_lst)
        msg = "determined ntpts={}".format(ntpts)
        FileLogger.info(msg)

        if ntpts==0:
            msg = "specified configuration does not yield any usable time-points."
            msg += "Cannot prepare for a retrieval!"
            FileLogger.fatal(msg)
            raise RuntimeError(msg)

        #-- store scheduled dates
        self.schedule_dct['date_utc'] = date_lst
        # print "-"*30
        # for x in date_lst:
        #     print x.strftime('%Y%m%dT%H:%M')
        # print "-"*30

        #-- determine
        #   - illumination-view geometry
        #   - simulation type (which Sentinel)
        #
        for i,d in enumerate(date_lst):
            #-- mission observation times
            for m,obs_pair in self.obs_dct.iteritems():
                geom = obs_pair.geom
                satid_lst = obs_pair.sat_id_lst
                try:
                    #-- index of actual date in observation dictionary for m=S1 or m=S2
                    im = geom.date_utc.index(d)
                    self.schedule_dct['sentinel'].append(m)
                    if m=='S1':
#                        self.schedule_dct['sim_typ'].append(1) #-- !!! set BIT1 !!!
                        if satid_lst[im]=='S1A':
                            self.schedule_dct['sim_typ'].append(1+2**2) #-- !!! set BIT1,3 !!!
                        else:  #S1B
                            self.schedule_dct['sim_typ'].append(1+2**3) #-- !!! set BIT1,4 !!!
                    elif m=='S2':
                        # self.schedule_dct['sim_typ'].append(2) #-- !!! set BIT2 !!!
                        if satid_lst[im]=='S2A':
                            self.schedule_dct['sim_typ'].append(2+2**4) #-- !!! set BIT2,5 !!!
                        else: #S2B
                            self.schedule_dct['sim_typ'].append(2+2**5) #-- !!! set BIT2,6 !!!
                    else:
                        msg = "unexpected m={} found in obs_dct".format(m)
                        FileLogger.fatal(msg)
                        raise RuntimeError(msg)
                    #-- add geometries
                    self.schedule_dct['sza'].append( geom.sza[im] )
                    self.schedule_dct['saa'].append( geom.saa[im] )
                    self.schedule_dct['vza'].append( geom.vza[im] )
                    self.schedule_dct['vaa'].append( geom.vaa[im] )
                except ValueError:
                    pass
            #-- extra target time
            try:
                i_x = self.date_lst_xtgt.index(d)
                self.schedule_dct['sim_typ'].append(0)
                self.schedule_dct['sza'].append( -999. )
                self.schedule_dct['saa'].append( -999. )
                self.schedule_dct['vza'].append( -999. )
                self.schedule_dct['vaa'].append( -999. )
            except ValueError:
                pass
        #-- informational logging
        sim_typ_a = np.array(self.schedule_dct['sim_typ'])
        n_s1a   = np.count_nonzero( (sim_typ_a&2**2)==2**2 )
        n_s1b   = np.count_nonzero( (sim_typ_a&2**3)==2**3 )
        n_s2a   = np.count_nonzero( (sim_typ_a&2**4)==2**4 )
        n_s2b   = np.count_nonzero( (sim_typ_a&2**5)==2**5 )
        n_other = np.count_nonzero( sim_typ_a==0 )
        msg = "common schedule yields time-points  S1A={} S1B={}".format(
            n_s1a, n_s1b)
        msg += " S2A={} S2B={} other={}".format(n_s2a,n_s2b,n_other)
        FileLogger.info(msg)


    def set_prior_priorunc_synthetic(self):
        """
        Function that sets either
        - the default state vector shipped with the prototype tool or
        - the generic prior for agricultural sites
        """

        lai_coeff_absunc = None
        statevec_absunc  = None

        #-- 
        if self.prior_inifile!=None:
            lai_coeff_absunc, statevec_absunc = self._setprior_from_inifile()
        elif self.use_generic_prior:
            self._setprior_generic_agriculture()
            statevec_absunc = self.generic_prior_unc
        else:
            #-- overall number of time-points in schedule
            npts = self.get_npts()

            #-- default prior file
            prior_file = os.path.join(ipt_dir_path, 'mni_stat_jules_2017.csv')

            #-- get signature simulator default state
            msg = "START reading state variables from file ***{}***...".format(prior_file)
            FileLogger.info(msg)
            state_inst = sv.get_state_csv(fname=prior_file, fmt='%Y-%m-%d %H:%M:%S' )
            msg = "...reading DONE"
            FileLogger.info(msg)

            #-- LAI,Canopy-Height,Soil-Moisture
            self.prstate     = np.empty((3,npts), dtype=np.float64)

            for i,date_utc in enumerate(self.schedule_dct['date_utc']):
                idx, timedelt = sv.find_nearest_date_idx(state_inst.date_utc, date_utc)
                # print "MVMV::nearest={} idx={} timedelt={}".format(
                #     state_inst.date_utc[idx], idx, timedelt)
                #-- LAI
                self.prstate[0,i]     = state_inst.lai[idx]
                #-- canopy-height
                self.prstate[1,i]     = state_inst.can_height[idx]
                #-- SM
                self.prstate[2,i]     = state_inst.soil_moisture[idx]

        #-- set uncertainty values
        self._set_priorunc(statevec_absunc=statevec_absunc, lai_coeff_absunc=lai_coeff_absunc)

    def set_prior_priorunc_general(self):
        """
        Function that sets the prior control vector by either
        - using an external states file
        - using the prior information provided by a prior ini file
        - applying generic prior for agricultural sites
        """

        #-- some configurations apply absolute uncertainties
        lai_coeff_absunc = None
        statevec_absunc  = None
        is_generic_prior = False

        #--
        if self.prior_states_file!=None:
            states_file = self.prior_states_file
            basename    = os.path.basename(states_file)
            if os.path.splitext(basename)[1]=='.nc':
                msg = "Prior state information will be read from ***{}***".format(states_file)
                FileLogger.info(msg)
                self._setprior_jules(states_file)
                msg = "...reading prior DONE"
                FileLogger.info(msg)
            elif os.path.splitext(basename)[1]=='.csv':
                msg = "Prior state information will be read from ***{}***".format(states_file)
                FileLogger.info(msg)
                self._setprior_csv(states_file)
                msg = "...reading prior DONE"
                FileLogger.info(msg)
            else:
                msg = "Unrecognised format of states file ***{}***. Cannot continue!".format(
                    states_file)
                FileLogger.fatal(msg)
                raise RuntimeError(msg)
                return
        elif self.prior_inifile!=None:
            lai_coeff_absunc, statevec_absunc = self._setprior_from_inifile()
        else:
            self._setprior_generic_agriculture()
            is_generic_prior = True
            statevec_absunc = self.generic_prior_unc

        #-- set uncertainty values
        self._set_priorunc( lai_coeff_absunc=lai_coeff_absunc,
                            statevec_absunc=statevec_absunc,
                            is_generic_prior=is_generic_prior  )

    def _setprior_from_inifile(self):
        from ConfigParser import SafeConfigParser
        """
        Function to set a temporally constant prior control vector with values and
        uncertainties read from a user specified ini file.
        For the LAI coefficient (the S1 parameter) only the uncertainty must
        be specified, since it's value is taken from the site.nml file.
        For lai, canopy height and soil moisture the value and the uncertainty
        must be specified.
        """

        #-- number of time-points
        npts = self.get_npts()

        #-- LAI,Canopy-Height,Soil-Moisture
        self.prstate     = np.empty((3,npts), dtype=np.float64)

        inifile = self.prior_inifile
        if not os.path.exists(inifile):
            msg = "user-specified prior ini file does not exist ***{}***".format(inifile)
            FileLogger.fatal(msg)
            raise RuntimeError(msg)
        else:
            parser = SafeConfigParser()
            parser.read(inifile)
            lai_coeff_absunc = float(parser.get('lai_coeff','uncertainty'))
            statevec_absunc = []
            for i,c in enumerate(['lai','canht','sm']):
                self.prstate[i,:] = float(parser.get(c,'value'))
                statevec_absunc.append(float(parser.get(c,'uncertainty')))
            return (lai_coeff_absunc,statevec_absunc)

    def _setprior_generic_agriculture(self):
        """
        Function to set a generic state vector for an agricultural site and
        the given temporal schedule. Following values are applied uniformily
        at all time points. These values are internally stored in
        lists 'self.generic_prior', 'self.generic_prior_unc', the
        values should be:
        LAI          : value=1.5, unc=5.  [m2/m2]
        canopy-height: value=1.5, unc=2.  [m]
        Soilmoisture : value=0.3, unc=0.2 [m3/m3]
        """

        #-- number of time-points
        npts = self.get_npts()

        #-- LAI,Canopy-Height,Soil-Moisture
        self.prstate     = np.empty((3,npts), dtype=np.float64)
        #-- LAI
        self.prstate[0,:]     = self.generic_prior[0]
        #-- canopy-height
        self.prstate[1,:]     = self.generic_prior[1]
        #-- soil moisture (volumetric)
        self.prstate[2,:]     = self.generic_prior[2]


    def _setprior_csv(self, csv_file):
        """
        Function that sets the state vector from csv file compliant with signature simulator.

        :param csv_ncfile: file name
        :type csv_ncfile:  string
        """

        #-- number of time-points
        npts = self.get_npts()

        #-- read state from CSV file
        fmt = self._guess_time_format(csv_file)
        state_inst = sv.get_state_csv(fname=csv_file, fmt=fmt)

        #-- LAI,Canopy-Height,Soil-Moisture
        self.prstate     = np.empty((3,npts), dtype=np.float64)

        for i,date_utc in enumerate(self.schedule_dct['date_utc']):
            idx, timedelt = sv.find_nearest_date_idx(state_inst.date_utc, date_utc)
            if timedelt.days>=1:
                msg = "for scheduled date ---{}--- ".format(date_utc.strftime('%Y-%m-%dT%H%M'))
                msg += "time nearest state differs by at least one day!"
                FileLogger.warn(msg)
            #-- LAI
            self.prstate[0,i]     = state_inst.lai[idx]
            #-- canopy-height
            self.prstate[1,i]     = state_inst.can_height[idx]
            #-- SM
            self.prstate[2,i]     = state_inst.soil_moisture[idx]


    def _setprior_jules(self, jules_ncfile):
        """
        Function that sets the state vector for the given site and
        temporal schedule by reading from JULES output in NetCDF format.

        :param jules_ncfile: file name of JULES generated NetCDF file for single spatial location.
        :type jules_ncfile:  string
        """

        #-- number of time-points
        npts = self.get_npts()

        #-- read state components from JULES
        #   IMPL-ATTENTION::pft_idx=5 should be fine for maize-crop and
        #                   thus suitable for Wallerfing, ***BUT*** must
        #                   be adjusted for other sites.
        state_inst = sv.get_jules_state(jules_ncfile, pft_idx=5)

        #-- LAI,Canopy-Height,Soil-Moisture
        self.prstate     = np.empty((3,npts), dtype=np.float64)

        for i,date_utc in enumerate(self.schedule_dct['date_utc']):
            idx, timedelt = sv.find_nearest_date_idx(state_inst.date_utc, date_utc)
            if timdelt.days>=1:
                msg = "for scheduled date ---{}--- ".format(date_utc.strftime('%Y-%m-%dT%H%M'))
                msg += "time nearest state differs by at least one day!"
                FileLogger.warn(msg)
            #-- LAI
            self.prstate[0,i]     = state_inst.lai[idx]
            #-- canopy-height
            self.prstate[1,i]     = state_inst.can_height[idx]
            #-- SM
            self.prstate[2,i]     = state_inst.soil_moisture[idx]


    def _set_priorunc(self, lai_coeff_absunc=None, statevec_absunc=None, is_generic_prior=False):
        """
        Function that sets the uncertainty of the (prior) control vector.
        """
        if not self.__dict__.has_key('prstate') or self.prstate is None:
            msg = "internal error, prior state does not yet exist!"
            FileLogger.fatal(msg)
            raise RuntimeError(msg)
        if lai_coeff_absunc!=None:
            self.lai_coeff_unc    = np.maximum( lai_coeff_absunc, self.lai_coeff_uncfloor )
            msg = "applied absolute uncertainty on lai coefficient, "
            msg += "lai_coeff_absunc={}".format(lai_coeff_absunc)
            FileLogger.info(msg)
        else:
            self.lai_coeff_unc    = np.maximum( self.lai_coeff*self.lai_coeff_relunc,
                                                self.lai_coeff_uncfloor               )
            msg = "applied relative uncertainty {} on lai coefficient.".format(
                self.lai_coeff_relunc)
            FileLogger.info(msg)

        #-- allocate space for uncertainty array
        self.prstate_unc = np.empty(self.prstate.shape, dtype=np.float64)

        #--   generic prior, user-supplied relative uncertainties
        if is_generic_prior and self.relunc_by_user:
            self.prstate_unc[0,:] = np.maximum( self.prstate[0,:]*self.lai_relunc,
                                                self.lai_uncfloor                  )
            self.prstate_unc[1,:] = np.maximum( self.prstate[1,:]*self.canht_relunc,
                                                self.canht_uncfloor                  )
            self.prstate_unc[2,:] = np.maximum( self.prstate[2,:]*self.sm_relunc,
                                                self.sm_uncfloor                  )
            msg = "applied relative uncertainty on state vector components, "
            msg += "lai_relunc={} canht_relunc={} sm_relunc={}".format(
                self.lai_relunc, self.canht_relunc, self.sm_relunc)
            FileLogger.info(msg)
        #--   either: prior information from ini-file or generic prior
        elif statevec_absunc!=None:
            lai_absunc,canht_absunc,sm_absunc = statevec_absunc
            self.prstate_unc[0,:] = np.maximum( lai_absunc,   self.lai_uncfloor )
            self.prstate_unc[1,:] = np.maximum( canht_absunc, self.canht_uncfloor )
            self.prstate_unc[2,:] = np.maximum( sm_absunc,    self.sm_uncfloor )
            msg = "applied absolute uncertainty on state vector components, "
            msg += "lai_absunc={} canht_absunc={} sm_absunc={}".format(
                lai_absunc, canht_absunc, sm_absunc)
            FileLogger.info(msg)
        #--   relative uncertainty (default or user-supplied)
        else:
            self.prstate_unc[0,:] = np.maximum( self.prstate[0,:]*self.lai_relunc,
                                                self.lai_uncfloor                  )
            self.prstate_unc[1,:] = np.maximum( self.prstate[1,:]*self.canht_relunc,
                                                self.canht_uncfloor                  )
            self.prstate_unc[2,:] = np.maximum( self.prstate[2,:]*self.sm_relunc,
                                                self.sm_uncfloor                  )
            msg = "applied relative uncertainty on state vector components, "
            msg += "lai_relunc={} canht_relunc={} sm_relunc={}".format(
                self.lai_relunc, self.canht_relunc, self.sm_relunc)
            FileLogger.info(msg)
        #-- logging
        msg = "uncertainty floor values were applied as follows: "
        msg += "lai_uncfloor={} canht_uncfloor={} sm_uncfloor={}".format(
            self.lai_uncfloor, self.canht_uncfloor, self.sm_uncfloor)
        FileLogger.info(msg)

    def set_dynmodel(self):
        from ConfigParser import SafeConfigParser

        dynmodunc_dct = self.dynmodunc_dct #- abbreviate only

        #-- read dynamical model uncertainties
        #   from external file (if given)
        if self.dynmodunc_inifile==None:
            msg = "Applying default values on dynamical model uncertainties."
            FileLogger.info(msg)
        else:
            msg = "START reading values on dynamical model uncertainties from ***{}***...".format(
                self.dynmodunc_inifile)
            FileLogger.info(msg)
            parser = SafeConfigParser()
            parser.read(self.dynmodunc_inifile)
            for state in ['lai','canht','sm']:
                ndays_lst = [ int(v) for v in parser.get(state,'ndays_lst').split(',')]
                unc_lst   = [ float(v) for v in parser.get(state,'unc_lst').split(',')]
                if len(ndays_lst)+1!=len(unc_lst):
                    msg = "inconistent sizes in ndays_lst/unc_lst for state '{}'".format(state)
                    msg += " #ndays_lst={} #unc_lst={}".format(len(ndays_lst),len(unc_lst))
                    raise RuntimeError(msg)
                #-- update values in internal dictionary
                dynmodunc_dct[state] = {'ndays_lst':ndays_lst, 'unc_lst':unc_lst}
            msg = "...reading DONE"
            FileLogger.info(msg)

        #-- debug
        msg = "dynmodunc_dct => {}".format(dynmodunc_dct)
        FileLogger.debug(msg)

        #-- get time-points
        timepts = self.schedule_dct['date_utc']
        npts    = len(timepts)

        #-- determine model compontents (a,b,munc) with
        #   a    = (a_lai,   a_canht,   a_sm)
        #   b    = (b_lai,   b_canht,   b_sm)
        #   munc = (munc_lai,munc_canht,munc_sm)
        #
        self.a_components    = np.empty((3,npts), dtype=np.float64)
        self.b_components    = np.empty((3,npts), dtype=np.float64)
        self.munc_components = np.empty((3,npts), dtype=np.float64)

        #-- computation of dyn model uncertainty,
        #   depending on state and temporal difference (in days)
        def modelunc_(state, dt_):
            ndays_lst = dynmodunc_dct[state]['ndays_lst']
            unc_lst   = dynmodunc_dct[state]['unc_lst']
            for i,nd in enumerate(ndays_lst):
                if dt_.days<nd:
                    return unc_lst[i]
            #-- temporal difference equal/above maximum in ndays
            return unc_lst[2]

        self.a_components[:,:]    = 1.
        self.b_components[:,:]    = 0.
        self.munc_components[:,:] = 1.
        for i in range(1,npts):
            delta_t = timepts[i] - timepts[i-1]
            self.munc_components[0,i] = modelunc_('lai',  delta_t)
            self.munc_components[1,i] = modelunc_('canht',delta_t)
            self.munc_components[2,i] = modelunc_('sm',   delta_t)

# %%%RetrievalSetup%%%


#=============================
#
#          mkdirp_smart
#
def mkdirp_smart(newdir):
    """
    works the way a good 'mkdir --parents' should :)
    - already exists, silently complete
    - regular file in the way, raise an exception
    - parent directory(ies) does not exist, make them as well

    :param newdir: string specifying directory name
    :type newdir:  string
    """
    import os

    if os.path.isdir(newdir):
        pass
    elif os.path.exists(newdir):
        msg = "a non-directory path with the same name as the desired " \
              "dir, =***{}***, already exists.".format(newdir)
        raise RuntimeError(msg)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            mkdirp_smart(head)
        if tail:
            try:
                os.mkdir(newdir)
            except OSError as exc:
                msg = "newdir={} could not be created on system (exc={})".format(
                    newdir, exc)
                raise RuntimeError(msg)
# ---mkdirp_smart---


#=============================
#
#         datestr_parse
#
def datestr_parse(date_str):
    """
    Turns string into a datetime.datetime object

    :param date_str: string describing a date in suitable format
    :type date_str:  string
    :return:         datetime.datetime object generated from the input string
    :rtype:          datetime.datetime instance
    """

    fmt_lst = ['%Y%m%dT%H:%M:%S', '%Y%m%dT%H:%M', '%Y%m%d']
    dtime = None
    def _get_dtime(fmt):
        try:
            return dt.datetime.strptime(date_str, fmt)
        except ValueError:
            return dtime

    for fmt in fmt_lst:
        dtime = _get_dtime(fmt)
        if dtime!=None:
            break

    if dtime is None:
        msg = "date_str ***{}*** could not properly be handled!".format(date_str)
        raise RuntimeError(msg)

    return dtime
# ---datestr_parse---


#=============================
#
#         datelst_get_month_aligned_bounds
#
def datelst_get_month_aligned_bounds(dates_):
    """Function that returns monthly aligned boundary dates for the given list/series of dates

    :param dates_: list/array of datetime objects
    :return: tuple of datetime objects (lower, upper bound)
    :rtype: tuple
    """
    dfirst = dates_[0]
    dlast  = dates_[-1]

    bound_lo = dt.datetime(dfirst.year, dfirst.month, 1)
    bound_hi = (dt.datetime(dates_[-1].year, dates_[-1].month, 1)+dt.timedelta(days=32))
    bound_hi.replace(day=1)
    bound_hi = bound_hi.replace(day=1) - dt.timedelta(seconds=1)

    return (bound_lo, bound_hi)
# ---datelst_get_month_aligned_bounds---


#=============================
#
#         set_site_params
#
def set_site_params(site_nml):
    """
    Function that sets values of site_param_dic from values specified in a fortran namelist file
    :param site_nml: path to fortran nml file
    :return: None
    """
    nml_dic = f90nml.read(site_nml)
    for key in nml_dic['site_params'].keys():
        if type(nml_dic['site_params'][key]) == float:
            self.site_param_dic[key] = nml_dic['site_params'][key]
        elif type(nml_dic['site_params'][key]) == list:
            self.site_param_dic[key] = nml_dic['site_params'][key][0]
# ---set_site_params---


#=============================
#
#         ncwrt_retrieval_config
#
def ncwrt_retrieval_config( retr_setup, outname=None ):
    """Function that writes the (common S1+S2) state-vector of the Sentinel Simulator
    including required ancillary information (i.e. illumination-view geometries) to a NetCDF
    file suitable as input to the retrieval system of the S^3 Project.

    :param retr_setup: Retrieval Setup container
    :type retr_setup: Instance of RetrievalSetup class
    :param outname: name of generated output file
    :type outname: string
    """

    #-- set name of file to be generated
    act_outname = outname if outname!=None else 'retrconfig.nc'
    msg = "Start writing configuration file ***{}***...".format(act_outname)
    FileLogger.info(msg)


    #-- compression settings
    zlev     = retr_setup.zlev
    use_zlib = retr_setup.use_zlib

    #--
    schedule_dct  = retr_setup.schedule_dct
    statevector   = retr_setup.prstate
    #-- turn list into array
    sim_typ       = np.array(schedule_dct['sim_typ'], dtype=np.int32)
    timepts       = schedule_dct['date_utc']
    nstvar,npts   = statevector.shape
    #-- overpass geometries SZA,SAA,VZA,VAA
    ivgeom = np.empty((npts,4), dtype=np.float64)
    ivgeom[:,0] = schedule_dct['sza']
    ivgeom[:,1] = schedule_dct['saa']
    ivgeom[:,2] = schedule_dct['vza']
    ivgeom[:,3] = schedule_dct['vaa']

    #-- temporal settings, create time-values (time-since)
    time_start, time_end = datelst_get_month_aligned_bounds(timepts)
    time_coverage_start  = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    time_coverage_end    = time_end.strftime('%Y-%m-%dT%H:%M:%S')
    ref_time = dt.datetime(timepts[0].year,1,1) #January1st in year of first point in time
    time_unit = 'seconds since {}'.format(ref_time.strftime('%Y-%m-%dT%H:%M:%S'))
    time_values = nc4.date2num(timepts, time_unit)

    #-- ensure directory exists
    mkdirp_smart(os.path.dirname(act_outname))

    #-- open file pointer
    ncfp = nc4.Dataset(act_outname, 'w')
    #-- add dimensions
    d1 = ncfp.createDimension('npoints',npts)
    d2 = ncfp.createDimension('ngeo',4)

    #-- time-value
    ncvar = ncfp.createVariable( 'time', statevector.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name','time')
    ncvar.setncattr('long_name','time')
    ncvar.setncattr('units', time_unit)
    ncvar[:] = time_values[:]

    #-- simulation type
    ncvar = ncfp.createVariable( 'sim_typ', sim_typ.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar[:] = sim_typ[:]
    ncvar.setncattr('long_name','simulation_type')
    ncvar.setncattr('comment', 'integer value which is to be bit-interpreted')
    ncvar.setncattr('nobits_set',  'time-point with other state')
    ncvar.setncattr('bit0_is_set', 'time-point for S1 simulation')
    ncvar.setncattr('bit1_is_set', 'time-point for S2 simulation')
    ncvar.setncattr('bit2_is_set', 'time-point for S1A simulation')
    ncvar.setncattr('bit3_is_set', 'time-point for S1B simulation')
    ncvar.setncattr('bit4_is_set', 'time-point for S2A simulation')
    ncvar.setncattr('bit5_is_set', 'time-point for S2B simulation')
    
    #-- illumination-view geometry
    ncvar = ncfp.createVariable( 'ivgeom', ivgeom.dtype, ('npoints','ngeo'),
                                 zlib=use_zlib, complevel=zlev              )
    ncvar.setncattr('sza','igeo: 0')
    ncvar.setncattr('saa','igeo: 1')
    ncvar.setncattr('vza','igeo: 2')
    ncvar.setncattr('vaa','igeo: 3')
    ncvar[:,:] = ivgeom[:,:]
    
    #-- global attributes
    ncfp.setncattr('creator_name',"The Inversion Lab, Hamburg, Germany")
    ncfp.setncattr('creator_email', "Michael.Vossbeck(at)Inversion-Lab.com")
    ncfp.setncattr('netcdf_libversion',"{}".format(nc4.__netcdf4libversion__))
    ncfp.setncattr('date_created',"{}".format(dt.datetime.utcnow().isoformat()))
    ncfp.setncattr('time_coverage_start',time_coverage_start)
    ncfp.setncattr('time_coverage_end',time_coverage_end)

    #-- close file pointer
    ncfp.close()

    # logging
    msg = "...writing ***{}*** DONE".format(act_outname)
    FileLogger.info(msg)
# ---ncwrt_retrieval_config---


#=============================
#
#         ncwrt_retrieval_prior
#
def ncwrt_retrieval_prior(retr_setup, outname=None):
    """Function that writes the prior state vector (LAI,Canopy-Height,Soilmoisture)
    including uncertainties to a NetCDF file suitable as input
    to the retrieval system of the S^3 Project.

    :param retr_setup: Retrieval Setup container
    :type retr_setup: Instance of RetrievalSetup class
    :param outname: name of generated output file
    :type outname: string
    """

    #-- set name of file to be generated
    act_outname = outname if outname!=None else 'retrprior.nc'
    msg = "Start writing configuration file ***{}***...".format(act_outname)
    FileLogger.info(msg)

    #-- compression settings
    zlev     = retr_setup.zlev
    use_zlib = retr_setup.use_zlib

    #--
    schedule_dct    = retr_setup.schedule_dct
    statevector     = retr_setup.prstate
    statevector_unc = retr_setup.prstate_unc
    timepts         = schedule_dct['date_utc']
    nstvar,npts     = statevector.shape

    #-- temporal settings, create time-values (time-since)
    time_start, time_end = datelst_get_month_aligned_bounds(timepts)
    time_coverage_start  = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    time_coverage_end    = time_end.strftime('%Y-%m-%dT%H:%M:%S')
    ref_time = dt.datetime(timepts[0].year,1,1) #January1st in year of first point in time
    time_unit = 'seconds since {}'.format(ref_time.strftime('%Y-%m-%dT%H:%M:%S'))
    time_values = nc4.date2num(timepts, time_unit)

    #-- ensure directory exists
    mkdirp_smart(os.path.dirname(act_outname))

    #-- open file pointer
    ncfp = nc4.Dataset(act_outname, 'w')
    #-- add dimensions
    d1 = ncfp.createDimension('npoints',npts)
    d2 = ncfp.createDimension('nparam_s1',1)

    #-- time-value
    ncvar = ncfp.createVariable( 'time', statevector.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name','time')
    ncvar.setncattr('long_name','time')
    ncvar.setncattr('units', time_unit)
    ncvar[:] = time_values[:]

    #-- S1 lai_coeff
    ncvar = ncfp.createVariable( 'param_s1', np.float64, ('nparam_s1',),
                                 zlib=use_zlib, complevel=zlev               )
    ncvar.setncattr('long_name', 'model parameter involved in simulation of S1 signal')
    ncvar[:] = retr_setup.lai_coeff
    #-- S1 parameter uncertainty
    ncvar = ncfp.createVariable( 'param_unc_s1', np.float64, ('nparam_s1',),
                                 zlib=use_zlib, complevel=zlev               )
    ncvar.setncattr('long_name', '1-sigma uncertainty of model parameter involved in simulation of S1 signal')
    #-- check consistency
    if retr_setup.lai_coeff_unc==None:
        msg = "uncertainty of LAI coefficient has not been set. This is an internal error!"
        FileLogger.fatal(msg)
        raise RuntimeError(msg)
    else:
        ncvar[:] = retr_setup.lai_coeff_unc

    #-- LAI
    ncvar = ncfp.createVariable( 'lai', statevector.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name', 'LAI')
    ncvar.setncattr('long_name', 'Leaf area index')
    ncvar.setncattr('units','m2/m2')
    ncvar[:] = statevector[0,:]
    #-- LAI uncertainty
    ncvar = ncfp.createVariable( 'lai_unc', statevector.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name', 'LAI_unc')
    ncvar.setncattr('long_name', '1-sigma uncertainty of Leaf area index')
    ncvar.setncattr('units','m2/m2')
    # comment = "uncertainty was derived as {} [%] relative uncertainty ".format(
    #     100.*retr_setup.lai_relunc)
    # comment += "and an uncertainty floor value of {} was applied.".format(
    #     retr_setup.bounds_dct['lai_lobound'])
    # ncvar.setncattr('comment', comment)
    ncvar[:] = statevector_unc[0,:]

    #-- canopy height
    ncvar = ncfp.createVariable( 'canht', statevector.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name', 'canht')
    ncvar.setncattr('long_name', 'Canopy height')
    ncvar.setncattr('units','m')
    ncvar[:] = statevector[1,:]
    #-- canopy height uncertainty
    ncvar = ncfp.createVariable( 'canht_unc', statevector.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name', 'canht_unc')
    ncvar.setncattr('long_name', '1-sigma uncertainty of Canopy height')
    ncvar.setncattr('units','m')
    # comment = "uncertainty was derived as {} [%] relative uncertainty ".format(
    #     100.*retr_setup.canht_relunc)
    # comment += "and an uncertainty floor value of {} was applied.".format(
    #     retr_setup.bounds_dct['canht_lobound'])
    # ncvar.setncattr('comment', comment)
    ncvar[:] = statevector_unc[1,:]

    #-- soil moisture
    ncvar = ncfp.createVariable( 'sm', statevector.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name', 'SM')
    ncvar.setncattr('long_name', 'Soil moisture')
    ncvar.setncattr('units','m3/m3')
    ncvar[:] = statevector[2,:]
    #-- soil moisture uncertainty
    ncvar = ncfp.createVariable( 'sm_unc', statevector.dtype, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name', 'SM_unc')
    ncvar.setncattr('long_name', '1-sigma uncertainty of Soil moisture')
    ncvar.setncattr('units','m3/m3')
    # comment = "uncertainty was derived as {} [%] relative uncertainty ".format(
    #     100.*retr_setup.sm_relunc)
    # comment += "and an uncertainty floor value of {} was applied.".format(
    #     retr_setup.bounds_dct['sm_lobound'])
    # ncvar.setncattr('comment', comment)
    ncvar[:] = statevector_unc[2,:]

    #-- global attributes
    ncfp.setncattr('creator_name',"The Inversion Lab, Hamburg, Germany")
    ncfp.setncattr('creator_email', "Michael.Vossbeck(at)Inversion-Lab.com")
    ncfp.setncattr('netcdf_libversion',"{}".format(nc4.__netcdf4libversion__))
    ncfp.setncattr('date_created',"{}".format(dt.datetime.utcnow().isoformat()))
    ncfp.setncattr('time_coverage_start',time_coverage_start)
    ncfp.setncattr('time_coverage_end',time_coverage_end)

    #-- close file pointer
    ncfp.close()

    # logging
    msg = "...writing ***{}*** DONE".format(act_outname)
    FileLogger.info(msg)
# ---ncwrt_retrieval_prior---


#=============================
#
#         ncwrt_retrieval_model
#
def ncwrt_retrieval_model(retr_setup, outname=None):
    """Function that writes the (common S1+S2) state-vector of the Sentinel Simulator
    including required ancillary information (i.e. illumination-view geometries) to a NetCDF
    file suitable as input to the retrieval system of the S^3 Project.

    :param retr_setup: Retrieval Setup container
    :type retr_setup: Instance of RetrievalSetup class
    :param outname: name of generated output file
    :type outname: string
    """

    #-- set name of file to be generated
    act_outname = outname if outname!=None else 'retrmodel.nc'
    msg = "Start writing configuration file ***{}***...".format(act_outname)
    FileLogger.info(msg)

    #-- compression settings
    zlev     = retr_setup.zlev
    use_zlib = retr_setup.use_zlib

    #-- retrieval settings
    schedule_dct    = retr_setup.schedule_dct
    timepts         = schedule_dct['date_utc']
    npts            = len(timepts)
    a_components    = retr_setup.a_components
    b_components    = retr_setup.b_components
    munc_components = retr_setup.munc_components

    #-- temporal settings, create time-values (time-since)
    time_start, time_end = datelst_get_month_aligned_bounds(timepts)
    time_coverage_start  = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    time_coverage_end    = time_end.strftime('%Y-%m-%dT%H:%M:%S')
    ref_time = dt.datetime(timepts[0].year,1,1) #January1st in year of first point in time
    time_unit = 'seconds since {}'.format(ref_time.strftime('%Y-%m-%dT%H:%M:%S'))
    time_values = nc4.date2num(timepts, time_unit)

    #-- ensure directory exists
    mkdirp_smart(os.path.dirname(act_outname))

    #-- open file pointer
    ncfp = nc4.Dataset(act_outname, 'w')
    #-- add dimensions
    d1 = ncfp.createDimension('npoints',npts)

    #-- time-value
    ncvar = ncfp.createVariable( 'time', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name','time')
    ncvar.setncattr('long_name','time')
    ncvar.setncattr('units', time_unit)
    ncvar[:] = time_values[:]

    unit_one = np.float64(1)
    #--     a-coefficient part
    # lai
    ncvar = ncfp.createVariable( 'a_lai', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'a coefficient associated to LAI in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = a_components[0,:]
    # canht
    ncvar = ncfp.createVariable( 'a_canht', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'a coefficient associated to canopy height in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = a_components[1,:]
    # soil moisture
    ncvar = ncfp.createVariable( 'a_sm', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'a coefficient associated to soil moisture in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = a_components[2,:]

    #--     offset-coefficient part
    # lai
    ncvar = ncfp.createVariable( 'offset_lai', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'offset associated to LAI in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = b_components[0,:]
    # canht
    ncvar = ncfp.createVariable( 'offset_canht', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'offset associated to canopy height in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = b_components[1,:]
    # soil moisture
    ncvar = ncfp.createVariable( 'offset_sm', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'offset associated to soil moisture in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = b_components[2,:]

    #--     model-uncertainty part
    # lai
    ncvar = ncfp.createVariable( 'munc_lai', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'model uncertainty associated to LAI in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = munc_components[0,:]
    # canht
    ncvar = ncfp.createVariable( 'munc_canht', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'model uncertainty associated to canopy height in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = munc_components[1,:]
    # soil moisture
    ncvar = ncfp.createVariable( 'munc_sm', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev      )
    ncvar.setncattr('long_name', 'model uncertainty associated to soil moisture in dynamic state model')
    ncvar.setncattr('units',unit_one)
    ncvar[:] = munc_components[2,:]


    #-- global attributes
    ncfp.setncattr('creator_name',"The Inversion Lab, Hamburg, Germany")
    ncfp.setncattr('creator_email', "Michael.Vossbeck(at)Inversion-Lab.com")
    ncfp.setncattr('netcdf_libversion',"{}".format(nc4.__netcdf4libversion__))
    ncfp.setncattr('date_created',"{}".format(dt.datetime.utcnow().isoformat()))
    ncfp.setncattr('time_coverage_start',time_coverage_start)
    ncfp.setncattr('time_coverage_end',time_coverage_end)

    #-- close file pointer
    ncfp.close()

    # logging
    msg = "...writing ***{}*** DONE".format(act_outname)
    FileLogger.info(msg)
# ---ncwrt_retrieval_model---


#=============================
#
#         ncwrt_retrieval_obs_s1
#
def ncwrt_retrieval_obs_s1(retr_setup, outname=None):
    """Function that writes observations from Sentinel1 ('VH','VV' polarisations)
    including uncertainties to a NetCDF file suitable as input
    to the retrieval system of the S^3 Project.

    :param retr_setup: Retrieval Setup container
    :type retr_setup: Instance of RetrievalSetup class
    :param outname: name of generated output file
    :type outname: string
    """

    #-- set name of file to be generated
    act_outname = outname if outname!=None else 'obs_s1.nc'
    msg = "Start writing configuration file ***{}***...".format(act_outname)
    FileLogger.info(msg)

    #-- compression settings
    zlev     = retr_setup.zlev
    use_zlib = retr_setup.use_zlib

    #-- retrieval settings
    s1_table   = retr_setup.obs_dct['S1']
    timepts    = s1_table.geom.date_utc
    npts       = len(timepts)
    s1_satid   = np.array(s1_table.sat_id_lst, dtype=str)
    s1_data    = s1_table.data
    s1_dataunc = s1_table.dataunc
    nt,npol = s1_data.shape

    #-- temporal settings, create time-values (time-since)
    time_start, time_end = datelst_get_month_aligned_bounds(timepts)
    time_coverage_start  = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    time_coverage_end    = time_end.strftime('%Y-%m-%dT%H:%M:%S')
    ref_time = dt.datetime(timepts[0].year,1,1) #January1st in year of first point in time
    time_unit = 'seconds since {}'.format(ref_time.strftime('%Y-%m-%dT%H:%M:%S'))
    time_values = nc4.date2num(timepts, time_unit)

    #-- ensure directory exists
    mkdirp_smart(os.path.dirname(act_outname))

    #-- open file pointer
    ncfp = nc4.Dataset(act_outname, 'w')
    #-- add dimensions
    d1 = ncfp.createDimension('npoints',npts)
    d2 = ncfp.createDimension('npol',npol)

    #-- time-value
    ncvar = ncfp.createVariable( 'time', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev           )
    ncvar.setncattr('standard_name','time')
    ncvar.setncattr('long_name','time')
    ncvar.setncattr('units', time_unit)
    ncvar[:] = time_values[:]

    unit_one = np.float64(1)
    # backscatter
    ncvar = ncfp.createVariable( 'backscatter', np.float64, ('npoints','npol'),
                                 zlib=use_zlib, complevel=zlev            )
    ncvar.setncattr('long_name', 'backscatter in VH and VV polarisation')
    comment = "VH is associated to npol=0, VV to npol=1."
    comment += " linear units are used (not [dB])."
    ncvar.setncattr('comment', comment)
    ncvar.setncattr('units',unit_one)
    ncvar.setncattr('missing_value', retr_setup.obs_fill_value)
    ncvar[:,:] = s1_data[:,:]

    # backscatter uncertainty
    ncvar = ncfp.createVariable( 'backscatter_unc', np.float64, ('npoints','npol'),
                                 zlib=use_zlib, complevel=zlev            )
    ncvar.setncattr('long_name', 'backscatter uncertainty in VH and VV polarisation')
    comment = "uniform uncertainty of {} [dB] was applied on the observed backscatter".format(
        retr_setup.s1_unc_db)
    ncvar.setncattr('comment', comment)
    ncvar.setncattr('units',unit_one)
    ncvar.setncattr('missing_value', retr_setup.obs_fill_value)
    ncvar[:,:] = s1_dataunc[:,:]

    # satellite identifier
    ncvar = ncfp.createVariable( 'satellite_id', str, ('npoints',),
                                 zlib=use_zlib, complevel=zlev            )
    ncvar.setncattr('long_name', 'satellite identifer')
    ncvar[:] = s1_satid[:]

    #-- global attributes
    ncfp.setncattr('creator_name',"The Inversion Lab, Hamburg, Germany")
    ncfp.setncattr('creator_email', "Michael.Vossbeck(at)Inversion-Lab.com")
    ncfp.setncattr('netcdf_libversion',"{}".format(nc4.__netcdf4libversion__))
    ncfp.setncattr('date_created',"{}".format(dt.datetime.utcnow().isoformat()))
    ncfp.setncattr('time_coverage_start',time_coverage_start)
    ncfp.setncattr('time_coverage_end',time_coverage_end)

    #-- close file pointer
    ncfp.close()

    # logging
    msg = "...writing ***{}*** DONE".format(act_outname)
    FileLogger.info(msg)

# ---ncwrt_retrieval_obs_s1---


#=============================
#
#         ncwrt_retrieval_obs_s2
#
def ncwrt_retrieval_obs_s2(retr_setup, outname=None):
    """Function that writes observations from Sentinel2 (BRF values)
    including uncertainties to a NetCDF file suitable as input
    to the retrieval system of the S^3 Project.

    :param retr_setup: Retrieval Setup container
    :type retr_setup1: Instance of RetrievalSetup class
    :param outname: name of generated output file
    :type outname: string
    """

    #-- set name of file to be generated
    act_outname = outname if outname!=None else 'obs_s2.nc'
    msg = "Start writing configuration file ***{}***...".format(act_outname)
    FileLogger.info(msg)

    #-- compression settings
    zlev     = retr_setup.zlev
    use_zlib = retr_setup.use_zlib

    #-- retrieval settings
    s2_table   = retr_setup.obs_dct['S2']
    timepts    = s2_table.geom.date_utc
    npts       = len(timepts)
    s2_satid   = np.array(s2_table.sat_id_lst, dtype=str)
    s2_data    = s2_table.data
    s2_dataunc = s2_table.dataunc
    nt,nbands  = s2_data.shape

    #-- temporal settings, create time-values (time-since)
    time_start, time_end = datelst_get_month_aligned_bounds(timepts)
    time_coverage_start  = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    time_coverage_end    = time_end.strftime('%Y-%m-%dT%H:%M:%S')
    ref_time = dt.datetime(timepts[0].year,1,1) #January1st in year of first point in time
    time_unit = 'seconds since {}'.format(ref_time.strftime('%Y-%m-%dT%H:%M:%S'))
    time_values = nc4.date2num(timepts, time_unit)

    #-- ensure directory exists
    mkdirp_smart(os.path.dirname(act_outname))

    #-- open file pointer
    ncfp = nc4.Dataset(act_outname, 'w')
    #-- add dimensions
    d1 = ncfp.createDimension('npoints',npts)
    d2 = ncfp.createDimension('nbands',nbands)

    #-- time-value
    ncvar = ncfp.createVariable( 'time', np.float64, ('npoints',),
                                 zlib=use_zlib, complevel=zlev     )
    ncvar.setncattr('standard_name','time')
    ncvar.setncattr('long_name','time')
    ncvar.setncattr('units', time_unit)
    ncvar[:] = time_values[:]

    #-- unit (in correct type)
    unit_one = np.array([1]).astype(s2_data.dtype)[0]

    # BRF
    ncvar = ncfp.createVariable( 'brf', np.float64, ('npoints','nbands'),
                                 zlib=use_zlib, complevel=zlev            )
    ncvar.setncattr('long_name', 'BRF top-of-canopy reflectances')
    ncvar.setncattr('units',unit_one)
    ncvar.setncattr('missing_value', retr_setup.obs_fill_value)
    ncvar[:,:] = s2_data[:,:]

    # BRF uncertainty
    ncvar = ncfp.createVariable( 'brf_unc', np.float64, ('npoints','nbands'),
                                 zlib=use_zlib, complevel=zlev                )
    ncvar.setncattr('long_name', 'BRF top-of-canopy reflectances')
    ncvar.setncattr('units',unit_one)
    ncvar.setncattr('missing_value', retr_setup.obs_fill_value)
    comment = "BRF uncertainties are derived as {:.2f}[%] relative uncertainty ".format(
        100.*retr_setup.s2_relunc)
    comment += "and an uncertainty floor value of {:.4f} is applied.".format(retr_setup.s2_uncfloor)
    ncvar.setncattr('comment', comment)
    ncvar[:,:] = s2_dataunc[:,:]

    # satellite identifier
    ncvar = ncfp.createVariable( 'satellite_id', str, ('npoints',),
                                 zlib=use_zlib, complevel=zlev            )
    ncvar.setncattr('long_name', 'satellite identifer')
    ncvar[:] = s2_satid[:]

    #-- global attributes
    ncfp.setncattr('creator_name',"The Inversion Lab, Hamburg, Germany")
    ncfp.setncattr('creator_email', "Michael.Vossbeck(at)Inversion-Lab.com")
    ncfp.setncattr('netcdf_libversion',"{}".format(nc4.__netcdf4libversion__))
    ncfp.setncattr('date_created',"{}".format(dt.datetime.utcnow().isoformat()))
    ncfp.setncattr('time_coverage_start',time_coverage_start)
    ncfp.setncattr('time_coverage_end',time_coverage_end)

    #-- close file pointer
    ncfp.close()

    # logging
    msg = "...writing ***{}*** DONE".format(act_outname)
    FileLogger.info(msg)

# ---ncwrt_retrieval_obs_s2---


#=============================
#
#         wrt_nml_retrieval_control
#
def wrt_nml_retrieval_control( retr_setup, outname=None ):
    """Function that writes the namelist file which controls the minimisation
    of the retrieval system.

    :param retr_setup: Retrieval Setup container
    :type retr_setup: Instance of RetrievalSetup class
    :param outname: name of generated output file
    :type outname: string
    """
    
    #-- set name of file to be generated
    act_outname = outname if outname!=None else 'retrctl.nml'
    msg = "Start writing configuration file ***{}***...".format(act_outname)
    FileLogger.info(msg)

    #-- ensure directory exists
    mkdirp_smart(os.path.dirname(act_outname))

    prior_ftn = ".true." if retr_setup.use_prior else ".false."
    prior_cmt = "! if .true. then J_p included in Eq. 1.5"
    state_ftn = ".true." if retr_setup.use_model else ".false."
    state_cmt = "! if .true. then J_m included in Eq. 1.5"
    grad_tol  = retr_setup.gtol
    grad_cmt  = "! stopping criterion for minimisation: rel. reduction in gradient norm"
    pr_pert   = retr_setup.prior_pert
    pr_cmt    = "! initial guess will be the prior times (1+prior_pert)"
    try:
        fp = open(act_outname,'w')
    except IOError:
        msg = "file ***{}*** could not be opened for writing".format(act_outname)
        FileLogger.fatal(msg)
        return
    
    fp.write("&RETRCTL" + '\n')
    fp.write("  retr_use_prior_term = {} {}".format(prior_ftn, prior_cmt) + '\n')
    fp.write("  retr_use_model_term = {} {}".format(state_ftn, state_cmt) + '\n')
    fp.write("  gradient_tol = {:e} {}".format(grad_tol, grad_cmt) + '\n')
    fp.write("  prior_pert = {:e} {}".format(pr_pert, pr_cmt) + '\n')
    fp.write("/" + '\n')
    
    #-- close handle
    fp.close()

    # logging
    msg = "...writing ***{}*** DONE".format(act_outname)
    FileLogger.info(msg)
# ---wrt_nml_retrieval_control---


#=============================
#
#         wrt_nml_controlvector_bounds
#
def wrt_nml_controlvector_bounds( retr_setup, outname=None):
    """Function that writes namelist file specifying the control vector bounds
    to be used within the minimisation.

    :param retr_setup: Retrieval Setup container
    :type retr_setup: Instance of RetrievalSetup class
    :param outname: name of generated output file
    :type outname: string
    """

    #-- set name of file to be generated
    act_outname = outname if outname!=None else 'retrctlvecbounds.nml'
    msg = "Start writing configuration file ***{}***...".format(act_outname)
    FileLogger.info(msg)

    #-- ensure directory exists
    mkdirp_smart(os.path.dirname(act_outname))

    try:
        fp = open(act_outname,'w')
    except IOError:
        msg = "file ***{}*** could not be opened for writing".format(act_outname)
        FileLogger.fatal(msg)
        return

    fp.write("&ctlvector_bounds" + '\n')
    for k,v in retr_setup.bounds_dct.iteritems():
        fp.write(" {} = {} \n".format(k,v))
    fp.write('/' + '\n')

    #-- close handle
    fp.close()

    # logging
    msg = "...writing ***{}*** DONE".format(act_outname)
    FileLogger.info(msg)
# ---wrt_nml_controlvector_bounds---


#=============================
#
#          pre_synthetic
#
def pre_synthetic(options):

    #-- dictionary used to pass collected information to
    #   RetrievalSetup initialiser
    setup_dct = {}

    #-- setup site file
    site_nml = options.site_nml
    if site_nml==None:
        if os.path.exists('site.nml'):
            msg = "'site.nml' from current working directory will be used."
            FileLogger.info(msg)
        else:
            mni_nml = os.path.join(ss_dir_path,'site.nml')
            shutil.copyfile(mni_nml, 'site.nml')
            msg = "'site.nml' was created as copy of Wallerfing default namelist file ***{}***".format(
            mni_nml)
            FileLogger.info(msg)
    elif not os.path.exists(site_nml):
        msg = "specified namelist file ***{}*** does not exist. Will not continue!".format(site_nml)
        FileLogger.fatal(msg)
        raise RuntimeError(msg)
    else:
        #-- site specification file must exist exactly as 'site.nml' in current working
        #   directory
        if site_nml=='site.nml':
            msg = "'site.nml' will be used. There is no need to specify this filename explicitly."
            FileLogger.info(msg)
        else:
            shutil.copyfile(site_nml, 'site.nml')
            msg = "Copied default namelist file ***{}*** to current working directory.".format(
                site_nml)
            FileLogger.info(msg)

    #-- load site file
    site_nml = 'site.nml' #-- this will be the file used by the prototype tool
    msg = "START reading namelist ***{}***...".format(site_nml)
    FileLogger.info(msg)
    nml_dct  = f90nml.read('site.nml')
    #-- get dictionary for the site
    site_dct = nml_dct['site_params']
    msg = "...reading DONE"
    FileLogger.info(msg)

    #-- copy attributes
    setup_dct['lon']       = site_dct['lon']
    setup_dct['lat']       = site_dct['lat']
    setup_dct['lai_coeff'] = site_dct['lai_coeff']

    #-- copy attributes
    setup_dct['lon']       = site_dct['lon']
    setup_dct['lat']       = site_dct['lat']
    setup_dct['lai_coeff'] = site_dct['lai_coeff']

    #-- retrieval system control
    setup_dct['use_prior'] = options.use_prior
    setup_dct['use_model'] = options.use_model
    setup_dct['prior_pert'] = 0.25

    #-- prior control vector settings
    setup_dct['use_generic_prior'] = options.use_generic_prior

    #-- prior inifile specifying a temporally constant prior control vector
    setup_dct['prior_inifile']     = options.prior_inifile

    #-- temporal range
    if options.time_start==None:
        setup_dct['time_start'] = TIME_START_SYN
    else:
        setup_dct['time_start'] = datestr_parse(options.time_start)
    if options.time_end==None:
        setup_dct['time_end'] = TIME_END_SYN
    else:
        time_end = datestr_parse(options.time_end)
        #-- take last second of day
        time_end = dt.datetime(time_end.year,time_end.month,time_end.day,23,59,59)
        setup_dct['time_end'] = time_end

    time_start = setup_dct['time_start']
    time_end   = setup_dct['time_end']
    num_days = (time_end-time_start).days
    # logging
    msg = "operating for temporal range {} - {} (num_days={})".format(
        time_start.isoformat(), time_end.isoformat(), num_days)
    FileLogger.info(msg)

    #-- extra target points
    setup_dct['target_select']   = options.target_select
    setup_dct['target_schedule'] = options.target_schedule

    #-- instantiate RetrievalSetup
    retrieval_setup = RetrievalSetup(**setup_dct)

    #-- load satellite overpasses
    #-- set missions
    mission_lst = options.mission_lst

    # logging
    msg = "selected mission_lst: {}".format(mission_lst)
    FileLogger.info(msg)
    #-- Sentinel-1
    s1_mission_lst = filter(lambda m: m.find('S1')==0, mission_lst)
    if len(s1_mission_lst)>0: #'S1A' in mission_lst or 'S1B' in mission_lst:
        geomfile = os.path.join(ipt_dir_path,'mni_s1_508_med_2017.csv')
        retrieval_setup.load_obs_csv(geomfile, mission_lst=s1_mission_lst, only_geom=True)
    #-- Sentinel-2
    s2_mission_lst = filter(lambda m: m.find('S2')==0, mission_lst)
    if len(s2_mission_lst)>0: #'S2A' in mission_lst or 'S2B' in mission_lst:
        geomfile = os.path.join(ipt_dir_path,'mni_s2_508_med_2017.csv')
        retrieval_setup.load_obs_csv(geomfile, mission_lst=s2_mission_lst, only_geom=True)

    #-- setup schedule (merging the missions and extra target points)
    retrieval_setup.setup_common_schedule()

    #-- setup prior and prior uncertainty control vector
    retrieval_setup.set_prior_priorunc_synthetic()

    #-- setup dynamical-model
    retrieval_setup.set_dynmodel()

    #-- write configuration file
    ncwrt_retrieval_config(retrieval_setup)

    #-- write prior file
    ncwrt_retrieval_prior(retrieval_setup)

    #-- write dynamical-model file
    ncwrt_retrieval_model(retrieval_setup)

    #-- write minimisation control file
    wrt_nml_retrieval_control(retrieval_setup)

    #-- write control vector bounds file
    wrt_nml_controlvector_bounds(retrieval_setup)

# ---pre_synthetic---


#=============================
#
#          pre_general
#
def pre_general(options):

    #-- logging
    msg = "start setup of general retrieval setup..."
    FileLogger.info(msg)

    #-- dictionary used to pass collected information to
    #   RetrievalSetup initialiser
    setup_dct = {}

    #-- setup site file
    site_nml = options.site_nml
    if site_nml==None:
        if os.path.exists('site.nml'):
            msg = "'site.nml' from current working directory will be used."
            FileLogger.info(msg)
        else:
            mni_nml = os.path.join(ss_dir_path,'site.nml')
            shutil.copyfile(mni_nml, 'site.nml')
            msg = "'site.nml' was created as copy of Wallerfing default namelist file ***{}***".format(
            mni_nml)
            FileLogger.info(msg)
    elif not os.path.exists(site_nml):
        msg = "specified namelist file ***{}*** does not exist. Will not continue!".format(site_nml)
        FileLogger.fatal(msg)
        raise RuntimeError(msg)
    else:
        #-- site specification file must exist exactly as 'site.nml' in current working
        #   directory
        if site_nml=='site.nml':
            msg = "'site.nml' will be used. There is no need to specify this filename explicitly."
            FileLogger.info(msg)
        else:
            shutil.copyfile(site_nml, 'site.nml')
            msg = "Copied default namelist file ***{}*** to current working directory.".format(
                site_nml)
            FileLogger.info(msg)

    #-- load site file
    if options.outdir!=None:
        site_nml = os.path.join(options.outdir,'site.nml')
        shutil.copyfile('site.nml', site_nml)
    else:
        site_nml = 'site.nml'
    msg = "START reading namelist ***{}***...".format(site_nml)
    FileLogger.info(msg)
    nml_dct  = f90nml.read(site_nml)
    #-- get dictionary for the site
    site_dct = nml_dct['site_params']
    msg = "...reading DONE"
    FileLogger.info(msg)
    
    #-- DEBUG
    if options.verbose:
        print "-"*30
        for k,v in site_dct.iteritems():
            print "{} => {}".format(k,v)
        print "-"*30

    #-- copy attributes
    setup_dct['lon']       = site_dct['lon']
    setup_dct['lat']       = site_dct['lat']
    setup_dct['lai_coeff'] = site_dct['lai_coeff']

    #-- extra target points
    setup_dct['target_select']   = options.target_select
    setup_dct['target_schedule'] = options.target_schedule

    #-- retrieval system control
    setup_dct['use_prior'] = options.use_prior
    setup_dct['use_model'] = options.use_model
    setup_dct['gtol']      = options.gtol
    setup_dct['prior_pert'] = 0.

    #-- uncertainty settings on observations
    if options.s1_unc!=None:
        setup_dct['s1_unc'] = options.s1_unc
    setup_dct['s1_vv_uncfloor'] = options.s1_vv_uncfloor
    setup_dct['s1_vh_uncfloor'] = options.s1_vh_uncfloor
    if options.s1_pol!=None:
        setup_dct['s1_pol'] = options.s1_pol
        msg = "Assimilation of S1 backscatter will be restricted to polarisations +++{}+++".format(
            options.s1_pol)
        FileLogger.info(msg)
    if options.s2_relunc!=None:
        setup_dct['s2_relunc'] = options.s2_relunc
    if options.s2_uncfloor!=None:
        setup_dct['s2_uncfloor'] = options.s2_uncfloor
    if options.s2_bnds!=None:
        s2_bnd_idxs = []
        for b in options.s2_bnds:
            s2_bnd_idxs.append(_S2_BND_LST.index(b))
        setup_dct['s2_bnds'] = np.array(s2_bnd_idxs)
        msg = "Assimilation of S2 data will be restricted to bands +++{}+++".format(options.s2_bnds)
        FileLogger.info(msg)
    #-- uncertainty settings on control-vector
    if options.ctlvec_relunc!=None:
        setup_dct['lai_coeff_relunc'] = options.ctlvec_relunc[0]
        setup_dct['lai_relunc']       = options.ctlvec_relunc[1]
        setup_dct['canht_relunc']     = options.ctlvec_relunc[2]
        setup_dct['sm_relunc']        = options.ctlvec_relunc[3]
    if options.ctlvec_uncfloor!=None:
        setup_dct['lai_coeff_uncfloor'] = options.ctlvec_uncfloor[0]
        setup_dct['lai_uncfloor']       = options.ctlvec_uncfloor[1]
        setup_dct['canht_uncfloor']     = options.ctlvec_uncfloor[2]
        setup_dct['sm_uncfloor']        = options.ctlvec_uncfloor[3]

    #-- dynamical model uncertainty inifile
    setup_dct['dynmodunc_inifile'] = options.dynmodunc_inifile

    #-- data file specifying prior states
    setup_dct['prior_states_file'] = options.states_file
    #-- prior inifile specifying a temporally constant prior control vector
    setup_dct['prior_inifile']     = options.prior_inifile

    #-- temporal range
    if options.time_start==None:
        setup_dct['time_start'] = options.time_start
    else:
        setup_dct['time_start'] = datestr_parse(options.time_start)
    if options.time_end==None:
        setup_dct['time_end'] = options.time_end
    else:
        time_end = datestr_parse(options.time_end)
        #-- take last second of day
        time_end = dt.datetime(time_end.year,time_end.month,time_end.day,23,59,59)
        setup_dct['time_end'] = time_end

    #-- instantiate RetrievalSetup
    retrieval_setup = RetrievalSetup(**setup_dct)

    #-- read observations
    obs_s1 = options.obs_s1
    obs_s2 = options.obs_s2
    no_obs_given = obs_s1==None and obs_s2==None
    if no_obs_given:
        msg = "No observational input files specified by user. "
        msg += "This means, that a retrieval may be tried with already existing observation files "
        msg += "***obs_s1.nc*** and/or ***obs_s2.nc*** or without any observations."
        FileLogger.warn(msg)
    if obs_s1!=None:
        retrieval_setup.load_obs_csv(obs_s1)
    if obs_s2!=None:
        retrieval_setup.load_obs_csv(obs_s2)

    #-- setup schedule
    #   NOTE: this merges the different missions and extra target points
    #         into a consecutive time-series
    retrieval_setup.setup_common_schedule()

    #-- setup prior control vector and its uncertainty
    #   NOTE: this must be done *AFTER* the schedule is completed!
    retrieval_setup.set_prior_priorunc_general()

    #-- setup dynamical-model
    retrieval_setup.set_dynmodel()

    #-- logging
    msg = "...general retrieval setup DONE"
    FileLogger.info(msg)

    #-- logging
    msg = "START writing input files for retrieval system..."
    FileLogger.info(msg)

    #-- write configuration file
    outname = 'retrconfig.nc'
    if options.outdir!=None:
        outname = os.path.join(options.outdir, outname)
    ncwrt_retrieval_config(retrieval_setup, outname=outname)

    #-- write prior file
    outname = 'retrprior.nc'
    if options.outdir!=None:
        outname = os.path.join(options.outdir, outname)
    ncwrt_retrieval_prior(retrieval_setup, outname=outname)

    #-- write S1 observations (if given)
    if retrieval_setup.has_obs_s1():
        outname = 'obs_s1.nc'
        if options.outdir!=None:
            outname = os.path.join(options.outdir, outname)
        ncwrt_retrieval_obs_s1(retrieval_setup, outname=outname )

    #-- write S2 observations (if given)
    if retrieval_setup.has_obs_s2():
        outname = 'obs_s2.nc'
        if options.outdir!=None:
            outname = os.path.join(options.outdir, outname)
        ncwrt_retrieval_obs_s2(retrieval_setup, outname=outname)

    #-- write dynamical-model file
    outname = 'retrmodel.nc'
    if options.outdir!=None:
        outname = os.path.join(options.outdir, outname)
    ncwrt_retrieval_model(retrieval_setup, outname=outname)

    #-- write minimisation control file
    outname = 'retrctl.nml'
    if options.outdir!=None:
        outname = os.path.join(options.outdir, outname)
    wrt_nml_retrieval_control(retrieval_setup, outname=outname)

    #-- write control vector bounds file
    outname = 'retrctlvecbounds.nml'
    if options.outdir!=None:
        outname = os.path.join(options.outdir, outname)
    wrt_nml_controlvector_bounds(retrieval_setup, outname=outname)

    #-- logging
    msg = "...writing input files for retrieval system DONE"
    FileLogger.info(msg)
# ---pre_general---


#=============================
#
#          create_argument_parser
#
def create_argument_parser(progname=None):
    """
    :param progname: name of main program/module
    :type progname:  string/None
    :return:         an ArgumentParser object
    :rtype:          argparse.ArgumentParser instance
    """
    import argparse

    #%%%%%%%%%%%%%%%%%%%
    #
    #               AliasedSubParsersAction
    #
    # Purpose:
    # - provide aliases for 'sub commands'
    #
    # NOTE:
    # - this feature is *not* available with python 2.7.6
    #   (it will be included (at least) since python 3.3.x)
    #
    #
    class AliasedSubParsersAction(argparse._SubParsersAction):
        class _AliasedPseudoAction(argparse.Action):
            def __init__(self, name, aliases, help):
                dest = name
                if aliases:
                    dest += ' (%s)' % ','.join(aliases)
                sup = super(AliasedSubParsersAction._AliasedPseudoAction, self)
                sup.__init__(option_strings=[], dest=dest, help=help) 

        def add_parser(self, name, **kwargs):
            if 'aliases' in kwargs:
                aliases = kwargs['aliases']
                del kwargs['aliases']
            else:
                aliases = []

            parser = super(AliasedSubParsersAction, self).add_parser(name, **kwargs)

            # Make the aliases work.
            for alias in aliases:
                self._name_parser_map[alias] = parser
            # Make the help text reflect them, first removing old help entry.
            if 'help' in kwargs:
                help = kwargs.pop('help')
                self._choices_actions.pop()
                pseudo_action = self._AliasedPseudoAction(name, aliases, help)
                self._choices_actions.append(pseudo_action)

            return parser
    # ---AliasSubParsersAction

    def _add_output_options(aparser):
        aparser.add_argument( '--outdir',
                              help="""whether generated file(s) should be located in extra directory. (ATTENTION: This option should only be used in order to inspect the list of files prepared by the preprocessing tool, running the retrieval tool expects these files in the current working directory! """ )
    # ---_add_output_options---

    def _add_common_options(aparser):
        aparser.add_argument( '--time_start',
                              help="""retrieval system will discard time-points before the date given here (fmt: YYYYMMDDTHH:MM). However further temporal restriction may result from invoking options --target_select or --target_schedule.""" )
        aparser.add_argument( '--time_end',
                              help="""retrieval system will discard time-points after the date given here (fmt: YYYYMMDDTHH:MM). Howeveer, further temporal restriction may result from invoking options --target_select or --target_schedule.""" )
        aparser.add_argument( '--target_select',
                              metavar=('target_tmin','target_tmax','target_tdelta'),
                              nargs=3,
                              help="selection of series of equally distributed extra target " \
                              "time points. This must be a comma separated triple of "\
                              "the form 'tmin,tmax,tdelta' where tmim,tmax are " \
                              "time-points with format 'YYYYMMDDTHH:MM' " \
                              "and tdelta must be specified as Xd or Xh with integer X. " \
                              "(Note that option --target_schedule has priority over setting " \
                              "made here!)" )
        aparser.add_argument( '--target_schedule',
                              help="specify name of a target schedule file from which extra " \
                                   "time points will be read. " \
                                   "This file must contain a single time-point per line " \
                                   "specified as 'YYYYMMDDTHH:MM:SS'.""" )
        

    #-- create main parser instance
    parser = argparse.ArgumentParser( prog=progname, usage=globals()['__doc__'] )

    #-------------------
    #     c o m m o n   o p t i o n s
    #
    parser.add_argument( '--verbose','-v',
                         action='store_true',
                         help="""dump some more debugging messages to terminal.""" )

    #----------------------------
    #     s u b c o m m a n d s
    #
    subparsers = parser.add_subparsers( title='Available Subcommands',
                                        metavar='command_list:',
                                        description='',
                                        dest='subcmds',
                                        help=''
                                        )
    #-------------------
    #      p r e _ s y n t h e t i c
    #
    xparser = subparsers.add_parser( 'pre_synthetic',
                                     help="""generation of all relevant control and input files for the running the retrieval system in synthetic mode by applying the default settings from both the signature simulator and the retrieval system envirionment.""" )
    #-- add common options
    _add_common_options(xparser)
    xparser.add_argument( '--site_nml',
                          help="""specification namelist file for selected site. If none is given a default namelist file suitable for the Wallerfing site will be used.""" )
    xparser.add_argument( '--use_generic_prior',
                          action='store_true',
                          help="""whether the generic prior suitable for agriculture sites should be used (By default a default prior from JULES shipped with the signature simulator is applied)""" )
    xparser.add_argument( '--prior_inifile',
                          help="""ini file to specify temporally constant prior control vector components and their uncertainties. (NOTE: this option has priority over --use_generic_prior)""" )
    xparser.add_argument( '--mission_lst',
                          nargs='+',
                          default=[],
                          choices=['S1A','S1B','S2A','S2B'],
                          help="""list of satellite missions for which synthetic observations should be generated (default: %(default)s)""" )
    xparser.add_argument( '--no_use_prior',
                          dest='use_prior',
                          action='store_false',
                          default=True,
                          help="""whether to discard the prior information for the retrieval""" )
    xparser.add_argument( '--no_use_model',
                          dest='use_model',
                          action='store_false',
                          help="""whether to discard the dynamical model on the state vector for the retrieval""" )
#    _add_output_options(xparser)

    #-------------------
    #      p r e _ g e n e r a l
    #
    xparser = subparsers.add_parser( 'pre_general',
                                     help="""generation of all relevant control and input files for the running the retrieval system for a general case according to user specifications.""" )
    #-- add common options
    _add_common_options(xparser)
    xparser.add_argument( '--site_nml',
                          help="""specification namelist file for selected site. By default an existing file 'site.nml' will be used or a default namelist file suitable for the Wallerfing site will be copied to 'site.nml'.""" )
    xparser.add_argument( '--states_file',
                          help="""File that specifies a time-series of state-variables for the selected site (LAI, Canopy-height, Soilmoisture). This file may be given as a csv file compliant with the signature simulator or a JULES generated NetCDF file for a single spatial location. In case of a missing states file a generic prior suitable for agricultural sites will be applied.""" )
    xparser.add_argument( '--obs_s1',
                          help="""signature simulator compliant data file providing time-series of Sentinel1 observations together with satellite overpass times and illumination-view geometries.""" )
    xparser.add_argument( '--obs_s2',
                          help="""signature simulator compliant data file providing time-series of Sentinel2A observations together with satellite overpass times and illumination-view geometries.""" )
    xparser.add_argument( '--no_use_prior',
                          dest='use_prior',
                          action='store_false',
                          default=True,
                          help="""whether to discard the prior information for the retrieval""" )
    xparser.add_argument( '--no_use_model',
                          dest='use_model',
                          action='store_false',
                          help="""whether to discard the dynamical model on the state vector for the retrieval""" )
    xparser.add_argument( '--gtol',
                          type=float,
                          default=1.e-5,
                          help="""gradient norm used as stopping-criterion within the minimisation routine.""" )
    xparser.add_argument( '--dynmodunc_inifile',
                          help="""ini file from where uncertainty specification for the 3 states in the dynamical model are read.""" )
    xparser.add_argument( '--prior_inifile',
                          help="""ini file to specify temporally constant prior control vector components and their uncertainties.""" )
    xparser.add_argument( '--s1_unc',
                          type=float,
                          help="""user-specified absolute uncertainty of S1 observations [dB]. (Note: this turns into relative uncertainty in linear-units which are internally used within the retrieval system.)""" )
    xparser.add_argument( '--s1_vv_uncfloor',
                          type=float,
                          help="""user-specified floor value for uncertainty of S1 VV-polarised observations (Note: value here is interpreted in linear units!)""" )
    xparser.add_argument( '--s1_vh_uncfloor',
                          type=float,
                          help="""user-specified floor value for uncertainty of S1 VH-polarised observations (Note: value here is interpreted in linear units!)""" )
    xparser.add_argument( '--s1_pol',
                          nargs='+',
                          choices=_S1_POL_LST,
                          help="""user-specified selection to restrict the assimilation to specific polarisations of the S1 observations (By default observations from both polarisations 'VH' and 'VV' will be used).""" )
    xparser.add_argument( '--s2_relunc',
                          type=float,
                          help="""user-specified relative uncertainty of S2 observations""" )
    xparser.add_argument( '--s2_uncfloor',
                          type=float,
                          help="""user-specified floor value for uncertainty of S2 observations""" )
    xparser.add_argument( '--s2_bnds',
                          nargs='+',
                          choices=_S2_BND_LST,
                          help="""user-specified selection to restrict the assimilation to using only a subset of the Sentinel2 observations (By default observations from all 13 bands will be used.)""" )
    xparser.add_argument( '--ctlvec_relunc',
                          type=float,
                          nargs=4,
                          metavar=('laicoeff_relunc','lai_relunc','canht_relunc','sm_relunc'),
                          help="""user-specified relative uncertainty of controlvector components (uniform in time)""" )
    xparser.add_argument( '--ctlvec_uncfloor',
                          type=float,
                          nargs=4,
                          metavar=('laicoeff_uncfloor','lai_uncfloor','canht_uncfloor','sm_uncfloor'),
                          help="""user-specified floor value for uncertainty of controlvector components (uniform in time)""" )
    _add_output_options(xparser)


    return parser
# ---create_argument_parser---


#=============================
#
#          main
#
def main(options):
    """
    :param options: command line arguments packed into a Namespace object.
    :type options:  namespace instance resulting from parsing the command line by ArgumentParser
    """

    #-------------------
    #     s y n t h e t i c   s e t u p
    #
    if options.subcmds=='pre_synthetic':
        msg = "user selection +++{}+++".format(options.subcmds)
        FileLogger.info(msg)

        pre_synthetic(options)

    #-------------------
    #     g e n e r a l   r e t r i e v a l   s y s t e m
    #
    if options.subcmds=='pre_general':
        msg = "user selection +++{}+++".format(options.subcmds)
        FileLogger.info(msg)

        pre_general(options)

# ---main---


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#                    M A I N
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':

    progname = os.path.basename(__file__) #determine filename

    #-----------------------------
    #     P R O G R A M   S T A R T
    #
    fmt = "%Y-%m-%dT%H:%M:%S.%f"
    ttstart = dt.datetime.now()
    FileLogger.info("{}::PROGRAM START::{}".format(progname, ttstart.strftime(fmt)))
    FileLogger.info( "command-line: {}".format(' '.join(sys.argv)) )

    #          p a r s e   c o m m a n d   l i n e
    #
    parser = create_argument_parser(progname=progname)
    options = parser.parse_args()

    #debug:
    if options.verbose:
        OUT_HDLR.setLevel( logging.DEBUG )
        FileLogger.setLevel( logging.DEBUG )


    main(options)
