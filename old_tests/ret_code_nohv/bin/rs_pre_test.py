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

available options are listed by invoking './rs_pre.py' on the command line


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
dir_path = os.path.dirname(os.path.realpath(__file__))
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
MISSION_DICT = OrderedDict( [('S1a', 'Sentinel-1a'), ('S1b', 'Sentinel-1b'),
                             ('S2a', 'Sentinel-2a')                        ] )


BOUNDS_DICT  = OrderedDict( [('lai_coeff_lobound', 1.e-5),
                             ('lai_coeff_hibound', 1.0),
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



class ObsTable(object):
    """Class to hold time-series of observations together with satellite overpasses"""
    def __init__(self):
        self.geom    = None  #-- should become instance of SensorGeometry class
        self.data    = None  #-- should after I/O yield a 2D numpy array (nt,ndata)
                             #   ndata=2 for S1, ndata=13 for S2
        self.dataunc = None

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
        #-- statevector and uncertainties
        self.lai_coeff   = kwargs.get('laicoeff',0.1)
        self.prstate     = None
        self.prstate_unc = None

        #-- observational *output* fill value as expected by retrieval system
        #   (ref. mo_sensimul.f90)
        self.obs_fill_value = -99999.

        #-- observational uncertainties
        #-- S1 uncertainty: value of 0.4 as provided by Bjoern Rommen
        self.s1_unc_db   = kwargs.get('s1_unc', 0.6)
        self.s2_relunc   = kwargs.get('s2_relunc', 0.05)
        self.s2_uncfloor = kwargs.get('s2_uncfloor',0.007)

        #-- control vector default relative uncertainties
        self.lai_coeff_relunc = kwargs.get('laicoeff_relunc', 0.5) # 100%
        self.lai_relunc       = kwargs.get('lai_relunc', 0.5)     #  50%
        self.canht_relunc     = kwargs.get('canht_relunc', 0.5)   #  50%
        self.sm_relunc        = kwargs.get('sm_relunc', 0.5)      #  50%

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
        self.dynmodunc_dct['lai']   = {'ndays_lst':[5,10],  'unc_lst':[0.2,0.4,1.]}
        #-- sigma_canht=1m           dt<10days
        #-- sigma_canht=2m   10days<=dt<30days
        #-- sigma_canht=5m   30days<=dt
        self.dynmodunc_dct['canht'] = {'ndays_lst':[10,30], 'unc_lst':[.25,.5,1.]}
        #-- sigma_sm=0.25           dt<1days
        #-- sigma_sm=0.5     1days<=dt<5days
        #-- sigma_sm=1.0     5days<=dt
        self.dynmodunc_dct['sm']    = {'ndays_lst':[1,5],   'unc_lst':[0.05,0.25,0.4]}

        #-- inversion settings
        self.use_prior = kwargs.get('use_prior', True)
        self.use_state = kwargs.get('use_state', True)
        self.gtol      = kwargs.get('gtol',     1.e0)
        self.prior_pert = kwargs.get('prior_pert', 0.)


        #-- NetCDF output generation
        self.zlev     = kwargs.get('zlev',4)
        self.use_zlib = True if self.zlev>0 else False
    # ---__init__---

    def get_npts(self):
        return len(self.schedule_dct['date_utc'])


    def has_obs_s1(self):
        return self.obs_dct.has_key('S1a') or self.obs_dct.has_key('S1b')

    def has_obs_s2(self):
        return self.obs_dct.has_key('S2a') # or self.obs_dct.has_key('S2b')

    def load_geom_synthetic(self, mission_lst):
        for m in mission_lst:
            if m not in ['S1b','S2a']:
                continue
            else:
                if m=='S1b':
                    geomfile = os.path.join(ipt_dir_path,'mni_geom_s1_2017.csv')
                else:
                    geomfile = os.path.join(ipt_dir_path,'mni_geom_s2_2017.csv')
                msg = "Start retrieving geometries for mission"
                msg += " ---{}--- from file ***{}***...".format(
                    MISSION_DICT[m], geomfile)
                FileLogger.info(msg)
                sensor_geom = satgeo.get_geom_csv(geomfile)
                msg = "...DONE (ngeom={})".format(len(sensor_geom.date_utc))
                FileLogger.info(msg)
                #-- potential adjustment to specified temporal domain
                if self.time_start!=None:
                    time_start = self.time_start
                else:
                    time_start = sensor_geom.date_utc[0]
                if self.time_end!=None:
                    time_end = self.time_end
                else:
                    time_end = sensor_geom.date_utc[-1]
                self.obs_dct[m] = ObsTable()
                self.obs_dct[m].geom = satgeo.SensorGeometry()
                self.obs_dct[m].geom.date_utc = []
                self.obs_dct[m].geom.vza      = []
                self.obs_dct[m].geom.vaa      = []
                self.obs_dct[m].geom.sza      = []
                self.obs_dct[m].geom.saa      = []
                for i,d in enumerate(sensor_geom.date_utc):
                    if d<time_start:
                        continue
                    elif d>time_end:
                        break
                    else:
                        self.obs_dct[m].geom.date_utc.append(d)
                        self.obs_dct[m].geom.vza.append(sensor_geom.vza[i])
                        self.obs_dct[m].geom.vaa.append(sensor_geom.vaa[i])
                        self.obs_dct[m].geom.sza.append(sensor_geom.sza[i])
                        self.obs_dct[m].geom.saa.append(sensor_geom.saa[i])
                msg = "after applying temporal restriction to time_start={} time_end={}".format(
                    time_start, time_end)
                msg += " remaining ngeom={}".format(len(self.obs_dct[m].geom.date_utc))
                FileLogger.info(msg)


    def load_obs_csv(self, csv_file, date_fmt="%Y/%m/%d %H:%M"):
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

        #-- first 7 columns are always:date, vza, vaa, sza, saa, lat, lon

        if ncol==9: #-- assume S1b sensor
            msg = "start reading S1 observations..."
            FileLogger.info(msg)
            #   date, vza, vaa, sza, saa, lat, lon, vh, vv
            vh_lst = []
            vv_lst = []
            self.obs_dct['S1b'] = ObsTable()
            self.obs_dct['S1b'].geom = satgeo.SensorGeometry()
            #-- abreviate
            sat_geom = self.obs_dct['S1b'].geom
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
                sat_geom.date_utc.append(act_date)
                sat_geom.vza.append( float(obs_data[i,1]) )
                sat_geom.vaa.append( float(obs_data[i,2]) )
                sat_geom.sza.append( float(obs_data[i,3]) )
                sat_geom.saa.append( float(obs_data[i,4]) )
                #-- VH,VV in 0-indexed columns 7,8:
                vh_lst.append( float(obs_data[i,7]) )
                vv_lst.append( float(obs_data[i,8]) )
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
            msg = "assuming {} [dB] S1 observational uncertainty.".format(self.s1_unc_db)
            FileLogger.info(msg)
            #-- XX_db = XX_db(XX)  = 10*log10(XX)
            #-- XX    = XX(XX_db)  = 10**(XX_db/10)
            #
            # for the uncertainty in linear/raw unit we apply conservative estimation:
            # 2*sXX = [ XX(XX_db+sXX_db) - XX(XX_db-sXX_db) ] (XX=VH,VV)
            #       = [ XX(XX_db)*10**(sXX_db/10.) - XX(XX_db)*10**(-sXX_db/10.)]
            #       = XX(XX_db)*[10**(sXX_db/10.) - 10**(-sXX_db/10.)]
            #       = XX * [10**(sXX_db/10.) - 10**(-sXX_db/10.)]
            svh = 0.5*vh*(10**(self.s1_unc_db/10.) - 10**(-1*self.s1_unc_db/10.))
            svv = 0.5*vv*(10**(self.s1_unc_db/10.) - 10**(-1*self.s1_unc_db/10.))
            msg = "determined VH uncertainty in linear units, min/max={}/{}".format(
                svh.min(), svh.max())
            FileLogger.info(msg)
            msg = "determined VV uncertainty in linear units, min/max={}/{}".format(
                svv.min(), svv.max())
            FileLogger.info(msg)
            #-- 
            nt_use = len(sat_geom.date_utc)
            self.obs_dct['S1b'].data = np.empty((nt_use,2), dtype=np.float64) #-- 'VH','VV'
            self.obs_dct['S1b'].data[:,0] = vh
            self.obs_dct['S1b'].data[:,1] = vv
            self.obs_dct['S1b'].dataunc = np.empty((nt_use,2), dtype=np.float64)
            self.obs_dct['S1b'].dataunc[:,0] = svh
            self.obs_dct['S1b'].dataunc[:,1] = svv
            #-- logging
            msg = "...reading S1 observations DONE"
            FileLogger.info(msg)
        else:
            #-- logging
            msg = "start reading S2 observations..."
            FileLogger.info(msg)
            #   date, vza, vaa, sza, saa, lat, lon, BRF1,...,BRF13
            self.obs_dct['S2a'] = ObsTable()
            self.obs_dct['S2a'].geom = satgeo.SensorGeometry()
            #-- abreviate
            sat_geom = self.obs_dct['S2a'].geom
            sat_geom.date_utc = []
            sat_geom.vza      = []
            sat_geom.vaa      = []
            sat_geom.sza      = []
            sat_geom.saa      = []
            brf_lst = [ [] for i in xrange(13) ] #-- prepare lists for 13 BRF bands
            for i,act_date in enumerate(date_lst):
                if act_date<time_start:
                    continue
                elif act_date>time_end:
                    break
                sat_geom.date_utc.append(act_date)
                sat_geom.vza.append( float(obs_data[i,1]) )
                sat_geom.vaa.append( float(obs_data[i,2]) )
                sat_geom.sza.append( float(obs_data[i,3]) )
                sat_geom.saa.append( float(obs_data[i,4]) )
                #-- BRFs start at 0-indexed column 7:
                for ib in xrange(13):
                    icol = ib+7
                    brf_lst[ib].append( float(obs_data[i, icol]) )
            #--
            nt_use = len(sat_geom.date_utc)
            brf_data = np.empty((nt_use,13), dtype=np.float64) #-- BRF1-13
            for ib in xrange(13):
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
            #-- set into structure
            self.obs_dct['S2a'].data    = brf_data
            self.obs_dct['S2a'].dataunc = brf_dataunc
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

        #-- store scheduled dates
        self.schedule_dct['date_utc'] = date_lst

        #-- determine
        #   - illumination-view geometry
        #   - simulation type (which Sentinel)
        #
        for i,d in enumerate(date_lst):
            #-- mission observation times
            for m,obs_pair in self.obs_dct.iteritems():
                geom = obs_pair.geom
                try:
                    im = geom.date_utc.index(d)
                    self.schedule_dct['sentinel'].append(m)
                    if m in ['S1a','S1b']:
                        self.schedule_dct['sim_typ'].append(1)
                    else:
                        self.schedule_dct['sim_typ'].append(2)
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

    def setprior_generic_agriculture(self):
        """
        Function to set a generic state vector for an agricultural site and
        the given temporal schedule. Following values are applied uniformily
        at all time points:
        LAI          : value=1.5, unc=5.  [m2/m2]
        canopy-height: value=1.5, unc=2.  [m]
        Soilmoisture : value=0.3, unc=0.2 [m3/m3]
        """

        #-- number of time-points
        npts = self.get_npts()

        #-- LAI,Canopy-Height,Soil-Moisture
        self.prstate     = np.empty((3,npts), dtype=np.float64)
        self.prstate_unc = np.empty((3,npts), dtype=np.float64)
        #-- LAI
        self.prstate[0,:]     = 2.5
        self.prstate_unc[0,:] = 5.
        #-- canopy-height
        self.prstate[1,:]     = 1.5
        self.prstate_unc[1,:] = 2.
        #-- soil moisture (volumetric)
        self.prstate[2,:]     = 0.3
        self.prstate_unc[2,:] = 0.2


    def setprior_synthetic(self):
        """
        Function that sets the default state vector shipped with the signature simulator.
        """

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


    def setprior_unc_default(self):
        """
        Function that sets the default uncertainty of the state vector.
        """
        if not self.__dict__.has_key('prstate') or self.prstate is None:
            msg = "internal error, prior state does not yet exist!"
            FileLogger.fatal(msg)
            raise RuntimeError(msg)

        self.prstate_unc = np.empty(self.prstate.shape, dtype=np.float64)

        #-- for the uncertainty we apply the settings
        #   already used within the prototype tool
        lai_relunc     = 0.5
        lai_uncfloor   = 3.0
        canht_relunc   = 0.05
        canht_uncfloor = 0.05
        sm_relunc      = 0.5
        sm_uncfloor    = 0.05
        self.prstate_unc[0,:] = np.maximum(self.prstate[0,:]*lai_relunc,   lai_uncfloor)
        self.prstate_unc[1,:] = np.maximum(self.prstate[1,:]*canht_relunc, canht_uncfloor)
        self.prstate_unc[2,:] = np.maximum(self.prstate[2,:]*sm_relunc,    sm_uncfloor)


    def setprior_csv(self, csv_file):
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
            # print "MVMV::nearest={} idx={} timedelt={}".format(
            #     state_inst.date_utc[idx], idx, timedelt)
            #-- LAI
            self.prstate[0,i]     = state_inst.lai[idx]
            #-- canopy-height
            self.prstate[1,i]     = state_inst.can_height[idx]
            #-- SM
            self.prstate[2,i]     = state_inst.soil_moisture[idx]


    def setprior_jules(self, jules_ncfile):
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
            # print "MVMV::nearest={} idx={} timedelt={}".format(
            #     state_inst.date_utc[idx], idx, timedelt)
            #-- LAI
            self.prstate[0,i]     = state_inst.lai[idx]
            #-- canopy-height
            self.prstate[1,i]     = state_inst.can_height[idx]
            #-- SM
            self.prstate[2,i]     = state_inst.soil_moisture[idx]


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
    ncvar.setncattr('long_name','state type')
    ncvar.setncattr('comment', 'integer value are to be bit-interpreted')
    ncvar.setncattr('nobits_set',  'time-point with other state')
    ncvar.setncattr('bit0_is_set', 'time-point for S1 simulation')
    ncvar.setncattr('bit1_is_set', 'time-point for S2 simulation')
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
    #-- 100% relative uncertainty
    ncvar[:] = retr_setup.lai_coeff*retr_setup.lai_coeff_relunc

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
    comment = "uncertainty was derived as {} [%] relative uncertainty ".format(
        100.*retr_setup.lai_relunc)
    comment += "and an uncertainty floor value of {} was applied.".format(
        retr_setup.bounds_dct['lai_lobound'])
    ncvar.setncattr('comment', comment)
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
    comment = "uncertainty was derived as {} [%] relative uncertainty ".format(
        100.*retr_setup.canht_relunc)
    comment += "and an uncertainty floor value of {} was applied.".format(
        retr_setup.bounds_dct['canht_lobound'])
    ncvar.setncattr('comment', comment)
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
    comment = "uncertainty was derived as {} [%] relative uncertainty ".format(
        100.*retr_setup.sm_relunc)
    comment += "and an uncertainty floor value of {} was applied.".format(
        retr_setup.bounds_dct['sm_lobound'])
    ncvar.setncattr('comment', comment)
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
    s1_table   = retr_setup.obs_dct['S1b']
    timepts    = s1_table.geom.date_utc
    npts       = len(timepts)
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
    s2_table   = retr_setup.obs_dct['S2a']
    timepts    = s2_table.geom.date_utc
    npts       = len(timepts)
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
    state_ftn = ".true." if retr_setup.use_state else ".false."
    state_cmt = "! if .true. then J_m included in Eq. 1.5"
    grad_tol  = retr_setup.gtol
    grad_cmt  = "! stopping criterion for minimisation: rel. reduction in gradient norm"
    pr_pert   = retr_setup.prior_pert
    pr_cmt    = "! initial guess will be a perturbed prior"
    try:
        fp = open(act_outname,'w')
    except IOError:
        msg = "file ***{}*** could not be opened for writing".format(act_outname)
        FileLogger.fatal(msg)
        return
    
    fp.write("&RETRCTL" + '\n')
    fp.write("  retr_use_prior_term = {} {}".format(prior_ftn, prior_cmt) + '\n')
    fp.write("  retr_use_state_term = {} {}".format(state_ftn, state_cmt) + '\n')
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
    if not os.path.exists(site_nml):
        msg = "specified namelist file ***{}*** does not exist.".format(site_nml)
        FileLogger.fatal(msg)
        return
    else:
        if site_nml!='site.nml':
            #-- 'site.nml' must exist in current working directory !
            msg = "Copying default namelist file ***{}*** to current working directory".format(site_nml)
            shutil.copyfile(site_nml, 'site.nml')

        #-- load site file
        msg = "START reading namelist ***{}***...".format(site_nml)
        FileLogger.info(msg)
        nml_dct  = f90nml.read(site_nml)
        #-- get dictionary for the site
        site_dct = nml_dct['site_params']
        #-- DEBUG
        # for k,v in site_dct.iteritems():
        #     print "{} => {}".format(k,v)
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

    #-- set missions
    if options.nos1:
        mission_lst = ['S2a']
    elif options.nos2:
        mission_lst = ['S1b']
    else:
        mission_lst = ['S1b', 'S2a']
    # logging
    msg = "selected mission_lst: {}".format(mission_lst)
    FileLogger.info(msg)

    #-- retrieval system control
    setup_dct['use_prior'] = options.use_prior
    setup_dct['use_state'] = options.use_state
    setup_dct['prior_pert'] = 0.25

    #-- temporal range
    if options.time_start==None:
        setup_dct['time_start'] = TIME_START_SYN
    else:
        setup_dct['time_start'] = datestr_parse(options.time_start)
    if options.time_end==None:
        setup_dct['time_end'] = TIME_END_SYN
    else:
        setup_dct['time_end'] = datestr_parse(options.time_end)
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

    #-- load satellite overpasses using pyorbital interface
    retrieval_setup.load_geom_synthetic(mission_lst)

    #-- setup schedule (merging the missions and extra target points)
    retrieval_setup.setup_common_schedule()

    #-- setup default prior
    if options.use_generic_prior:
        retrieval_setup.setprior_generic_agriculture()
    else:
        retrieval_setup.setprior_synthetic()

    #-- setup default prior uncertainty
    retrieval_setup.setprior_unc_default()

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
    if not os.path.exists(site_nml):
        msg = "specified namelist file ***{}*** does not exist.".format(site_nml)
        FileLogger.fatal(msg)
        return
    else:
        if site_nml!='site.nml':
            #-- 'site.nml' must exist in current working directory !
            msg = "Copying default namelist file ***{}*** to current working directory".format(site_nml)
            if options.outdir!=None:
                mkdirp_smart(options.outdir)
                shutil.copyfile(site_nml, os.path.join(options.outdir,'site.nml'))
            else:
                shutil.copyfile(site_nml, 'site.nml')

        #-- load site file
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
    setup_dct['use_state'] = options.use_state
    if options.s1_unc!=None:
        setup_dct['s1_unc'] = options.s1_unc
    if options.s2_relunc!=None:
        setup_dct['s2_relunc'] = options.s2_relunc
    if options.ctlvec_relunc!=None:
        setup_dct['lai_coeff_relunc'] = options.ctlvec_relunc[0]
        setup_dct['lai_relunc']       = options.ctlvec_relunc[1]
        setup_dct['canht_relunc']     = options.ctlvec_relunc[2]
        setup_dct['sm_relunc']        = options.ctlvec_relunc[3]

    setup_dct['gtol']      = options.gtol
    setup_dct['prior_pert'] = 0.

    #-- temporal range
    if options.time_start==None:
        setup_dct['time_start'] = options.time_start
    else:
        setup_dct['time_start'] = datestr_parse(options.time_start)
    if options.time_end==None:
        setup_dct['time_end'] = options.time_end
    else:
        setup_dct['time_end'] = datestr_parse(options.time_end)

    #-- instantiate RetrievalSetup
    retrieval_setup = RetrievalSetup(**setup_dct)

    #-- read observations
    obs_s1 = options.obs_s1
    obs_s2 = options.obs_s2
    no_obs_given = obs_s1==None and obs_s2==None
    if no_obs_given:
        msg = "No observational input files specified by user. Cannot continue!"
        FileLogger.fatal(msg)
        return
    if obs_s1!=None:
        retrieval_setup.load_obs_csv(obs_s1)
    if obs_s2!=None:
        retrieval_setup.load_obs_csv(obs_s2)

    #-- setup schedule
    #   NOTE: this merges the different missions and extra target points
    #         into a consecutive time-series
    retrieval_setup.setup_common_schedule()

    #-- setup default prior (state-variables)
    #   NOTE: this must be done *AFTER* the schedule is completed!
    #
    states_file = options.states_file
    if states_file==None:
        msg = "No states file specified by user, applying generic state suitable for agriculture!"
        FileLogger.info(msg)
        retrieval_setup.setprior_generic_agriculture()
    else:
        basename    = os.path.basename(states_file)
        if os.path.splitext(basename)[1]=='.nc':
            msg = "Prior state information will be read from ***{}***".format(states_file)
            FileLogger.info(msg)
            retrieval_setup.setprior_jules(states_file)
            msg = "...reading prior DONE"
            FileLogger.info(msg)
            msg = "! ! ! Will apply default uncertainties only - this must be changed ! ! !"
            FileLogger.warn(msg)
            retrieval_setup.setprior_unc_default()
        elif os.path.splitext(basename)[1]=='.csv':
            msg = "Prior state information will be read from ***{}***".format(states_file)
            FileLogger.info(msg)
            retrieval_setup.setprior_csv(states_file)
            msg = "...reading prior DONE"
            FileLogger.info(msg)
            msg = "! ! ! Will apply default uncertainties only - this must be changed ! ! !"
            FileLogger.warn(msg)
            retrieval_setup.setprior_unc_default()
        else:
            msg = "Unrecognised format of states file ***{}***. Cannot continue!".format(states_file)
            FileLogger.fatal(msg)
            return

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
        # aparser.add_argument( '--outname',
        #                       help="""name of generated file""" )
        aparser.add_argument( '--outdir',
                              help="""whether generated file(s) should be located in extra directory""" )
    # ---_add_output_options---

    def _add_common_options(aparser):
        aparser.add_argument( '--time_start',
                              help="""restrict retrieval system to time-points after this time (fmt: YYYYMMDDTHH:MM). NOTE, settings given here will overrule specifications given by options --target_select or --target_schedule!""" )
        aparser.add_argument( '--time_end',
                              help="""restrict retrieval system to time-points before this time (fmt: YYYYMMDDTHH:MM). NOTE, settings given here will overrule specifications given by options --target_select or --target_schedule!""" )
        aparser.add_argument( '--target_select',
                              metavar=('target_tmin','target_tmax','target_tdelta'),
                              nargs=3,
                              help="selection of series of equally distributed extra target " \
                              "time points, must be comma separated triple of tmin,tmax,tdelta " \
                              "with tmim,tmax given in format " \
                              "YYYYMMDDTHH:MM and tdelta specified as Xd or Xh with integer X. " \
                              "ATTENTION: If also a target schedule file was passed via option --target_schedule settings made here will be ignored!" )
        aparser.add_argument( '--target_schedule',
                              help="""specify name of a target schedule file from which extra time points will be read. The format of this file must be single time-point per line formatted as 'YmdTH:M:S'.""" )
        

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
                          default=os.path.join(ss_dir_path,'site.nml'),
                          help="""specification namelist file for selected site. If none is given a default namelist file suitable for the Wallerfing site will be used.""" )
    xparser.add_argument( '--use_generic_prior',
                          action='store_true',
                          help="""whether a generic prior suitable for agriculture sites should be used (By default a default prior shipped with the signature simulator is applied)""" )
    xparser.add_argument( '--nos1',
                          action='store_true',
                          help="""exclude points with S1 overpass time when generating the retrieval config file""" )
    xparser.add_argument( '--nos2',
                          action='store_true',
                          help="""exclude points with S2 overpass time when generating the retrieval config file""" )
    xparser.add_argument( '--no_use_prior',
                          dest='use_prior',
                          action='store_false',
                          default=True,
                          help="""whether to discard the prior information for the retrieval""" )
    xparser.add_argument( '--no_use_state',
                          dest='use_state',
                          action='store_false',
                          help="""whether to discard the dynamical model on the state vector for the retrieval""" )
    _add_output_options(xparser)

    #-------------------
    #      p r e _ g e n e r a l
    #
    xparser = subparsers.add_parser( 'pre_general',
                                     help="""generation of all relevant control and input files for the running the retrieval system for a general case according to user specifications.""" )
    #-- add common options
    _add_common_options(xparser)
    xparser.add_argument( '--site_nml',
                          default=os.path.join(ss_dir_path,'site.nml'),
                          help="""specification namelist file for selected site. If none is given a default namelist file suitable for the Wallerfing site will be used.""" )
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
    xparser.add_argument( '--no_use_state',
                          dest='use_state',
                          action='store_false',
                          help="""whether to discard the dynamical model on the state vector for the retrieval""" )
    xparser.add_argument( '--gtol',
                          type=float,
                          default=1.e0,
                          help="""gradient norm used as stopping-criterion within the minimisation routine.""" )
    xparser.add_argument( '--dynmodunc_inifile',
                          help="""ini file from where uncertainty specification for the 3 states in the dynamical model are read.""" )
    xparser.add_argument( '--s1_unc',
                          type=float,
                          help="""user-specified absolute uncertainty of S1 observations [dB]""" )
    xparser.add_argument( '--s2_relunc',
                          type=float,
                          help="""user-specified relative uncertainty of S2 observations""" )
    xparser.add_argument( '--ctlvec_relunc',
                          type=float,
                          nargs=4,
                          metavar=('laicoeff_relunc','lai_relunc','canht_relunc','sm_relunc'),
                          help="""user-specified relative uncertainty of controlvector components (uniform in time)""" )
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
