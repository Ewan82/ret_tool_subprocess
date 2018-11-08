import numpy as np
import subprocess
import os
import sys
import itertools as itt
import matplotlib.pyplot as plt
import glob
import netCDF4 as nc
import plot_ret as pr


def find_nearest(array, value):
    """
    Find nearest value in an array
    :param array: array of values
    :param value: value for which to find nearest element
    :return: nearest value in array, index of nearest value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def run_ret_tool(out_dir, run_ret=True, **kwargs):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    print dir_path
    ret_tool_dic = {'time_start': '20170323', 'time_end': '20170720', 'site_nml': dir_path+'/site.nml', 'states_file':
                    False, 'obs_s1': False, 'obs_s2': False, 'no_use_prior': False, 'no_use_states': False, 'gtol':
                    1e0, 'dynmodunc_inifile': dir_path+'/mod_ini/dynmod_default.ini', 's1_unc': 0.4, 's2_relunc': 0.05,
                    's2_uncfloor': 0.001, 'ctlvec_relunc': [1.0, 0.5, 0.5, 0.5],
                    'ctlvec_uncfloor': [1e-6, 1e-3, 1e-3, 1e-3]}
    if kwargs is not None:
        for key, value in kwargs.iteritems():
            if value is not False:
                ret_tool_dic[key] = value
    if run_ret is True:
        subprocess.call(['cp', '-R', dir_path+'/ret_code', out_dir])
        os.chdir(out_dir)
        subprocess.call(['make', 'setup'])
        cmd = ['bin/rs_pre.py', 'pre_general']
        for key in ret_tool_dic.keys():
            if ret_tool_dic[key] is not False:
                if type(ret_tool_dic[key]) is str:
                    print 'str, '+key
                    cmd.append('--'+key)
                    cmd.append(ret_tool_dic[key])
                elif type(ret_tool_dic[key]) is float:
                    print 'float, '+key
                    cmd.append('--'+key)
                    cmd.append(str(ret_tool_dic[key]))
                elif type(ret_tool_dic[key]) is list:
                    cmd.append('--'+key)
                    for x in ret_tool_dic[key]:
                        print x
                        cmd.append(str(x))
                else:
                    cmd.append('--'+key)
        subprocess.call(cmd)
        subprocess.call(['make', 'retrieval'])
        subprocess.call(['make', 'mba'])
        os.chdir(dir_path)

    return ret_tool_dic


def save_plot_err(_dir, sname='mni_508_med_test.png', point='508_med', field_obs=True):
    #fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(14, 6),)
    #pr.plot_var('sm', _dir, axes=ax[0], point=point, field_obs_on=field_obs)
        #fig.savefig(_dir + '/sm.png', bbox_inches='tight')
    #pr.plot_var('lai', _dir, axes=ax[1], point=point, field_obs_on=field_obs)
        #fig.savefig(_dir + '/lai.png', bbox_inches='tight')
    #pr.plot_var('canht', _dir, axes=ax[2], point=point, field_obs_on=field_obs)
    #fig.autofmt_xdate()
    #fig.savefig(_dir + '/'+sname, bbox_inches='tight')

    #fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6),)
    #pr.plot_refl(_dir, point, axes=ax[0])
    #pr.plot_backscat(_dir, point, axes=ax[1])
    #fig.autofmt_xdate()
    #fig.savefig(_dir + '/r_b_'+sname, bbox_inches='tight')

    #fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 8),)
    #pr.plot_refl_mod(_dir, point, axes=ax[0])
    #pr.plot_backscat_mod(_dir, point, axes=ax[1])
    #fig.autofmt_xdate()
    #fig.savefig(_dir + '/r_b_mod_'+sname, bbox_inches='tight')

    sm_stats = pr.save_stats('sm', _dir, point=point)
    lai_stats = pr.save_stats('lai', _dir, point=point)
    canht_stats = pr.save_stats('canht', _dir, point=point)
    stats_file = _dir + '/stats_tayl.txt'
    lines = []
    lines.append('experiment stats \n')
    lines.append('state_var: prior rmse, post rmse, prior ubrmse, post ubrmse, prior corrcoef, post corrcoef, prior taylor, post taylor \n')
    lines.append('lai: '+str(lai_stats)+' \n')
    lines.append('sm: ' + str(sm_stats) + ' \n')
    lines.append('canht: ' + str(canht_stats) + ' \n')
    f = open(stats_file, 'w')
    for line in lines:
        f.write(line)
    f.close()
    return 'stats saved! :-)'


def make_exps_508med_arr():
    time_start = ['20170323']
    time_end = ['20170720']
    states_file = ['/home/users/if910917/projects/ret_tool_subprocess/state_files/retr_jules_prior.csv',
                   False]
    obs_s1 = ['/home/users/if910917/projects/ret_tool_subprocess/S1/mni_s1_508_med_hv.csv']  #, False]
    obs_s2 = ['/home/users/if910917/projects/ret_tool_subprocess/S2/mni_s2_508_med_b4b8.csv',]
              #'/home/users/if910917/projects/ret_tool_subprocess/S2/mni_s2_508_med_allbands.csv', False]
    #no_use_prior = [False, True]
    #no_use_states = [False, True]
    dynmodunc_inifile = ['/home/users/if910917/projects/ret_tool_subprocess/mod_ini/dynmod_default.ini',
                         #'/home/users/if910917/projects/ret_tool_subprocess/mod_ini/dynmod_tight.ini',
                         '/home/users/if910917/projects/ret_tool_subprocess/mod_ini/dynmod_vtight.ini']
    s1_unc = [0.4, 0.8, 1.6]
    s2_uncfloor = [0.001, 0.01, 0.02]
    s2_relunc = [0.05]
    ctlvec_relunc = [[1.0, 0.5, 0.5, 0.5], [0.1, 0.1, 0.1, 0.1]]
    ctlvec_uncfloor = [[1e-6, 1e-3, 1e-3, 1e-3], [0.05, 0.5, 0.1, 0.05],]  # [0.05, 3.0, 1.5, 0.1]]
    exp_list = [time_start, time_end, states_file, obs_s1, obs_s2, dynmodunc_inifile,
                s1_unc, s2_uncfloor, ctlvec_relunc, ctlvec_uncfloor, s2_relunc]
    it_exp_list = list(itt.product(*exp_list))
    return it_exp_list


def make_exps_508med_arr_test(point='508_med', date_start='20170323', date_end='20170720'):
    time_start = [date_start]  # ['20170324']
    time_end = [date_end]  # ['20170717']
    states_file = ['/home/users/if910917/projects/ret_tool_subprocess/state_files/mni_retr_jules_prior.csv',
                   False]
    obs_s1 = ['/home/users/if910917/projects/ret_tool_subprocess/S1/mni_s1_'+point+'_hv.csv',
              '/home/users/if910917/projects/ret_tool_subprocess/S1/mni_s1_' + point + '_vv.csv',
              '/home/users/if910917/projects/ret_tool_subprocess/S1/mni_s1_'+point+'_hvvv.csv',
              #'/home/users/if910917/projects/ret_tool_subprocess/S1/mni_s1_' + point + '_2017.csv',
              False
              ]
    obs_s2 = ['/home/users/if910917/projects/ret_tool_subprocess/S2/mni_s2_'+point+'_b4b8.csv',
              '/home/users/if910917/projects/ret_tool_subprocess/S2/mni_s2_'+point+'_b4b5b6b7b8b8a.csv',
              '/home/users/if910917/projects/ret_tool_subprocess/S2/mni_s2_' + point + '_allbands.csv',
              #'/home/users/if910917/projects/ret_tool_subprocess/S2/mni_s2_' + point + '_2017.csv',
              False
              ]
    #no_use_prior = [False, True]
    #no_use_states = [False, True]
    dynmodunc_inifile = [#'/home/users/if910917/projects/ret_tool_subprocess/mod_ini/dynmod_default.ini',
                         #'/home/users/if910917/projects/ret_tool_subprocess/mod_ini/dynmod_tight.ini',
                         #'/home/users/if910917/projects/ret_tool_subprocess/mod_ini/dynmod_test.ini',
                         '/home/users/if910917/projects/ret_tool_subprocess/mod_ini/dynmod_vtight.ini']
    s1_unc = [1.6]
    s2_uncfloor = [0.02]
    s2_relunc = [0.05]
    ctlvec_relunc = [[0.01, 0.5, 0.05, 0.5],]  # [0.01, 0.5, 0.05, 0.5]]
    ctlvec_uncfloor = [[0.001, 3.0, 0.1, 0.05],]  # [0.001, 3.0, 1.0, 0.05]]
    exp_list = [time_start, time_end, states_file, obs_s1, obs_s2, dynmodunc_inifile,
                s1_unc, s2_uncfloor, ctlvec_relunc, ctlvec_uncfloor, s2_relunc]
    it_exp_list = list(itt.product(*exp_list))
    exp_list_final = [a for a in it_exp_list if a[3:5] != (False, False)]
    return exp_list_final


def run_retr_exp(exp_no, exp_tup, point='508_med', run_ret=True, field_obs=True):
    if type(exp_tup[2]) is str:
        prior = 'jules'
    else:
        prior = 'flat'
    try:
        s1_str = exp_tup[3].split('_')[-1][:-4]
    except(AttributeError):
        s1_str = 'nos1'
    try:
        s2_str = exp_tup[4].split('_')[-1][:-4]
    except(AttributeError):
        s2_str = 'nos2'
    save_dir = '/export/cloud/nceo/users/if910917/s3_exps/mni/'+point+'/'+prior+'_'+s1_str+'_'+s2_str
    ret_dic = run_ret_tool(save_dir, run_ret=run_ret, time_start=exp_tup[0], time_end=exp_tup[1], states_file=exp_tup[2],
                           obs_s1=exp_tup[3], obs_s2=exp_tup[4], dynmodunc_inifile=exp_tup[5], s1_unc=exp_tup[6],
                           s2_uncfloor=exp_tup[7], ctlvec_relunc=exp_tup[8], ctlvec_uncfloor=exp_tup[9],
                           s2_relunc=exp_tup[10])
    save_plot_err(save_dir, sname='mni_'+point+'_'+prior+'_'+s1_str+'_'+s2_str+'.png', point=point, field_obs=field_obs)
    stats_file = save_dir + '/exp_' + str(exp_no) + '_setup.txt'
    lines = []
    for key in ret_dic.keys():
        lines.append(key+': '+str(ret_dic[key])+'\n')
    f = open(stats_file, 'w')
    for line in lines:
        f.write(line)
    f.close()
    return 'experiment '+str(exp_no)+' done! :-)'


def exp_run_setup(point='508_med', run_ret=True, field_ob=True, date_start=None, date_end=None):
    """
    Runs JULES for specified lat lon
    :param lat_lon: tuple containing latitude and longitude coordinate
    :return: na
    """
    exp_list = make_exps_508med_arr_test()
    run_dir = '/home/if910917/qSub_runMe/'
    for x in xrange(len(exp_list)):
        run_file = run_dir+'exp_'+str(x)+'_'+ point +'.bash'
        lines = []
        lines.append('cd ' + os.getcwd() + '\n')
        lines.append('module load python/canopy-1.7.2 \n')
        lines.append('python run_retr.py ' + str(x) + ' ' + point + ' ' + str(run_ret) + ' ' + str(field_ob) +
                     ' ' + date_start + ' ' + date_end + '\n')
        f = open(run_file, 'w')
        for line in lines:
            f.write(line)
        f.close()
    return 'exp list written!'


def find_good_exp(_dir):

    sm_stats = pr.save_stats('sm', _dir)
    if sm_stats[-1] > 0.45:
        print _dir+' SM is a gooden! '+str(sm_stats[-1])
    lai_stats = pr.save_stats('lai', _dir)
    if lai_stats[-1] > 0.8:
        print _dir+' LAI is a gooden! '+str(lai_stats[-1])
    canht_stats = pr.save_stats('canht', _dir)
    if canht_stats[-1] > 0.7:
        print _dir+' Canht is a gooden! '+str(canht_stats[-1])
    return None


if __name__ == "__main__":
    exp_list = make_exps_508med_arr_test(sys.argv[2], sys.argv[5], sys.argv[6])
    exp_no = int(sys.argv[1])
    point = sys.argv[2]
    run_ret = sys.argv[3] == 'True'
    field_ob = sys.argv[4] == 'True'
    run_retr_exp(exp_no, exp_list[exp_no], point, run_ret=run_ret, field_obs=field_ob)
    print 'ran experiment '+str(exp_no)
