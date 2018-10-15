import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.mlab as mlab
import datetime as dt


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


def plot_var(var, dir='ret_code', site='ita', point='sg01', axes=None):
    prior_val = nc.Dataset(dir+'/controlvector_prior.nc', 'r')
    post = nc.Dataset(dir+'/controlvector_post.nc', 'r')
    #field_laican = mlab.csv2rec('/home/users/if910917/projects/ret_tool_subprocess/state_files/mni_lai_canht_field_508_'
    #                     'mean.csv', comments='%')
    if axes is None:
        fig, ax = plt.subplots()
        ret_val = fig
    else:
        ax = axes
        ret_val = ax
    var_dict = {'lai': r'Leaf area index (m$^{2}$ m$^{-2}$)', 'canht': 'Canopy height (m)',
                'sm': r'Soil moisture (m$^{3}$ m$^{-3}$)'}
    sat_times = nc.num2date(post.variables['time'][:], post.variables['time'].units)
    ax.plot(sat_times, prior_val.variables[var][:], '+', label='Prior')
    #ax.plot(times, np.array([np.nan]*len(times)), 'X')
    #else:
    #    field_times = field_laican['_date'][:]
    #    field_ob = field_laican[var][:]
    ax.errorbar(sat_times[post.variables['sim_typ'][:] == 9], post.variables[var][post.variables['sim_typ'][:] == 9],
                #yerr= post.variables[var+'_unc'][post.variables['sim_typ'][:] == 1],
                fmt='o', label='Retrieval output S1')
    ax.errorbar(sat_times[post.variables['sim_typ'][:] == 34], post.variables[var][post.variables['sim_typ'][:] == 34],
                #yerr=post.variables[var+'_unc'][post.variables['sim_typ'][:] == 2],
                fmt='o', label='Retrieval output S2')
    if var == 'sm':
        field_sm = mlab.csv2rec(
            '/home/users/if910917/projects/ret_tool_subprocess/state_files/' + site + '_sm_' + point + '.csv',
            comments='%', converterd={0: lambda x: dt.datetime.strptime(x, '%d/%m/%Y %H:%M')})
        t_idx = np.array([find_nearest(field_sm['date'], x)[1] for x in sat_times])
        print t_idx
        field_times = field_sm['date'][t_idx]
        field_ob = field_sm['sm'][t_idx]
        ax.plot(field_times, field_ob, '*', label='Field obs')
    ax.set_xlabel('Date')
    ax.set_ylabel(var_dict[var])
    if axes is None:
        fig.autofmt_xdate()
    plt.legend(frameon=True, fancybox=True, framealpha=0.5)
    return ret_val


def save_stats(var, dir, site='ita', point='sg01'):
    pri = nc.Dataset(dir+'/controlvector_prior.nc', 'r')
    post = nc.Dataset(dir+'/controlvector_post.nc', 'r')
    field_sm = mlab.csv2rec(
        '/home/users/if910917/projects/ret_tool_subprocess/state_files/' + site + '_sm_' + point + '.csv',
        comments='%', converterd={0: lambda x: dt.datetime.strptime(x, '%d/%m/%Y %H:%M')})
    sat_times = nc.num2date(post.variables['time'][:], post.variables['time'].units)
    if var == 'sm':
        t_idx = np.array([find_nearest(field_sm['date'], x)[1] for x in sat_times])
        #print t_idx
        field_times = field_sm['date'][t_idx]
        field_ob = field_sm['sm'][t_idx]
        post_obs = post.variables[var][:]
        pri_obs = pri.variables[var][:]

    innov = [((field_ob[i] - np.mean(field_ob)) - (post_obs[i] - np.mean(post_obs)))**2 for i in xrange(len(post_obs))]
    pos_ubrmse = np.sqrt(np.sum(innov) / len(post_obs))
    innov = [((field_ob[i] - np.mean(field_ob)) - (pri_obs[i] - np.mean(pri_obs)))**2 for i in xrange(len(post_obs))]
    pri_ubrmse = np.sqrt(np.sum(innov) / len(post_obs))

    pos_corrc = np.corrcoef(post_obs, field_ob)[0,1]
    pri_corrc = np.corrcoef(pri_obs, field_ob)[0, 1]

    innov = [(post_obs[i] - field_ob[i])**2 for i in xrange(len(post_obs))]
    pos_rmse = np.sqrt(np.sum(innov) / len(post_obs))
    innov = [(pri_obs[i] - field_ob[i])**2 for i in xrange(len(pri_obs))]
    pri_rmse = np.sqrt(np.sum(innov) / len(pri_obs))
    pri.close()
    post.close()
    return pri_rmse, pos_rmse, pri_ubrmse, pos_ubrmse, pri_corrc, pos_corrc