import numpy as np
import matplotlib.pyplot as mp
import csv

nights = ['171012', '171013', '171014','171017', '171018', '171019',
    '171020', '171021', '171022', '171023','171024','171026','171027','171028','171029','171102','171103','171104','171105','171106','171107']

stars = np.loadtxt('periods.txt', dtype=[('star', 'S10'), ('period', float)])

t0 = 8037.0

dtype = np.dtype([('frame_num', int), ('object', 'S10'),
    ('exptime', float), ('ra_off', float), ('dec_off', float),
    ('date_utc', 'S10'), ('time_utc', 'S10'), ('date_lt', 'S10'),
    ('time_lt', 'S10'), ('jd', float), ('epoch', int), ('ra', 'S10'),
    ('dec', 'S11'), ('airmass', float)])

objects = np.array([], dtype='S10')
times = np.array([], dtype=float)
night_num = np.array([], dtype='S6')



for ind, night in enumerate(nights):
    obs = np.loadtxt(night+'/obslog/obslog-nodark.txt', dtype=dtype, skiprows=1)
    objects = np.append(objects, obs['object'])
    times = np.append(times, obs['jd'])
    night_num = np.append(night_num, np.repeat(night, len(obs['object'])))

uniq_objects = np.unique(objects)
uniq_objects = np.delete(uniq_objects, np.argwhere(uniq_objects == 'twflat'))
#uniq_objects = np.delete(uniq_objects, np.argwhere(uniq_objects == 'DARK'))


fig, ax1 = mp.subplots(1,1)
star_num = np.arange(len(uniq_objects))
ax1.set_yticks(star_num)
ax1.set_yticklabels(uniq_objects)

colors=['r', 'b', 'g', 'm', 'c', 'k', 'xkcd:brown', 'xkcd:orange',
    'xkcd:purple', 'xkcd:teal', 'xkcd:lavender', 'xkcd:salmon',
    'xkcd:mustard', 'xkcd:mint','xkcd:sky','xkcd:maroon','xkcd:lime','xkcd:poop','deepskyblue','xkcd:grey','xkcd:tan']

for ind, obj in enumerate(uniq_objects):
#    print(obj)
    period = stars['period'][stars['star'] == obj]
    if np.isnan(period) == 1 : period = [0.5]

    # phase observations
    times_this_star = times[objects == obj]
    night_this_star = night_num[objects == obj]
    n_obs = len(times_this_star)
    print(obj)
#    print(n_obs/5)
    phase = np.mod((times_this_star - t0)/period, 1)
    for idx, night in enumerate(nights):

        phase_per_night = phase[night_this_star == night]
        y = np.repeat(ind, len(phase_per_night))
        ax1.plot(phase_per_night, y, 'o', color=colors[idx])



mp.show()
