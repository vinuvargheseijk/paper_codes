import numpy as np

f = 2.0

def pSpike(rate, T):
    spike_times = []
    t = 0
    while t < T:
       interval = -np.log(np.random.rand()) / rate
       t = t + interval
       if t < T:
          spike_times.append(t)
    return spike_times

def sine_series(T, rate, dt):
    time_series = np.arange(0, T, dt)
    trate = rate * np.sin((f * np.pi * 1 * time_series + 1) / 2.0)
    return trate
    
def IpSpike(rate, T, dt):
    trate = sine_series(T, rate, dt)
    n_bins = len(trate)
    spikes = np.random.rand(n_bins) < trate * dt
    spike_times = np.nonzero(spikes)[0] * dt
    print("CV :", calc_cv(np.diff(spike_times)))
    return spike_times

def calc_cv(isi):
    return np.std(isi) / np.mean(isi)
