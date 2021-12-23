from model_func import *

cols = np.array([[246, 56, 29, 195], [82, 161, 242, 255], [126, 126, 126, 233]])/255


######## simulate one node active (fig.1b)
input_time1 = 50
n_units = 5
tau_i = 32
n_times = 1000+input_time1*2
all_alpha_time = run_netw_simul(n_times = n_times, n_units=n_units, tau_i = tau_i,
		input_time1=input_time1, input_nodes1 = np.array([0]), impulse_amp1 = 5, impulse_dur1 = 10)

timeps = np.arange(-input_time1, n_times-input_time1)

# plot time course of the simulation for all units
fig, ax = pl.subplots(1, 1)
for unit_i in range(3):
	ax.plot(timeps, all_alpha_time[unit_i, 0, :], color=cols[unit_i])
ax.set_ylim([-1, 13])
ax.set_xlim([-input_time1, n_times-input_time1*2])
ax.set_xticks(np.arange(0, n_times, 200))
ax.grid()
fig.set_size_inches([6.88, 2.17])
ax.set_title('One unit active in a 5-unit network (Fig. 1b)')
pl.tight_layout()



######## simulate two pairs of nodes active (fig.1b)
input_time1 = 50
input_time2 = 100
n_units = 5
tau_i = 32
n_times = 1000+input_time1*2
all_alpha_time = run_netw_simul(n_times = n_times, n_units=n_units, tau_i = tau_i,
		input_time1=input_time1, input_nodes1 = np.array([0, 1]), impulse_amp1 = 5, impulse_dur1 = 10,
		input_time2=input_time2, input_nodes2 = np.array([2, 3]), impulse_amp2 = 5, impulse_dur2 = 10)

# plot time course of the simulation for all units
fig, axs = pl.subplots(n_units, 1)
for i in np.arange(n_units):
	for unit_i in range(3):
		axs[i].plot(timeps, all_alpha_time[unit_i, i, :], color=cols[unit_i])
	axs[i].set_ylim([-1, 13])
	axs[i].set_xlim([-input_time1, n_times-input_time1*2])
	axs[i].set_xticks(np.arange(0, n_times, 200))
	axs[i].grid()

fig.set_size_inches([5.88, 6.17])
pl.tight_layout()
pl.suptitle('Fig. 3')



######## simulate three nodes active out-of-phase (fig.4a)
input_time1 = 50
n_units = 5
tau_i = 32
n_times = 1000+input_time1*2
all_alpha_time = run_netw_simul(n_times = n_times, n_units=n_units, tau_i = tau_i,
		input_time1=input_time1, input_nodes1 = np.array([0]), impulse_amp1 = 5, impulse_dur1 = 10,
		input_time2=250, input_nodes2 = np.array([1]), impulse_amp2 = 5, impulse_dur2 = 10,
		input_time3=450, input_nodes3 = np.array([2]), impulse_amp3 = 5, impulse_dur3 = 10)

timeps = np.arange(-input_time1, n_times-input_time1)

# plot time course of the simulation for all units
fig, axs = pl.subplots(n_units, 1)
for i in np.arange(n_units):
	for unit_i in range(3):
		axs[i].plot(timeps, all_alpha_time[unit_i, i, :], color=cols[unit_i])
	axs[i].set_ylim([-1, 13])
	axs[i].set_xlim([-input_time1, n_times-input_time1*2])
	axs[i].set_xticks(np.arange(0, n_times, 200))
	axs[i].grid()

fig.set_size_inches([5.88, 6.17])
axs[0].set_title('Slower freq., 3 nodes active out-of-phase (Fig. 4a)')
pl.tight_layout()





######## simulate faster frequency and impossibility of three units active out-of-phase (fig.4b)
n_units = 5
tau_i = 20
n_times = 1000+input_time1*2
all_alpha_time = run_netw_simul(n_times = n_times, n_units=n_units, tau_i = tau_i,
		input_time1=50, input_nodes1 = np.array([0]), impulse_amp1 = 8, impulse_dur1 = 10,
		input_time2=200, input_nodes2 = np.array([1]), impulse_amp2 = 8, impulse_dur2 = 10,
		input_time3=335, input_nodes3 = np.array([2]), impulse_amp3 = 8, impulse_dur3 = 10)

timeps = np.arange(-input_time1, n_times-input_time1)

# plot time course of the simulation for all units
fig, axs = pl.subplots(n_units, 1)
for i in np.arange(n_units):
	for unit_i in range(3):
		axs[i].plot(timeps, all_alpha_time[unit_i, i, :], color=cols[unit_i])
	axs[i].set_ylim([-1, 13])
	axs[i].set_xlim([-input_time1, n_times-input_time1*2])
	axs[i].set_xticks(np.arange(0, n_times, 200))
	axs[i].grid()
fig.set_size_inches([5.88, 6.17])
axs[0].set_title('Faster freq., 3 nodes NOT active out-of-phase (Fig. 4b)')
pl.tight_layout()





