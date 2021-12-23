# model of binding through synchrony for working memory based on Pina et al., Plos Computational Biology (2018)

import numpy as np
import scipy
from matplotlib import pyplot as pl

def frate_func(x, beta):
	return np.sqrt(x/(1 - np.exp(-beta * x)))

def u_der(u, u_tilde, v_tilde, n_tilde, a_ee, a_ei, a_en, theta_e, input_sig, beta):
	return (-u + frate_func((a_ee * u_tilde) - (a_ei * v_tilde) + (a_en * n_tilde) - theta_e + input_sig, beta))

def v_der(v, u, n, a_ie, a_ii, a_in, theta_i, beta):
	return -v + frate_func((a_ie * u) - (a_ii * v) + (a_in * n) - theta_i, beta)

def n_der(n, a_n, u, p):
	return -n + a_n*(u**p)*(1-n)


def alpha_tilde(all_alphas, unit_i, c_z, N):
	other_alphas_inds = np.arange(N)
	other_alphas_inds = other_alphas_inds[other_alphas_inds != unit_i]

	other_alphas_sum = all_alphas[other_alphas_inds].sum()
	return (all_alphas[unit_i] + c_z * other_alphas_sum)/(1 + c_z * (N - 1))



def run_netw_simul(n_times=1000, n_units=5, tau_i = 24, tau_n = 244, c_e = .001,
	c_ei = .03, a_ee = 14, a_ei = 9, a_en = 4, theta_e = 6, a_ie = 20, a_ii = 8, a_in = .1,
	theta_i = 5, a_n = 2, beta = 1, p = 2,
	input_time1=None, input_nodes1 = None, impulse_amp1 = None, impulse_dur1 = None,
	input_time2=None, input_nodes2 = None, impulse_amp2 = None, impulse_dur2 = None,
	input_time3=None, input_nodes3 = None, impulse_amp3 = None, impulse_dur3 = None,
	input_time4=None, input_nodes4 = None, impulse_amp4 = None, impulse_dur4 = None):

	all_alpha_time = np.zeros(shape = [3, n_units, n_times])
	# need to do this otherwise if we keep all to 0 the network starts oscillating by itself...
	all_alpha_time[:, :, 0] = np.array([.093, .2, .02])[:, np.newaxis]
	all_alpha_tilde_time = np.zeros(shape = [3, n_units, n_times])

	for t in np.arange(0, n_times - 1):
		if input_time1 != None:
			input_sig = np.zeros(n_units)
			if (t >= input_time1) & (t < (input_time1+impulse_dur1)):
				input_sig[input_nodes1] = impulse_amp1

		if input_time2 != None:
			# add out-of-phase pop
			if (t >= input_time2) & (t < (input_time2+impulse_dur2)):
				input_sig[input_nodes2] = impulse_amp2

		if input_time3 != None:
			# add out-of-phase pop
			if (t >= input_time3) & (t < (input_time3+impulse_dur3)):
				input_sig[input_nodes3] = impulse_amp3

		if input_time4 != None:
			# add out-of-phase pop
			if (t >= input_time4) & (t < (input_time4+impulse_dur4)):
				input_sig[input_nodes4] = impulse_amp4

		# compute tildes
		for alpha_ind in np.arange(3):
			if alpha_ind == 1:
				c_z = c_ei
			else:
				c_z = c_e

			for unit_n in np.arange(n_units):
				all_alpha_tilde_time[alpha_ind, unit_n, t + 1] = alpha_tilde(all_alpha_time[alpha_ind, :, t], unit_n, c_z, n_units)

		for alpha_ind in np.arange(3):
			for unit_n in np.arange(n_units):
				u, v, n = all_alpha_time[:, unit_n, t]
				if alpha_ind == 0:
					u = all_alpha_time[alpha_ind, unit_n, t]
					u_tilde, v_tilde, n_tilde = all_alpha_tilde_time[:, unit_n, t]
					all_alpha_time[alpha_ind, unit_n, t + 1] = all_alpha_time[alpha_ind, unit_n, t] +\
						u_der(u, u_tilde, v_tilde, n_tilde, a_ee, a_ei, a_en, theta_e, input_sig[unit_n], beta)

				elif alpha_ind == 1:
					all_alpha_time[alpha_ind, unit_n, t + 1] = all_alpha_time[alpha_ind, unit_n, t] + v_der(v, u, n, a_ie, a_ii, a_in, theta_i, beta)/tau_i

				elif alpha_ind == 2:
					all_alpha_time[alpha_ind, unit_n, t + 1] = all_alpha_time[alpha_ind, unit_n, t] + n_der(n, a_n, u, p)/tau_n


	return all_alpha_time
