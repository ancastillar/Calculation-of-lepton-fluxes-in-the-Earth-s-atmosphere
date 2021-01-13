import os, sys

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.integrate import trapz
from scipy import interpolate

from MCEq.core import MCEqRun
import CRFluxModels as pm
from mceq_config import config

from optparse import OptionParser

### options ###

parser = OptionParser()
parser.add_option("-d", "--day", dest="day", default="0", help="Day for AIRS data. From 0 to 364. Default: 0")
parser.add_option("-y", "--year", dest="year", default="2012", help="Year for AIRS data. Default: 2012")
parser.add_option("-m", "--model", dest="model", default="SIBYLL-2.3c", help="Hadronic interaction model: SIBYLL-2.1, SIBYLL-2.1, QGSJet-II-04, EPOS-LHC. Default: SIBYLL-2.3c")
parser.add_option("-e", "--emin", dest="emin", default="1.0", help="Minimum lepton energy in GeV. Default: 1.0")
parser.add_option("-p", "--prim", dest="prim", default="H3a", help="Primary spectrum: H3a, H4a, GST3, GST4, polygonato, polygonato-mod. Default: H3a")
(options, args) = parser.parse_args()

print("Options: ", options)

model = options.model
day = int(options.day) + 1
year = int(options.year)
emin = float(options.emin)
prim = options.prim

data_dir = '/data/user/dsoldin/MCEq/hadr_models/data/'
Aeffdir = '/data/user/dsoldin/seasonal_variations/data/effective_areas/'
AeffFile = Aeffdir + 'Aeff.dat'

coszen = np.linspace(0.05, 1., 20)


def get_season(day):
    if day<32:
    	return 'January'
    elif (day>=32 and day<60):
    	return 'February'
    elif (day>=60 and day<91):
    	return 'March'
    elif (day>=91 and day<121):
    	return 'April'
    elif (day>=121 and day<152):
    	return 'May'
    elif (day>=152 and day<182):
    	return 'June'
    elif (day>=182 and day<213):
    	return 'July'
    elif (day>=213 and day<244):
    	return 'August'
    elif (day>=244 and day<274):
    	return 'September'
    elif (day>=274 and day<305):
    	return 'October'
    elif (day>=305 and day<335):
    	return 'November'
    elif (day>=335 and day<366):
    	return 'December'
    else:
    	print 'day larger than year...'
    	print '...exit!'
    	sys.exit()
    	

def ReadAeff( infile, row ):
	infile = open(infile, 'r')
	dummy= []
	for line in infile:
		line = line.strip()
		columns = line.split()
		if row==10 or row==0:
			dummy.append(float(columns[row]))
		else:
			dummy.append(float(columns[row][:-1]))

	infile.close()
	return np.array(dummy)
	
### lepton production channels ###

channels = [ 'total_mu+', 'total_mu-', 'conv_mu+', 'conv_mu-', 'pi_mu+', 'pi_mu-', 'k_mu+', 'k_mu-', 'pr_mu+', 'pr_mu-',
			'total_numu', 'total_antinumu', 'conv_numu', 'conv_antinumu', 'pi_numu', 'pi_antinumu', 'k_numu', 'k_antinumu', 'pr_numu', 'pr_antinumu',
			'total_nue', 'total_antinue', 'conv_nue', 'conv_antinue', 'pi_nue', 'pi_antinue', 'k_nue', 'k_antinue', 'pr_nue', 'pr_antinue']
channels2 = [ 'obs_numu', 'obs_antinumu', 'obs_nue', 'obs_antinue' ]

spec = { 'H3a': [pm.HillasGaisser2012, 'H3a'], 'H4a': [pm.HillasGaisser2012, 'H4a'], 'GST3': [pm.GaisserStanevTilav, '3-gen'], 'GST4': [pm.GaisserStanevTilav, '4-gen'], 'polygonato': [pm.PolyGonato, False], 'polygonato-mod': [pm.PolyGonato, False] }


##############################
#          run MCEq          #
##############################
z0=1
for zen in coszen:

	mag = 3.
	
	if emin ==1.0:
		outfile = data_dir + '/%(p)s/%(z)1.2f/%(m)s_%(y)i_%(d)03i.p' % dict(p=prim, z=zen, m=model, y=year, d=day)
	else:
		outfile = data_dir + '/%(p)s/%(z)1.2f/%(m)s_%(e)3.0f_%(y)i_%(d)03i.p' % dict(p=prim, z=zen, e=emin, m=model, y=year, d=day)
	
	print outfile
	### spectrum ###
	
	config["h_obs"] = 2834.

	pmodel = (spec[prim][0], spec[prim][1])  
	mceq_run = MCEqRun(**dict(interaction_model=model,
							primary_model=pmodel,
							theta_deg=np.degrees(np.arccos(zen)),
							**config))

	dens_args = ('SouthPole', get_season(day), True, {"year": year ,"day": day} )
	mceq_run.set_density_model(('AIRS', dens_args))
	e_grid = mceq_run.e_grid
	
	mceq_run.solve()
	
	if not z0%2==0:
		E_Aeff = ReadAeff(AeffFile, 0)
		Aeff = ReadAeff(AeffFile, int(z0/2.))
	else:
		E_Aeff = 0
		Aeff = 0

	total = {}
	total_Aeff = {}
	total_sum = {}
	total_sum_Aeff = {}
	for id in channels:
		y = mceq_run.get_solution(id, mag)
		total[id] = y
		y0 = mceq_run.get_solution(id, 0)
		total_sum[id] = trapz(y0[e_grid>=emin], e_grid[e_grid>=emin])
		if not z0%2==0:
			spec_spline = interpolate.UnivariateSpline(e_grid[e_grid>=emin], y[e_grid>=emin], s=0)
			spec_Aeff = spec_spline(10**E_Aeff) * 100.**2 * Aeff * 2.*np.pi 
			total_Aeff[id] = spec_Aeff
			total_sum_Aeff[id] = trapz(spec_Aeff[10**E_Aeff>=emin], E_Aeff[10**E_Aeff>=emin])
		
	
	mceq_run.set_obs_particles(["mu+", "mu-", "pi_mu+", "pi_mu-", "k_mu+", "k_mu-", "pr_mu+", "pr_mu-"])
	mceq_run.solve()

	for id in channels2:
		y = mceq_run.get_solution(id, mag)
		total[id] = y
		y0 = mceq_run.get_solution(id, 0)
		total_sum[id] = trapz(y0[e_grid>=emin], e_grid[e_grid>=emin])
		if not z0%2==0:
			spec_spline = interpolate.UnivariateSpline(e_grid[e_grid>=emin], y[e_grid>=emin], s=0)
			spec_Aeff = spec_spline(10**E_Aeff) * 100.**2 * Aeff * 2.*np.pi 
			total_Aeff[id] = spec_Aeff
			total_sum_Aeff[id] = trapz(spec_Aeff[E_Aeff>=emin], E_Aeff[E_Aeff>=emin])


	### production depth ###

	config["h_obs"] = 1.
 
	mceq_run = MCEqRun(**dict(interaction_model=model,
							primary_model=pmodel,
							theta_deg=np.degrees(np.arccos(zen)),
							**config))

	mceq_run.set_density_model(('AIRS', dens_args))

	Xvec = np.logspace(-2, np.log10(mceq_run.density_model.max_X),1000)
	Hvec = mceq_run.density_model.s_lX2h(np.log(Xvec))/1e2

	mceq_run._calculate_integration_path(int_grid=Xvec, grid_var='X',force=True)
	mceq_run.solve(int_grid=Xvec, grid_var="X")

	total_X = {}
	total_X_Aeff = {}
	for id in channels:
		out = np.empty_like(Xvec)
		out_Aeff = np.empty_like(Xvec)
		for idx, X in enumerate(Xvec):
			y = (mceq_run.get_solution(id, 0, grid_idx=idx))
			out[idx] = trapz(y[e_grid>=emin], e_grid[e_grid>=emin])
			if not z0%2==0:
				spec_spline = interpolate.UnivariateSpline(e_grid[e_grid>=emin], y[e_grid>=emin], s=0)
				spec_Aeff = spec_spline(10**E_Aeff) * 100.**2 * Aeff * 2.*np.pi 
				out_Aeff[idx] = trapz(spec_Aeff[10**E_Aeff>=emin], E_Aeff[10**E_Aeff>=emin])
		total_X[id] = out
		if not z0%2==0:
			total_X_Aeff[id] = out_Aeff

	mceq_run.set_obs_particles(["mu+", "mu-", "pi_mu+", "pi_mu-", "k_mu+", "k_mu-", "pr_mu+", "pr_mu-"])
	mceq_run._calculate_integration_path(int_grid=Xvec, grid_var='X',force=True)
	mceq_run.solve(int_grid=Xvec, grid_var="X")
	
	for id in channels2:
		out = np.empty_like(Xvec)
		out_Aeff = np.empty_like(Xvec)
		for idx, X in enumerate(Xvec):
			y = (mceq_run.get_solution(id, 0, grid_idx=idx))
			out[idx] = trapz(y[e_grid>=emin], e_grid[e_grid>=emin])
			if not z0%2==0:
				spec_spline = interpolate.UnivariateSpline(e_grid[e_grid>=emin], y[e_grid>=emin], s=0)
				spec_Aeff = spec_spline(10**E_Aeff) * 100.**2 * Aeff * 2.*np.pi 
				out_Aeff[idx] = trapz(spec_Aeff[10**E_Aeff>=emin], E_Aeff[10**E_Aeff>=emin])
		total_X[id] = out
		if not z0%2==0:
			total_X_Aeff[id] = out_Aeff
	

	### write outfile ###

	out = {}
	out['E_grid'] = e_grid
	out['X_grid'] = Xvec
	out['H_grid'] = Hvec
	out['spectrum'] = total
	out['spectrum_Aeff'] = total_Aeff
	out['depth'] = total_X
	out['depth_Aeff'] = total_X_Aeff
	out['flux'] = total_sum
	out['flux_Aeff'] = total_sum_Aeff

	outfile = open(outfile, 'w')
	pickle.dump(out, outfile)
	outfile.close()
	
	z0+=1
