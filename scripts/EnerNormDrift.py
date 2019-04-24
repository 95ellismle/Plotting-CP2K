allStats = []
for p in all_p:
    allStats.append({
       'ener_cons': np.mean(p.ener_drift_per_rep['Total']),
       'norm': p.norm_drift,
       'scaling': p.run_inp_params['SCALING_FACTOR'],
       'NS': p.run_inp_params['NUCLEAR_TIMESTEP'],
       'EperN': p.run_inp_params['ELECTRONIC_PARTIAL_STEP'],
       'ES': p.run_inp_params['NUCLEAR_TIMESTEP']/float(p.run_inp_params['ELECTRONIC_PARTIAL_STEP'])
                     })


allES = []
allNS = []
allEnerCons = []
allNorm = []
allEperN = []
for j, i in enumerate(allStats):
    allES.append(i['ES'])
    allNS.append(i['NS'])
    allEperN.append(i['EperN'])
    allEnerCons.append(i['ener_cons'])
    allNorm.append(i['norm'])
allES = np.array(allES)
allNS = np.array(allNS)
allEnerCons = np.array(allEnerCons)
allNorm = np.array(allNorm)
allEperN = np.array(allEperN)

mask = allNS == 0.1
xdata = allES[mask]
ydata = allNorm[mask]

f, ax = plt.subplots()
ax.plot(xdata, ydata, 'ko')

fit = np.polyfit(xdata, ydata, 2)
xfit = np.arange(0, 0.101, 0.001)
ax.plot(xfit, np.polyval(fit, xfit))

ax.set_ylabel(r"Norm Drift [ps$^{-1}$]")
ax.set_xlabel(r"Electronic Timestep [fs]")

ax.set_title("Nuclear Timestep = 0.1fs")

plt.show()
