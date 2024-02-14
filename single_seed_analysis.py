import tfs
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde


LIST_PARAMETERS = ['DBETX_MIN','DBETY_MIN','DBETX_RMS','DBETY_RMS','DBETX_MAX','DBETY_MAX']
NBINS = 100
SEED = 0

data = tfs.read(f"permutation.tfs")
print(f"{np.unique(data.SEED)=}")

APPENDIX = "_AllSeed"
if SEED >= 0:
    data=data[data.SEED == SEED]
    APPENDIX = f"_Seed{SEED}"
print(f"{len(data)=}")

unique, counts = np.unique(data['PERM_Q1'], return_counts=True) # dict(zip(unique, counts))
select_value = unique[counts == max(counts)][0]
mask_Q2 = (data['PERM_Q1'] == select_value)

unique, counts = np.unique(data['PERM_Q2'], return_counts=True) # dict(zip(unique, counts))
select_value = unique[counts == max(counts)][0]
mask_Q1 = (data['PERM_Q2'] == select_value)


if True:
    for kk1 in ['MIN', 'RMS', 'MAX']:
        fig, ax = plt.subplots(1,2,figsize=(14,7),sharey=True)
        print(f' * Hist for DBETX and DBETY {kk1}')
 
        # Make the plot
        norm2 = np.sqrt(data.loc[mask_Q1,f'DBETX_{kk1}']**2 + data.loc[mask_Q1,f'DBETY_{kk1}']**2)
        ax[0].hist(data.loc[mask_Q1,f'DBETX_{kk1}'], bins=50, density=True, label=f'DBETX_{kk1}', alpha=0.6)
        ax[0].hist(data.loc[mask_Q1,f'DBETY_{kk1}'], bins=50, density=True, label=f'DBETY_{kk1}', alpha=0.6)
        ax[0].hist(norm2, bins=50, density=True, label=r'$\sqrt{both}$', alpha=0.6)
        
        norm2 = np.sqrt(data.loc[mask_Q2,f'DBETX_{kk1}']**2 + data.loc[mask_Q2,f'DBETY_{kk1}']**2)
        ax[1].hist(data.loc[mask_Q2,f'DBETX_{kk1}'], bins=50, density=True, label=f'DBETX_{kk1}', alpha=0.6)
        ax[1].hist(data.loc[mask_Q2,f'DBETY_{kk1}'], bins=50, density=True, label=f'DBETY_{kk1}', alpha=0.6)
        ax[1].hist(norm2, bins=50, density=True, label=r'$\sqrt{both}$', alpha=0.6)
        
        ax[0].set_xlabel(f'Beta-Beating {kk1}')
        ax[0].set_ylabel(f'Density')
        ax[0].set_title(f'Scan for Q1-Q3')
        
        ax[1].set_xlabel(f'Beta-Beating {kk1}')
        ax[1].set_title(f'Scan for Q2')
        
        ax[0].legend()
        ax[1].legend()
        
        
        fig.savefig(f"plot/hist_DBET_{kk1}{APPENDIX}.pdf",bbox_inches='tight')
        
        
        fig, ax = plt.subplots(1,2,figsize=(14,7),sharey=True)
 
        # Make the plot
        norm2 = np.sqrt(data.loc[mask_Q1,f'DBETX_{kk1}']**2 + data.loc[mask_Q1,f'DBETY_{kk1}']**2)
        ax[0].hist(data.loc[mask_Q1,f'DBETX_{kk1}'], bins=50, density=True, cumulative=True, label=f'DBETX_{kk1}', alpha=0.6)
        ax[0].hist(data.loc[mask_Q1,f'DBETY_{kk1}'], bins=50, density=True, cumulative=True, label=f'DBETY_{kk1}', alpha=0.6)
        ax[0].hist(norm2, bins=50, density=True, cumulative=True, label=r'$\sqrt{both}$', alpha=0.6)
        
        norm2 = np.sqrt(data.loc[mask_Q2,f'DBETX_{kk1}']**2 + data.loc[mask_Q2,f'DBETY_{kk1}']**2)
        ax[1].hist(data.loc[mask_Q2,f'DBETX_{kk1}'], bins=50, density=True, cumulative=True, label=f'DBETX_{kk1}', alpha=0.6)
        ax[1].hist(data.loc[mask_Q2,f'DBETY_{kk1}'], bins=50, density=True, cumulative=True, label=f'DBETY_{kk1}', alpha=0.6)
        ax[1].hist(norm2, bins=50, density=True, cumulative=True, label=r'$\sqrt{both}$', alpha=0.6)
        
        ax[0].set_xlabel(f'Beta-Beating {kk1}')
        ax[0].set_ylabel(f'Cumulative Density')
        ax[0].set_title(f'Scan for Q1-Q3')
        
        ax[1].set_xlabel(f'Beta-Beating {kk1}')
        ax[1].set_title(f'Scan for Q2')
        
        ax[0].legend()
        ax[1].legend()
        
        
        fig.savefig(f"plot/hist_cumul_DBET_{kk1}{APPENDIX}.pdf",bbox_inches='tight')
    plt.close()
        
	
    for ii1,kk1 in enumerate(LIST_PARAMETERS):
        for kk2 in LIST_PARAMETERS[ii1+1:]:
            #plt.scatter(data[kk2],data[kk1])
        
            fig, ax = plt.subplots(1,2,figsize=(14,7),sharey=True)
            print(f' * Correlation between {kk2} - {kk1}')
            
            # select Data (for Q1-Q3 scan)
            x = data.loc[mask_Q1,kk2]
            y = data.loc[mask_Q1,kk1]
    
            # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
            #k = gaussian_kde([x,y])
            #xi, yi = np.mgrid[x.min():x.max():NBINS*1j, y.min():y.max():NBINS*1j]
            #zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
            # Make the plot
            #im = ax[0].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            #fig.colorbar(im, ax=ax)
            ax[0].hist2d(x, y, bins=(NBINS, NBINS), cmap=plt.cm.seismic)
        
            ax[0].set_xlabel(kk2)
            ax[0].set_ylabel(kk1)
            ax[0].set_title(f'Scan for Q1-Q3 (r = {np.corrcoef(x, y)[0,1]:6E})')
            
            
            
            # select Data (for Q2 scan)
            x = data.loc[mask_Q2,kk2]
            y = data.loc[mask_Q2,kk1]
    
            # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
            #k = gaussian_kde([x,y])
            #xi, yi = np.mgrid[x.min():x.max():NBINS*1j, y.min():y.max():NBINS*1j]
            #zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
            # Make the plot
            #im = ax[0].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            #fig.colorbar(im, ax=ax)
            ax[1].hist2d(x, y, bins=(NBINS, NBINS), cmap=plt.cm.seismic)
        
            ax[1].set_xlabel(kk2)
            #ax[1].set_ylabel(kk1)
            ax[1].set_title(f'Scan for Q2 (r = {np.corrcoef(x, y)[0,1]:6E})')
            
            
        
            #fig = plt.gcf()
            #fig.set_size_inches(8, 8)
            fig.savefig(f"plot/scatter_{kk2}-{kk1}{APPENDIX}.pdf",bbox_inches='tight')
        
