import geopandas as gpd
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt

mat = scipy.io.loadmat('../input/yaoresults/p_est.mat')

beats = gpd.read_file('../input/apd-beats-zones-neighborhoods/APD_BEAT_Jan2020_region/APD_BEAT_Jan2020_region.shp')
# p_mat = pd.read_csv('../input/yaoresults/p_est.mat')
p_mat = pd.DataFrame(mat['p_est'],columns=['p'])

items = list(beats.BEAT.values)
beat_list = list(dict.fromkeys(items))
beat_list.remove(50)
beat_list.sort()
p_mat['BEAT']=beat_list
beats = pd.merge(beats,p_mat,on='BEAT')

fig, ax = plt.subplots(figsize = (10,10))
beats.plot(ax=ax,cmap='viridis',column='p',legend=True)
ax.set_title('Discovery Rate By Beat')
plt.show()
plt.savefig('p_values.png')
