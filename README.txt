Author: Andrew Bergman 
Project: Optimizing mGPS 15 CTS 

GITHUB: https://github.com/AndrewBergman1/optimizing_mGPS/tree/main/mGPS_proj

****PLEASE OBSERVE THAT ALL DEPENDENCIES STATED BELOW ARE COMPATIBLE WITH R 4.3.2!****


-------------------------------------------------------------------------------------------------------------------
mGPS 
-------------------------------------------------------------------------------------------------------------------
Any script you find in soil_scripts/mGPS/ and subway_scripts/mGPS are mGPS scripts


Description: 
The scripts run different variants of mGPS.  The soil_scripts load the soil data, preprocesses it, performs feature selection with either RFE or BORUTA, then runs mGPS (1 fold! Not CV). 


The map plots are generated using leaflet. Leaflet is directly compatible with the code. 

Dependencies:
Leaflet 
Boruta
caret
xgboost
doParallel
tibble
smotefamily
geosphere
rnaturalearth
rnaturalearthdata
sf
sp
imbalance

-------------------------------------------------------------------------------------------------------------------
SMOTE 
-------------------------------------------------------------------------------------------------------------------

soil_scripts/SMOTE/smote.R 

Description: 
The script loads a dataset then generates synthetic samples and add it to the dataset. Different SMOTE settings can be tested, PCAs comparing the original and smoted datasets can be generated, pie charts illustrating the sample balance in your target variable can be generated.

Dependencies:

leaflet
UBL
caret
FactoMinR
ggfortify
facctoextra
rgeodata
cluster
sf
sp
rnaturalearth
rnaturalearthdata


-------------------------------------------------------------------------------------------------------------------
Interpolation
-------------------------------------------------------------------------------------------------------------------

soil_scripts/INTERPOLATION/interpolation.R

Description:
The script tests different interpolation methods and utilizes the best methods for each feature to interpolate synthetetic samples to the dataset.

Dependencies:
gstat
sp
caret
leaflet
sf
mgcb
MBA
nngeo
RANN
fields

-------------------------------------------------------------------------------------------------------------------
Aggregating species by family prefix
-------------------------------------------------------------------------------------------------------------------

subway_scripts/aggregate_by_prefix.R

Description:
The script identifies the prefixes in the features, then aggregates the features by the prefixes.

Dependencies:
None to declare

-------------------------------------------------------------------------------------------------------------------
Clustering data
-------------------------------------------------------------------------------------------------------------------
soil_scripts/K_MEANS_CLUSTERING/cluster_data.R

Description:
The script identifies the optimal number of clusters for each country, clusters the samples accordingly.

Dependencies:
Cluster
ClusterR

