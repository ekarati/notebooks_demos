---
title: "The Boston Light Swim temperature analysis with Python"
layout: notebook

---

In the past we demonstrated how to perform a CSW catalog search with [`OWSLib`](https://ioos.github.io/notebooks_demos//notebooks/2016-12-19-exploring_csw),
and how to obtain near real-time data with [`pyoos`](https://ioos.github.io/notebooks_demos//notebooks/2016-10-12-fetching_data).
In this notebook we will use both to find all observations and model data around the Boston Harbor to access the sea water temperature.


This workflow is part of an example to advise swimmers of the annual [Boston lighthouse swim](http://bostonlightswim.org/) of the Boston Harbor water temperature conditions prior to the race. For more information regarding the workflow presented here see [Signell, Richard P.; Fernandes, Filipe; Wilcox, Kyle.   2016. "Dynamic Reusable Workflows for Ocean Science." *J. Mar. Sci. Eng.* 4, no. 4: 68](http://dx.doi.org/10.3390/jmse4040068).

<div class="prompt input_prompt">
In&nbsp;[1]:
</div>

```python
import warnings

# Suppresing warnings for a "pretty output."
warnings.simplefilter('ignore')
```

This notebook is quite big and complex,
so to help us keep things organized we'll define a cell with the most important options and switches.

Below we can define the date,
bounding box, phenomena `SOS` and `CF` names and units,
and the catalogs we will search.

<div class="prompt input_prompt">
In&nbsp;[2]:
</div>

```python
%%writefile config.yaml

# Specify a YYYY-MM-DD hh:mm:ss date or integer day offset.
# If both start and stop are offsets they will be computed relative to datetime.today() at midnight.
# Use the dates commented below to reproduce the last Boston Light Swim event forecast.
date:
    start: -5 # 2016-8-16 00:00:00
    stop: +4 # 2016-8-29 00:00:00

run_name: 'latest'

# Boston harbor.
region:
    bbox: [-71.3, 42.03, -70.57, 42.63]
    # Try the bounding box below to see how the notebook will behave for a different region.
    #bbox: [-74.5, 40, -72., 41.5]
    crs: 'urn:ogc:def:crs:OGC:1.3:CRS84'

sos_name: 'sea_water_temperature'

cf_names:
    - sea_water_temperature
    - sea_surface_temperature
    - sea_water_potential_temperature
    - equivalent_potential_temperature
    - sea_water_conservative_temperature
    - pseudo_equivalent_potential_temperature

units: 'celsius'

catalogs:
    - https://data.ioos.us/csw
```
<div class="output_area"><div class="prompt"></div>
<pre>
    Overwriting config.yaml

</pre>
</div>
We'll print some of the search configuration options along the way to keep track of them.

<div class="prompt input_prompt">
In&nbsp;[3]:
</div>

```python
import os
import shutil
from datetime import datetime
from ioos_tools.ioos import parse_config

config = parse_config('config.yaml')

# Saves downloaded data into a temporary directory.
save_dir = os.path.abspath(config['run_name'])
if os.path.exists(save_dir):
    shutil.rmtree(save_dir)
os.makedirs(save_dir)

fmt = '{:*^64}'.format
print(fmt('Saving data inside directory {}'.format(save_dir)))
print(fmt(' Run information '))
print('Run date: {:%Y-%m-%d %H:%M:%S}'.format(datetime.utcnow()))
print('Start: {:%Y-%m-%d %H:%M:%S}'.format(config['date']['start']))
print('Stop: {:%Y-%m-%d %H:%M:%S}'.format(config['date']['stop']))
print('Bounding box: {0:3.2f}, {1:3.2f},'
      '{2:3.2f}, {3:3.2f}'.format(*config['region']['bbox']))
```
<div class="output_area"><div class="prompt"></div>
<pre>
    Saving data inside directory /home/filipe/IOOS/notebooks_demos/notebooks/latest
    *********************** Run information ************************
    Run date: 2019-02-06 19:46:20
    Start: 2019-02-01 00:00:00
    Stop: 2019-02-10 00:00:00
    Bounding box: -71.30, 42.03,-70.57, 42.63

</pre>
</div>
We already created an `OWSLib.fes` filter [before](https://ioos.github.io/notebooks_demos//notebooks/2016-12-19-exploring_csw).
The main difference here is that we do not want the atmosphere model data,
so we are filtering out all the `GRIB-2` data format.

<div class="prompt input_prompt">
In&nbsp;[4]:
</div>

```python
def make_filter(config):
    from owslib import fes
    from ioos_tools.ioos import fes_date_filter
    kw = dict(wildCard='*', escapeChar='\\',
              singleChar='?', propertyname='apiso:AnyText')

    or_filt = fes.Or([fes.PropertyIsLike(literal=('*%s*' % val), **kw)
                      for val in config['cf_names']])

    not_filt = fes.Not([fes.PropertyIsLike(literal='GRIB-2', **kw)])

    begin, end = fes_date_filter(config['date']['start'],
                                 config['date']['stop'])
    bbox_crs = fes.BBox(config['region']['bbox'],
                        crs=config['region']['crs'])
    filter_list = [fes.And([bbox_crs, begin, end, or_filt, not_filt])]
    return filter_list


filter_list = make_filter(config)
```

In the cell below we ask the catalog for all the returns that match the filter and have an OPeNDAP endpoint.

<div class="prompt input_prompt">
In&nbsp;[5]:
</div>

```python
from ioos_tools.ioos import service_urls, get_csw_records
from owslib.csw import CatalogueServiceWeb


dap_urls = []
print(fmt(' Catalog information '))
for endpoint in config['catalogs']:
    print('URL: {}'.format(endpoint))
    try:
        csw = CatalogueServiceWeb(endpoint, timeout=120)
    except Exception as e:
        print('{}'.format(e))
        continue
    csw = get_csw_records(csw, filter_list, esn='full')
    OPeNDAP = service_urls(csw.records, identifier='OPeNDAP:OPeNDAP')
    odp = service_urls(csw.records, identifier='urn:x-esri:specification:ServiceType:odp:url')
    dap = OPeNDAP + odp
    dap_urls.extend(dap)

    print('Number of datasets available: {}'.format(len(csw.records.keys())))

    for rec, item in csw.records.items():
        print('{}'.format(item.title))
    if dap:
        print(fmt(' DAP '))
        for url in dap:
            print('{}.html'.format(url))
    print('\n')

# Get only unique endpoints.
dap_urls = list(set(dap_urls))
```
<div class="output_area"><div class="prompt"></div>
<pre>
    ********************* Catalog information **********************
    URL: https://data.ioos.us/csw
    Number of datasets available: 33
    ROMS doppio Real-Time Operational PSAS Forecast System Version 1 FMRC Averages
    ROMS doppio Real-Time Operational PSAS Forecast System Version 1 FMRC History
    urn:ioos:station:NOAA.NOS.CO-OPS:8443970 station, Boston, MA
    A01 Accelerometer - Waves
    A01 Directional Waves (waves.mstrain Experimental)
    A01 Met - Meteorology
    A01 Optics - Chlorophyll / Turbidity
    A01 Optode - Oxygen
    A01 SBE16 Oxygen
    A01 Sbe37 - CTD
    BOSTON 16 NM East of Boston, MA
    Buoy A01 - Massachusetts Bay
    COAWST Modeling System: USEast: ROMS-WRF-SWAN coupled model (aka CNAPS)
    Coupled Northwest Atlantic Prediction System (CNAPS)
    Department of Physical Oceanography, School of Marine Sciences, University of Maine A01 Accelerometer Buoy Sensor
    Department of Physical Oceanography, School of Marine Sciences, University of Maine A01 Met Buoy Sensor
    Department of Physical Oceanography, School of Marine Sciences, University of Maine A01 Optode 51m Buoy Sensor
    Department of Physical Oceanography, School of Marine Sciences, University of Maine A01 Sbe37 1m Buoy Sensor
    Department of Physical Oceanography, School of Marine Sciences, University of Maine A01 Sbe37 20m Buoy Sensor
    Directional wave and sea surface temperature measurements collected in situ by Datawell DWR-M3 directional buoy located near CLATSOP SPIT, OR from 2018/02/07 22:03:45 to 2019/02/05 18:00:25.
    Directional wave and sea surface temperature measurements collected in situ by Datawell DWR-M3 directional buoy located near GRAYS HARBOR, WA from 2018/06/28 19:03:45 to 2019/02/05 18:00:25.
    G1SST, 1km blended SST
    Global SST & Sea Ice Analysis, L4 OSTIA, 0.05 deg daily (METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2)
    NDBC Standard Meteorological Buoy Data, 1970-present
    NECOFS (FVCOM) - Scituate - Latest Forecast
    NECOFS Massachusetts (FVCOM) - Boston - Latest Forecast
    NECOFS Massachusetts (FVCOM) - Massachusetts Coastal - Latest Forecast
    NERACOOS Gulf of Maine Ocean Array: Realtime Buoy Observations: A01 Massachusetts Bay: A01 ACCELEROMETER Massachusetts Bay
    NERACOOS Gulf of Maine Ocean Array: Realtime Buoy Observations: A01 Massachusetts Bay: A01 CTD1m Massachusetts Bay
    NERACOOS Gulf of Maine Ocean Array: Realtime Buoy Observations: A01 Massachusetts Bay: A01 CTD20m Massachusetts Bay
    NERACOOS Gulf of Maine Ocean Array: Realtime Buoy Observations: A01 Massachusetts Bay: A01 MET Massachusetts Bay
    NERACOOS Gulf of Maine Ocean Array: Realtime Buoy Observations: A01 Massachusetts Bay: A01 OPTICS3m Massachusetts Bay
    NERACOOS Gulf of Maine Ocean Array: Realtime Buoy Observations: A01 Massachusetts Bay: A01 OPTODE51m Massachusetts Bay
    ***************************** DAP ******************************
    http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/avg/Averages_Best.html
    http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best.html
    http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/036p1_rt.nc.html
    http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/162p1_rt.nc.html
    http://thredds.secoora.org/thredds/dodsC/AOOS_OSTIA.nc.html
    http://thredds.secoora.org/thredds/dodsC/G1_SST_GLOBAL.nc.html
    http://thredds.secoora.org/thredds/dodsC/SECOORA_NCSU_CNAPS.nc.html
    http://www.neracoos.org/thredds/dodsC/UMO/DSG/SOS/A01/Accelerometer/HistoricRealtime.html
    http://www.neracoos.org/thredds/dodsC/UMO/DSG/SOS/A01/CTD1m/HistoricRealtime.html
    http://www.neracoos.org/thredds/dodsC/UMO/DSG/SOS/A01/CTD20m/HistoricRealtime.html
    http://www.neracoos.org/thredds/dodsC/UMO/DSG/SOS/A01/Met/HistoricRealtime.html
    http://www.neracoos.org/thredds/dodsC/UMO/DSG/SOS/A01/OPTICS_S3m/HistoricRealtime.html
    http://www.neracoos.org/thredds/dodsC/UMO/DSG/SOS/A01/OPTODE51m/HistoricRealtime.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.accelerometer.realtime.nc.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.met.realtime.nc.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.optode.realtime.51m.nc.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.sbe37.realtime.1m.nc.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.sbe37.realtime.20m.nc.html
    http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_BOSTON_FORECAST.nc.html
    http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc.html
    http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_SCITUATE_FORECAST.nc.html
    
    

</pre>
</div>
We found some models, and observations from NERACOOS there.
However, we do know that there are some buoys from NDBC and CO-OPS available too.
Also, those NERACOOS observations seem to be from a [CTD](http://www.neracoos.org/thredds/dodsC/UMO/DSG/SOS/A01/CTD1m/HistoricRealtime/Agg.ncml.html) mounted at 65 meters below the sea surface. Rendering them useless from our purpose.

So let's use the catalog only for the models by filtering the observations with `is_station` below.
And we'll rely `CO-OPS` and `NDBC` services for the observations.

<div class="prompt input_prompt">
In&nbsp;[6]:
</div>

```python
from timeout_decorator import TimeoutError
from ioos_tools.ioos import is_station

# Filter out some station endpoints.
non_stations = []
for url in dap_urls:
    url = f'{url}#fillmismatch'
    try:
        if not is_station(url):
            non_stations.append(url)
    except (IOError, OSError, RuntimeError, TimeoutError) as e:
        print('Could not access URL {}.html\n{!r}'.format(url, e))

dap_urls = non_stations

print(fmt(' Filtered DAP '))
for url in dap_urls:
    print('{}.html'.format(url))
```
<div class="output_area"><div class="prompt"></div>
<pre>
    ************************* Filtered DAP *************************
    http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc#fillmismatch.html
    http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_BOSTON_FORECAST.nc#fillmismatch.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.optode.realtime.51m.nc#fillmismatch.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.accelerometer.realtime.nc#fillmismatch.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.met.realtime.nc#fillmismatch.html
    http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/avg/Averages_Best#fillmismatch.html
    http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best#fillmismatch.html
    http://thredds.secoora.org/thredds/dodsC/AOOS_OSTIA.nc#fillmismatch.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.sbe37.realtime.1m.nc#fillmismatch.html
    http://thredds.secoora.org/thredds/dodsC/G1_SST_GLOBAL.nc#fillmismatch.html
    http://thredds.secoora.org/thredds/dodsC/SECOORA_NCSU_CNAPS.nc#fillmismatch.html
    http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_SCITUATE_FORECAST.nc#fillmismatch.html
    http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.sbe37.realtime.20m.nc#fillmismatch.html

</pre>
</div>
Now we can use `pyoos` collectors for `NdbcSos`,

<div class="prompt input_prompt">
In&nbsp;[7]:
</div>

```python
from pyoos.collectors.ndbc.ndbc_sos import NdbcSos

collector_ndbc = NdbcSos()

collector_ndbc.set_bbox(config['region']['bbox'])
collector_ndbc.end_time = config['date']['stop']
collector_ndbc.start_time = config['date']['start']
collector_ndbc.variables = [config['sos_name']]

ofrs = collector_ndbc.server.offerings
title = collector_ndbc.server.identification.title
print(fmt(' NDBC Collector offerings '))
print('{}: {} offerings'.format(title, len(ofrs)))
```
<div class="output_area"><div class="prompt"></div>
<pre>
    ******************* NDBC Collector offerings *******************
    National Data Buoy Center SOS: 1037 offerings

</pre>
</div>
<div class="prompt input_prompt">
In&nbsp;[8]:
</div>

```python
import pandas as pd
from ioos_tools.ioos import collector2table

ndbc = collector2table(collector=collector_ndbc,
                       config=config,
                       col='sea_water_temperature (C)')

if ndbc:
    data = dict(
        station_name=[s._metadata.get('station_name') for s in ndbc],
        station_code=[s._metadata.get('station_code') for s in ndbc],
        sensor=[s._metadata.get('sensor') for s in ndbc],
        lon=[s._metadata.get('lon') for s in ndbc],
        lat=[s._metadata.get('lat') for s in ndbc],
        depth=[s._metadata.get('depth') for s in ndbc],
    )

table = pd.DataFrame(data).set_index('station_code')
table
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>station_name</th>
      <th>sensor</th>
      <th>lon</th>
      <th>lat</th>
      <th>depth</th>
    </tr>
    <tr>
      <th>station_code</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>44013</th>
      <td>BOSTON 16 NM East of Boston, MA</td>
      <td>urn:ioos:sensor:wmo:44013::watertemp1</td>
      <td>-70.651</td>
      <td>42.346</td>
      <td>0.6</td>
    </tr>
  </tbody>
</table>
</div>



and `CoopsSos`.

<div class="prompt input_prompt">
In&nbsp;[9]:
</div>

```python
from pyoos.collectors.coops.coops_sos import CoopsSos

collector_coops = CoopsSos()

collector_coops.set_bbox(config['region']['bbox'])
collector_coops.end_time = config['date']['stop']
collector_coops.start_time = config['date']['start']
collector_coops.variables = [config['sos_name']]

ofrs = collector_coops.server.offerings
title = collector_coops.server.identification.title
print(fmt(' Collector offerings '))
print('{}: {} offerings'.format(title, len(ofrs)))
```
<div class="output_area"><div class="prompt"></div>
<pre>
    ********************* Collector offerings **********************
    NOAA.NOS.CO-OPS SOS: 1230 offerings

</pre>
</div>
<div class="prompt input_prompt">
In&nbsp;[10]:
</div>

```python
coops = collector2table(collector=collector_coops,
                        config=config,
                        col='sea_water_temperature (C)')

if coops:
    data = dict(
        station_name=[s._metadata.get('station_name') for s in coops],
        station_code=[s._metadata.get('station_code') for s in coops],
        sensor=[s._metadata.get('sensor') for s in coops],
        lon=[s._metadata.get('lon') for s in coops],
        lat=[s._metadata.get('lat') for s in coops],
        depth=[s._metadata.get('depth') for s in coops],
    )

table = pd.DataFrame(data).set_index('station_code')
table
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>station_name</th>
      <th>sensor</th>
      <th>lon</th>
      <th>lat</th>
      <th>depth</th>
    </tr>
    <tr>
      <th>station_code</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>44013</th>
      <td>BOSTON 16 NM East of Boston, MA</td>
      <td>urn:ioos:sensor:wmo:44013::watertemp1</td>
      <td>-70.651</td>
      <td>42.346</td>
      <td>0.6</td>
    </tr>
  </tbody>
</table>
</div>



We will join all the observations into an uniform series, interpolated to 1-hour interval, for the model-data comparison.

This step is necessary because the observations can be 7 or 10 minutes resolution,
while the models can be 30 to 60 minutes.

<div class="prompt input_prompt">
In&nbsp;[11]:
</div>

```python
data = ndbc + coops

index = pd.date_range(start=config['date']['start'].replace(tzinfo=None),
                      end=config['date']['stop'].replace(tzinfo=None),
                      freq='1H')

# Preserve metadata with `reindex`.
observations = []
for series in data:
    _metadata = series._metadata
    series.index = series.index.tz_localize(None)
    series.index = series.index.tz_localize(None)
    obs = series.reindex(index=index, limit=1, method='nearest')
    obs._metadata = _metadata
    observations.append(obs)
```

In this next cell we will save the data for quicker access later.

<div class="prompt input_prompt">
In&nbsp;[12]:
</div>

```python
import iris
from ioos_tools.tardis import series2cube

attr = dict(
    featureType='timeSeries',
    Conventions='CF-1.6',
    standard_name_vocabulary='CF-1.6',
    cdm_data_type='Station',
    comment='Data from http://opendap.co-ops.nos.noaa.gov'
)


cubes = iris.cube.CubeList(
    [series2cube(obs, attr=attr) for obs in observations]
)

outfile = os.path.join(save_dir, 'OBS_DATA.nc')
iris.save(cubes, outfile)
```

Taking a quick look at the observations:

<div class="prompt input_prompt">
In&nbsp;[13]:
</div>

```python
%matplotlib inline

ax = pd.concat(data).plot(figsize=(11, 2.25))
```


![png](data:image/png;base64,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
)


Now it is time to loop the models we found above,

<div class="prompt input_prompt">
In&nbsp;[14]:
</div>

```python
from iris.exceptions import (CoordinateNotFoundError, ConstraintMismatchError,
                             MergeError)
from ioos_tools.ioos import get_model_name
from ioos_tools.tardis import quick_load_cubes, proc_cube, is_model, get_surface

print(fmt(' Models '))
cubes = dict()
for k, url in enumerate(dap_urls):
    print('\n[Reading url {}/{}]: {}'.format(k+1, len(dap_urls), url))
    try:
        cube = quick_load_cubes(url, config['cf_names'],
                                callback=None, strict=True)
        if is_model(cube):
            cube = proc_cube(cube,
                             bbox=config['region']['bbox'],
                             time=(config['date']['start'],
                                   config['date']['stop']),
                             units=config['units'])
        else:
            print('[Not model data]: {}'.format(url))
            continue
        cube = get_surface(cube)
        mod_name = get_model_name(url)
        cubes.update({mod_name: cube})
    except (RuntimeError, ValueError,
            ConstraintMismatchError, CoordinateNotFoundError,
            IndexError) as e:
        print('Cannot get cube for: {}\n{}'.format(url, e))
```
<div class="output_area"><div class="prompt"></div>
<pre>
    **************************** Models ****************************
    
    [Reading url 1/13]: http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc#fillmismatch
    
    [Reading url 2/13]: http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_BOSTON_FORECAST.nc#fillmismatch
    
    [Reading url 3/13]: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.optode.realtime.51m.nc#fillmismatch
    [Not model data]: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.optode.realtime.51m.nc#fillmismatch
    
    [Reading url 4/13]: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.accelerometer.realtime.nc#fillmismatch
    Cannot get cube for: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.accelerometer.realtime.nc#fillmismatch
    Cannot find ['sea_water_temperature', 'sea_surface_temperature', 'sea_water_potential_temperature', 'equivalent_potential_temperature', 'sea_water_conservative_temperature', 'pseudo_equivalent_potential_temperature'] in http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.accelerometer.realtime.nc#fillmismatch.
    
    [Reading url 5/13]: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.met.realtime.nc#fillmismatch
    Cannot get cube for: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.met.realtime.nc#fillmismatch
    Cannot find ['sea_water_temperature', 'sea_surface_temperature', 'sea_water_potential_temperature', 'equivalent_potential_temperature', 'sea_water_conservative_temperature', 'pseudo_equivalent_potential_temperature'] in http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.met.realtime.nc#fillmismatch.
    
    [Reading url 6/13]: http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/avg/Averages_Best#fillmismatch
    
    [Reading url 7/13]: http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best#fillmismatch
    
    [Reading url 8/13]: http://thredds.secoora.org/thredds/dodsC/AOOS_OSTIA.nc#fillmismatch
    [Not model data]: http://thredds.secoora.org/thredds/dodsC/AOOS_OSTIA.nc#fillmismatch
    
    [Reading url 9/13]: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.sbe37.realtime.1m.nc#fillmismatch
    [Not model data]: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.sbe37.realtime.1m.nc#fillmismatch
    
    [Reading url 10/13]: http://thredds.secoora.org/thredds/dodsC/G1_SST_GLOBAL.nc#fillmismatch
    Cannot get cube for: http://thredds.secoora.org/thredds/dodsC/G1_SST_GLOBAL.nc#fillmismatch
    Found more than one time coordinates!
    
    [Reading url 11/13]: http://thredds.secoora.org/thredds/dodsC/SECOORA_NCSU_CNAPS.nc#fillmismatch
    
    [Reading url 12/13]: http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_SCITUATE_FORECAST.nc#fillmismatch
    
    [Reading url 13/13]: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.sbe37.realtime.20m.nc#fillmismatch
    [Not model data]: http://www.neracoos.org/thredds/dodsC/UMO/Realtime/SOS/A01/DSG_A0140.sbe37.realtime.20m.nc#fillmismatch

</pre>
</div>
Next, we will match them with the nearest observed time-series. The `max_dist=0.08` is in degrees, that is roughly 8 kilometers.

<div class="prompt input_prompt">
In&nbsp;[15]:
</div>

```python
import iris
from iris.pandas import as_series
from ioos_tools.tardis import (make_tree, get_nearest_water,
                               add_station, ensure_timeseries, remove_ssh)

for mod_name, cube in cubes.items():
    fname = '{}.nc'.format(mod_name)
    fname = os.path.join(save_dir, fname)
    print(fmt(' Downloading to file {} '.format(fname)))
    try:
        tree, lon, lat = make_tree(cube)
    except CoordinateNotFoundError as e:
        print('Cannot make KDTree for: {}'.format(mod_name))
        continue
    # Get model series at observed locations.
    raw_series = dict()
    for obs in observations:
        obs = obs._metadata
        station = obs['station_code']
        try:
            kw = dict(k=10, max_dist=0.08, min_var=0.01)
            args = cube, tree, obs['lon'], obs['lat']
            try:
                series, dist, idx = get_nearest_water(*args, **kw)
            except RuntimeError as e:
                print('Cannot download {!r}.\n{}'.format(cube, e))
                series = None
        except ValueError as e:
            status = 'No Data'
            print('[{}] {}'.format(status, obs['station_name']))
            continue
        if not series:
            status = 'Land   '
        else:
            raw_series.update({station: series})
            series = as_series(series)
            status = 'Water  '
        print('[{}] {}'.format(status, obs['station_name']))
    if raw_series:  # Save cube.
        for station, cube in raw_series.items():
            cube = add_station(cube, station)
            cube = remove_ssh(cube)
        try:
            cube = iris.cube.CubeList(raw_series.values()).merge_cube()
        except MergeError as e:
            print(e)
        ensure_timeseries(cube)
        try:
            iris.save(cube, fname)
        except AttributeError:
            # FIXME: we should patch the bad attribute instead of removing everything.
            cube.attributes = {}
            iris.save(cube, fname)
        del cube
    print('Finished processing [{}]'.format(mod_name))
```
<div class="output_area"><div class="prompt"></div>
<pre>
     Downloading to file /home/filipe/IOOS/notebooks_demos/notebooks/latest/Forecasts-NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc#fillmismatch.nc 
    [Water  ] BOSTON 16 NM East of Boston, MA
    Finished processing [Forecasts-NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc#fillmismatch]
     Downloading to file /home/filipe/IOOS/notebooks_demos/notebooks/latest/Forecasts-NECOFS_FVCOM_OCEAN_BOSTON_FORECAST.nc#fillmismatch.nc 
    [No Data] BOSTON 16 NM East of Boston, MA
    Finished processing [Forecasts-NECOFS_FVCOM_OCEAN_BOSTON_FORECAST.nc#fillmismatch]
     Downloading to file /home/filipe/IOOS/notebooks_demos/notebooks/latest/roms_doppio_2017_da_avg-Averages_Best#fillmismatch.nc 
    [Water  ] BOSTON 16 NM East of Boston, MA
    Finished processing [roms_doppio_2017_da_avg-Averages_Best#fillmismatch]
     Downloading to file /home/filipe/IOOS/notebooks_demos/notebooks/latest/roms_doppio_2017_da-History_Best#fillmismatch.nc 
    [Water  ] BOSTON 16 NM East of Boston, MA
    Finished processing [roms_doppio_2017_da-History_Best#fillmismatch]
     Downloading to file /home/filipe/IOOS/notebooks_demos/notebooks/latest/SECOORA_NCSU_CNAPS.nc#fillmismatch.nc 
    [Water  ] BOSTON 16 NM East of Boston, MA
    Finished processing [SECOORA_NCSU_CNAPS.nc#fillmismatch]
     Downloading to file /home/filipe/IOOS/notebooks_demos/notebooks/latest/Forecasts-NECOFS_FVCOM_OCEAN_SCITUATE_FORECAST.nc#fillmismatch.nc 
    [No Data] BOSTON 16 NM East of Boston, MA
    Finished processing [Forecasts-NECOFS_FVCOM_OCEAN_SCITUATE_FORECAST.nc#fillmismatch]

</pre>
</div>
Now it is possible to compute some simple comparison metrics. First we'll calculate the model mean bias:

$$ \text{MB} = \mathbf{\overline{m}} - \mathbf{\overline{o}}$$

<div class="prompt input_prompt">
In&nbsp;[16]:
</div>

```python
from ioos_tools.ioos import stations_keys


def rename_cols(df, config):
    cols = stations_keys(config, key='station_name')
    return df.rename(columns=cols)
```

<div class="prompt input_prompt">
In&nbsp;[17]:
</div>

```python
from ioos_tools.ioos import load_ncs
from ioos_tools.skill_score import mean_bias, apply_skill

dfs = load_ncs(config)

df = apply_skill(dfs, mean_bias, remove_mean=False, filter_tides=False)
skill_score = dict(mean_bias=df.to_dict())

# Filter out stations with no valid comparison.
df.dropna(how='all', axis=1, inplace=True)
df = df.applymap('{:.2f}'.format).replace('nan', '--')
```

And the root mean squared rrror of the deviations from the mean:
$$ \text{CRMS} = \sqrt{\left(\mathbf{m'} - \mathbf{o'}\right)^2}$$

where: $\mathbf{m'} = \mathbf{m} - \mathbf{\overline{m}}$ and $\mathbf{o'} = \mathbf{o} - \mathbf{\overline{o}}$

<div class="prompt input_prompt">
In&nbsp;[18]:
</div>

```python
from ioos_tools.skill_score import rmse

dfs = load_ncs(config)

df = apply_skill(dfs, rmse, remove_mean=True, filter_tides=False)
skill_score['rmse'] = df.to_dict()

# Filter out stations with no valid comparison.
df.dropna(how='all', axis=1, inplace=True)
df = df.applymap('{:.2f}'.format).replace('nan', '--')
```

The next 2 cells make the scores "pretty" for plotting.

<div class="prompt input_prompt">
In&nbsp;[19]:
</div>

```python
import pandas as pd

# Stringfy keys.
for key in skill_score.keys():
    skill_score[key] = {str(k): v for k, v in skill_score[key].items()}

mean_bias = pd.DataFrame.from_dict(skill_score['mean_bias'])
mean_bias = mean_bias.applymap('{:.2f}'.format).replace('nan', '--')

skill_score = pd.DataFrame.from_dict(skill_score['rmse'])
skill_score = skill_score.applymap('{:.2f}'.format).replace('nan', '--')
```

<div class="prompt input_prompt">
In&nbsp;[20]:
</div>

```python
import folium
from ioos_tools.ioos import get_coordinates


def make_map(bbox, **kw):
    line = kw.pop('line', True)
    layers = kw.pop('layers', True)
    zoom_start = kw.pop('zoom_start', 5)

    lon = (bbox[0] + bbox[2]) / 2
    lat = (bbox[1] + bbox[3]) / 2
    m = folium.Map(width='100%', height='100%',
                   location=[lat, lon], zoom_start=zoom_start)

    if layers:
        url = 'http://oos.soest.hawaii.edu/thredds/wms/hioos/satellite/dhw_5km'
        w = folium.WmsTileLayer(
            url,
            name='Sea Surface Temperature',
            fmt='image/png',
            layers='CRW_SST',
            attr='PacIOOS TDS',
            overlay=True,
            transparent=True)
        w.add_to(m)

    if line:
        p = folium.PolyLine(
            get_coordinates(bbox),
            color='#FF0000',
            weight=2,
            opacity=0.9,
        )
        p.add_to(m)
    return m
```

<div class="prompt input_prompt">
In&nbsp;[21]:
</div>

```python
bbox = config['region']['bbox']

m = make_map(
    bbox,
    zoom_start=11,
    line=True,
    layers=True
)
```

The cells from `[20]` to `[25]` create a [`folium`](https://github.com/python-visualization/folium) map with [`bokeh`](http://bokeh.pydata.org/en/latest/) for the time-series at the observed points.

Note that we did mark the nearest model cell location used in the comparison.

<div class="prompt input_prompt">
In&nbsp;[22]:
</div>

```python
all_obs = stations_keys(config)

from glob import glob
from operator import itemgetter

import iris
from folium.plugins import MarkerCluster

iris.FUTURE.netcdf_promote = True

big_list = []
for fname in glob(os.path.join(save_dir, '*.nc')):
    if 'OBS_DATA' in fname:
        continue
    cube = iris.load_cube(fname)
    model = os.path.split(fname)[1].split('-')[-1].split('.')[0]
    lons = cube.coord(axis='X').points
    lats = cube.coord(axis='Y').points
    stations = cube.coord('station_code').points
    models = [model]*lons.size
    lista = zip(models, lons.tolist(), lats.tolist(), stations.tolist())
    big_list.extend(lista)

big_list.sort(key=itemgetter(3))
df = pd.DataFrame(big_list, columns=['name', 'lon', 'lat', 'station'])
df.set_index('station', drop=True, inplace=True)
groups = df.groupby(df.index)


locations, popups = [], []
for station, info in groups:
    sta_name = all_obs[station]
    for lat, lon, name in zip(info.lat, info.lon, info.name):
        locations.append([lat, lon])
        popups.append('[{}]: {}'.format(name, sta_name))

MarkerCluster(locations=locations, popups=popups, name='Cluster').add_to(m);
```

Here we use a dictionary with some models we expect to find so we can create a better legend for the plots. If any new models are found, we will use its filename in the legend as a default until we can go back and add a short name to our library.

<div class="prompt input_prompt">
In&nbsp;[23]:
</div>

```python
titles = {
    'coawst_4_use_best': 'COAWST_4',
    'global': 'HYCOM',
    'NECOFS_GOM3_FORECAST': 'NECOFS_GOM3',
    'NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST': 'NECOFS_MassBay',
    'OBS_DATA': 'Observations'
}
```

<div class="prompt input_prompt">
In&nbsp;[24]:
</div>

```python
from bokeh.resources import CDN
from bokeh.plotting import figure
from bokeh.embed import file_html
from bokeh.models import HoverTool
from itertools import cycle
from bokeh.palettes import Category20

from folium import IFrame

# Plot defaults.
colors = Category20[20]
colorcycler = cycle(colors)
tools = 'pan,box_zoom,reset'
width, height = 750, 250


def make_plot(df, station):
    p = figure(
        toolbar_location='above',
        x_axis_type='datetime',
        width=width,
        height=height,
        tools=tools,
        title=str(station)
    )
    for column, series in df.iteritems():
        series.dropna(inplace=True)
        if not series.empty:
            if 'OBS_DATA' not in column:
                bias = mean_bias[str(station)][column]
                skill = skill_score[str(station)][column]
                line_color = next(colorcycler)
                kw = dict(alpha=0.65, line_color=line_color)
            else:
                skill = bias = 'NA'
                kw = dict(alpha=1, color='crimson')
            line = p.line(
                x=series.index,
                y=series.values,
                legend='{}'.format(titles.get(column, column)),
                line_width=5,
                line_cap='round',
                line_join='round',
                **kw
            )
            p.add_tools(HoverTool(tooltips=[('Name', '{}'.format(titles.get(column, column))),
                                            ('Bias', bias),
                                            ('Skill', skill)],
                                  renderers=[line]))
    return p


def make_marker(p, station):
    lons = stations_keys(config, key='lon')
    lats = stations_keys(config, key='lat')

    lon, lat = lons[station], lats[station]
    html = file_html(p, CDN, station)
    iframe = IFrame(html, width=width+40, height=height+80)

    popup = folium.Popup(iframe, max_width=2650)
    icon = folium.Icon(color='green', icon='stats')
    marker = folium.Marker(location=[lat, lon],
                           popup=popup,
                           icon=icon)
    return marker
```

<div class="prompt input_prompt">
In&nbsp;[25]:
</div>

```python
dfs = load_ncs(config)

for station in dfs:
    sta_name = all_obs[station]
    df = dfs[station]
    if df.empty:
        continue
    p = make_plot(df, station)
    marker = make_marker(p, station)
    marker.add_to(m)

folium.LayerControl().add_to(m)

m
```




<div style="width:100%;"><div style="position:relative;width:100%;height:0;padding-bottom:60%;"><iframe src="data:text/html;charset=utf-8;base64,PCFET0NUWVBFIGh0bWw+CjxoZWFkPiAgICAKICAgIDxtZXRhIGh0dHAtZXF1aXY9ImNvbnRlbnQtdHlwZSIgY29udGVudD0idGV4dC9odG1sOyBjaGFyc2V0PVVURi04IiAvPgogICAgPHNjcmlwdD5MX1BSRUZFUl9DQU5WQVM9ZmFsc2U7IExfTk9fVE9VQ0g9ZmFsc2U7IExfRElTQUJMRV8zRD1mYWxzZTs8L3NjcmlwdD4KICAgIDxzY3JpcHQgc3JjPSJodHRwczovL2Nkbi5qc2RlbGl2ci5uZXQvbnBtL2xlYWZsZXRAMS4zLjQvZGlzdC9sZWFmbGV0LmpzIj48L3NjcmlwdD4KICAgIDxzY3JpcHQgc3JjPSJodHRwczovL2FqYXguZ29vZ2xlYXBpcy5jb20vYWpheC9saWJzL2pxdWVyeS8xLjExLjEvanF1ZXJ5Lm1pbi5qcyI+PC9zY3JpcHQ+CiAgICA8c2NyaXB0IHNyYz0iaHR0cHM6Ly9tYXhjZG4uYm9vdHN0cmFwY2RuLmNvbS9ib290c3RyYXAvMy4yLjAvanMvYm9vdHN0cmFwLm1pbi5qcyI+PC9zY3JpcHQ+CiAgICA8c2NyaXB0IHNyYz0iaHR0cHM6Ly9jZG5qcy5jbG91ZGZsYXJlLmNvbS9hamF4L2xpYnMvTGVhZmxldC5hd2Vzb21lLW1hcmtlcnMvMi4wLjIvbGVhZmxldC5hd2Vzb21lLW1hcmtlcnMuanMiPjwvc2NyaXB0PgogICAgPGxpbmsgcmVsPSJzdHlsZXNoZWV0IiBocmVmPSJodHRwczovL2Nkbi5qc2RlbGl2ci5uZXQvbnBtL2xlYWZsZXRAMS4zLjQvZGlzdC9sZWFmbGV0LmNzcyIvPgogICAgPGxpbmsgcmVsPSJzdHlsZXNoZWV0IiBocmVmPSJodHRwczovL21heGNkbi5ib290c3RyYXBjZG4uY29tL2Jvb3RzdHJhcC8zLjIuMC9jc3MvYm9vdHN0cmFwLm1pbi5jc3MiLz4KICAgIDxsaW5rIHJlbD0ic3R5bGVzaGVldCIgaHJlZj0iaHR0cHM6Ly9tYXhjZG4uYm9vdHN0cmFwY2RuLmNvbS9ib290c3RyYXAvMy4yLjAvY3NzL2Jvb3RzdHJhcC10aGVtZS5taW4uY3NzIi8+CiAgICA8bGluayByZWw9InN0eWxlc2hlZXQiIGhyZWY9Imh0dHBzOi8vbWF4Y2RuLmJvb3RzdHJhcGNkbi5jb20vZm9udC1hd2Vzb21lLzQuNi4zL2Nzcy9mb250LWF3ZXNvbWUubWluLmNzcyIvPgogICAgPGxpbmsgcmVsPSJzdHlsZXNoZWV0IiBocmVmPSJodHRwczovL2NkbmpzLmNsb3VkZmxhcmUuY29tL2FqYXgvbGlicy9MZWFmbGV0LmF3ZXNvbWUtbWFya2Vycy8yLjAuMi9sZWFmbGV0LmF3ZXNvbWUtbWFya2Vycy5jc3MiLz4KICAgIDxsaW5rIHJlbD0ic3R5bGVzaGVldCIgaHJlZj0iaHR0cHM6Ly9yYXdjZG4uZ2l0aGFjay5jb20vcHl0aG9uLXZpc3VhbGl6YXRpb24vZm9saXVtL21hc3Rlci9mb2xpdW0vdGVtcGxhdGVzL2xlYWZsZXQuYXdlc29tZS5yb3RhdGUuY3NzIi8+CiAgICA8c3R5bGU+aHRtbCwgYm9keSB7d2lkdGg6IDEwMCU7aGVpZ2h0OiAxMDAlO21hcmdpbjogMDtwYWRkaW5nOiAwO308L3N0eWxlPgogICAgPHN0eWxlPiNtYXAge3Bvc2l0aW9uOmFic29sdXRlO3RvcDowO2JvdHRvbTowO3JpZ2h0OjA7bGVmdDowO308L3N0eWxlPgogICAgCiAgICA8bWV0YSBuYW1lPSJ2aWV3cG9ydCIgY29udGVudD0id2lkdGg9ZGV2aWNlLXdpZHRoLAogICAgICAgIGluaXRpYWwtc2NhbGU9MS4wLCBtYXhpbXVtLXNjYWxlPTEuMCwgdXNlci1zY2FsYWJsZT1ubyIgLz4KICAgIDxzdHlsZT4jbWFwXzg4ZGI1MGJmYjliYTRjNzE5NmE3YzEzOGVmMGNjZjg4IHsKICAgICAgICBwb3NpdGlvbjogcmVsYXRpdmU7CiAgICAgICAgd2lkdGg6IDEwMC4wJTsKICAgICAgICBoZWlnaHQ6IDEwMC4wJTsKICAgICAgICBsZWZ0OiAwLjAlOwogICAgICAgIHRvcDogMC4wJTsKICAgICAgICB9CiAgICA8L3N0eWxlPgogICAgPHNjcmlwdCBzcmM9Imh0dHBzOi8vY2RuanMuY2xvdWRmbGFyZS5jb20vYWpheC9saWJzL2xlYWZsZXQubWFya2VyY2x1c3Rlci8xLjEuMC9sZWFmbGV0Lm1hcmtlcmNsdXN0ZXIuanMiPjwvc2NyaXB0PgogICAgPGxpbmsgcmVsPSJzdHlsZXNoZWV0IiBocmVmPSJodHRwczovL2NkbmpzLmNsb3VkZmxhcmUuY29tL2FqYXgvbGlicy9sZWFmbGV0Lm1hcmtlcmNsdXN0ZXIvMS4xLjAvTWFya2VyQ2x1c3Rlci5jc3MiLz4KICAgIDxsaW5rIHJlbD0ic3R5bGVzaGVldCIgaHJlZj0iaHR0cHM6Ly9jZG5qcy5jbG91ZGZsYXJlLmNvbS9hamF4L2xpYnMvbGVhZmxldC5tYXJrZXJjbHVzdGVyLzEuMS4wL01hcmtlckNsdXN0ZXIuRGVmYXVsdC5jc3MiLz4KPC9oZWFkPgo8Ym9keT4gICAgCiAgICAKICAgIDxkaXYgY2xhc3M9ImZvbGl1bS1tYXAiIGlkPSJtYXBfODhkYjUwYmZiOWJhNGM3MTk2YTdjMTM4ZWYwY2NmODgiID48L2Rpdj4KPC9ib2R5Pgo8c2NyaXB0PiAgICAKICAgIAogICAgCiAgICAgICAgdmFyIGJvdW5kcyA9IG51bGw7CiAgICAKCiAgICB2YXIgbWFwXzg4ZGI1MGJmYjliYTRjNzE5NmE3YzEzOGVmMGNjZjg4ID0gTC5tYXAoCiAgICAgICAgJ21hcF84OGRiNTBiZmI5YmE0YzcxOTZhN2MxMzhlZjBjY2Y4OCcsIHsKICAgICAgICBjZW50ZXI6IFs0Mi4zMywgLTcwLjkzNV0sCiAgICAgICAgem9vbTogMTEsCiAgICAgICAgbWF4Qm91bmRzOiBib3VuZHMsCiAgICAgICAgbGF5ZXJzOiBbXSwKICAgICAgICB3b3JsZENvcHlKdW1wOiBmYWxzZSwKICAgICAgICBjcnM6IEwuQ1JTLkVQU0czODU3LAogICAgICAgIHpvb21Db250cm9sOiB0cnVlLAogICAgICAgIH0pOwoKICAgIAogICAgCiAgICB2YXIgdGlsZV9sYXllcl9kZDNiOGU0ZDdhODA0ODkxOWQ3MWNlN2Q0ZTRjMTYwNSA9IEwudGlsZUxheWVyKAogICAgICAgICdodHRwczovL3tzfS50aWxlLm9wZW5zdHJlZXRtYXAub3JnL3t6fS97eH0ve3l9LnBuZycsCiAgICAgICAgewogICAgICAgICJhdHRyaWJ1dGlvbiI6IG51bGwsCiAgICAgICAgImRldGVjdFJldGluYSI6IGZhbHNlLAogICAgICAgICJtYXhOYXRpdmVab29tIjogMTgsCiAgICAgICAgIm1heFpvb20iOiAxOCwKICAgICAgICAibWluWm9vbSI6IDAsCiAgICAgICAgIm5vV3JhcCI6IGZhbHNlLAogICAgICAgICJvcGFjaXR5IjogMSwKICAgICAgICAic3ViZG9tYWlucyI6ICJhYmMiLAogICAgICAgICJ0bXMiOiBmYWxzZQp9KS5hZGRUbyhtYXBfODhkYjUwYmZiOWJhNGM3MTk2YTdjMTM4ZWYwY2NmODgpOwogICAgCiAgICAgICAgICAgIHZhciBtYWNyb19lbGVtZW50XzY1YzcwNDc2MTcxZjQzMTY5Yjc2ZThlYmZkMjY1OWZiID0gTC50aWxlTGF5ZXIud21zKAogICAgICAgICAgICAgICAgJ2h0dHA6Ly9vb3Muc29lc3QuaGF3YWlpLmVkdS90aHJlZGRzL3dtcy9oaW9vcy9zYXRlbGxpdGUvZGh3XzVrbScsCiAgICAgICAgICAgICAgICB7CiAgImF0dHJpYnV0aW9uIjogIlBhY0lPT1MgVERTIiwKICAiZm9ybWF0IjogImltYWdlL3BuZyIsCiAgImxheWVycyI6ICJDUldfU1NUIiwKICAic3R5bGVzIjogIiIsCiAgInRyYW5zcGFyZW50IjogdHJ1ZSwKICAidmVyc2lvbiI6ICIxLjEuMSIKfQogICAgICAgICAgICAgICAgKS5hZGRUbyhtYXBfODhkYjUwYmZiOWJhNGM3MTk2YTdjMTM4ZWYwY2NmODgpOwoKICAgICAgICAKICAgIAogICAgICAgICAgICAgICAgdmFyIHBvbHlfbGluZV9mMTEwYTAzMGM4NTc0MDAyYTk1YmU0MTlkMmI3NGYwOCA9IEwucG9seWxpbmUoCiAgICAgICAgICAgICAgICAgICAgW1s0Mi4wMywgLTcxLjNdLCBbNDIuMDMsIC03MC41N10sIFs0Mi42MywgLTcwLjU3XSwgWzQyLjYzLCAtNzEuM10sIFs0Mi4wMywgLTcxLjNdXSwKICAgICAgICAgICAgICAgICAgICB7CiAgImJ1YmJsaW5nTW91c2VFdmVudHMiOiB0cnVlLAogICJjb2xvciI6ICIjRkYwMDAwIiwKICAiZGFzaEFycmF5IjogbnVsbCwKICAiZGFzaE9mZnNldCI6IG51bGwsCiAgImZpbGwiOiBmYWxzZSwKICAiZmlsbENvbG9yIjogIiNGRjAwMDAiLAogICJmaWxsT3BhY2l0eSI6IDAuMiwKICAiZmlsbFJ1bGUiOiAiZXZlbm9kZCIsCiAgImxpbmVDYXAiOiAicm91bmQiLAogICJsaW5lSm9pbiI6ICJyb3VuZCIsCiAgIm5vQ2xpcCI6IGZhbHNlLAogICJvcGFjaXR5IjogMC45LAogICJzbW9vdGhGYWN0b3IiOiAxLjAsCiAgInN0cm9rZSI6IHRydWUsCiAgIndlaWdodCI6IDIKfQogICAgICAgICAgICAgICAgICAgICkKICAgICAgICAgICAgICAgICAgICAuYWRkVG8obWFwXzg4ZGI1MGJmYjliYTRjNzE5NmE3YzEzOGVmMGNjZjg4KTsKICAgICAgICAgICAgCiAgICAKICAgICAgICAgICAgdmFyIG1hcmtlcl9jbHVzdGVyXzZkMzE4MmM5MzcwMjQ1MzI5YWUzYzA1YWI2YzlkMWY4ID0gTC5tYXJrZXJDbHVzdGVyR3JvdXAoe30pOwogICAgICAgICAgICBtYXBfODhkYjUwYmZiOWJhNGM3MTk2YTdjMTM4ZWYwY2NmODguYWRkTGF5ZXIobWFya2VyX2NsdXN0ZXJfNmQzMTgyYzkzNzAyNDUzMjlhZTNjMDVhYjZjOWQxZjgpOwogICAgICAgICAgICAKICAgIAogICAgICAgIHZhciBtYXJrZXJfYzExM2M4ZDNjNDhmNGVhMDljYzIzMTI5ZDcwNGY1ZjAgPSBMLm1hcmtlcigKICAgICAgICAgICAgWzQyLjMzNzIxMDQwODYxOTQ3LCAtNzAuNjg5ODkwNjYwMDk5NDJdLAogICAgICAgICAgICB7CiAgICAgICAgICAgICAgICBpY29uOiBuZXcgTC5JY29uLkRlZmF1bHQoKQogICAgICAgICAgICAgICAgfQogICAgICAgICAgICApLmFkZFRvKG1hcmtlcl9jbHVzdGVyXzZkMzE4MmM5MzcwMjQ1MzI5YWUzYzA1YWI2YzlkMWY4KTsKICAgICAgICAKICAgIAogICAgICAgICAgICB2YXIgcG9wdXBfMWM2YWQxOWY3ZWUzNGU1M2I0Y2Y5NDhiYmNjODViM2MgPSBMLnBvcHVwKHttYXhXaWR0aDogJzMwMCcKICAgICAgICAgICAgCiAgICAgICAgICAgIH0pOwoKICAgICAgICAgICAgCiAgICAgICAgICAgICAgICB2YXIgaHRtbF83MTI3MmI5YjVlYzk0MjE4OWFlMzM3OGI3OTUwNzgyNCA9ICQoYDxkaXYgaWQ9Imh0bWxfNzEyNzJiOWI1ZWM5NDIxODlhZTMzNzhiNzk1MDc4MjQiIHN0eWxlPSJ3aWR0aDogMTAwLjAlOyBoZWlnaHQ6IDEwMC4wJTsiPltIaXN0b3J5X0Jlc3QjZmlsbG1pc21hdGNoXTogQk9TVE9OIDE2IE5NIEVhc3Qgb2YgQm9zdG9uLCBNQTwvZGl2PmApWzBdOwogICAgICAgICAgICAgICAgcG9wdXBfMWM2YWQxOWY3ZWUzNGU1M2I0Y2Y5NDhiYmNjODViM2Muc2V0Q29udGVudChodG1sXzcxMjcyYjliNWVjOTQyMTg5YWUzMzc4Yjc5NTA3ODI0KTsKICAgICAgICAgICAgCgogICAgICAgICAgICBtYXJrZXJfYzExM2M4ZDNjNDhmNGVhMDljYzIzMTI5ZDcwNGY1ZjAuYmluZFBvcHVwKHBvcHVwXzFjNmFkMTlmN2VlMzRlNTNiNGNmOTQ4YmJjYzg1YjNjKQogICAgICAgICAgICA7CgogICAgICAgICAgICAKICAgICAgICAKICAgIAogICAgICAgIHZhciBtYXJrZXJfMjcxYmI5NDE3NTA3NGVjMzg3YWY1YjFmNzhkZTM4ZDQgPSBMLm1hcmtlcigKICAgICAgICAgICAgWzQyLjM0MTMyNzY2NzIzNjMzLCAtNzAuNjQ4MzE1NDI5Njg3NV0sCiAgICAgICAgICAgIHsKICAgICAgICAgICAgICAgIGljb246IG5ldyBMLkljb24uRGVmYXVsdCgpCiAgICAgICAgICAgICAgICB9CiAgICAgICAgICAgICkuYWRkVG8obWFya2VyX2NsdXN0ZXJfNmQzMTgyYzkzNzAyNDUzMjlhZTNjMDVhYjZjOWQxZjgpOwogICAgICAgIAogICAgCiAgICAgICAgICAgIHZhciBwb3B1cF9kOTIwMjAzYjI2YjM0NWYwOTM1MDk4NThiMjllNDE5MCA9IEwucG9wdXAoe21heFdpZHRoOiAnMzAwJwogICAgICAgICAgICAKICAgICAgICAgICAgfSk7CgogICAgICAgICAgICAKICAgICAgICAgICAgICAgIHZhciBodG1sXzRkNjdhNTEzZDQ3YTQ1NGJiN2IyYjcyMDAzNTg4M2M0ID0gJChgPGRpdiBpZD0iaHRtbF80ZDY3YTUxM2Q0N2E0NTRiYjdiMmI3MjAwMzU4ODNjNCIgc3R5bGU9IndpZHRoOiAxMDAuMCU7IGhlaWdodDogMTAwLjAlOyI+W05FQ09GU19GVkNPTV9PQ0VBTl9NQVNTQkFZX0ZPUkVDQVNUXTogQk9TVE9OIDE2IE5NIEVhc3Qgb2YgQm9zdG9uLCBNQTwvZGl2PmApWzBdOwogICAgICAgICAgICAgICAgcG9wdXBfZDkyMDIwM2IyNmIzNDVmMDkzNTA5ODU4YjI5ZTQxOTAuc2V0Q29udGVudChodG1sXzRkNjdhNTEzZDQ3YTQ1NGJiN2IyYjcyMDAzNTg4M2M0KTsKICAgICAgICAgICAgCgogICAgICAgICAgICBtYXJrZXJfMjcxYmI5NDE3NTA3NGVjMzg3YWY1YjFmNzhkZTM4ZDQuYmluZFBvcHVwKHBvcHVwX2Q5MjAyMDNiMjZiMzQ1ZjA5MzUwOTg1OGIyOWU0MTkwKQogICAgICAgICAgICA7CgogICAgICAgICAgICAKICAgICAgICAKICAgIAogICAgICAgIHZhciBtYXJrZXJfMDNlYjY4Mzc0MTgxNGQ1NWExOTNlZGEzZjM5OWI0M2UgPSBMLm1hcmtlcigKICAgICAgICAgICAgWzQyLjMzNzIxMDQwODYxOTQ3LCAtNzAuNjg5ODkwNjYwMDk5NDJdLAogICAgICAgICAgICB7CiAgICAgICAgICAgICAgICBpY29uOiBuZXcgTC5JY29uLkRlZmF1bHQoKQogICAgICAgICAgICAgICAgfQogICAgICAgICAgICApLmFkZFRvKG1hcmtlcl9jbHVzdGVyXzZkMzE4MmM5MzcwMjQ1MzI5YWUzYzA1YWI2YzlkMWY4KTsKICAgICAgICAKICAgIAogICAgICAgICAgICB2YXIgcG9wdXBfZDVmMzJkMjNjODVlNDNhZDhmZWJmN2MwZjFmM2Q2ODUgPSBMLnBvcHVwKHttYXhXaWR0aDogJzMwMCcKICAgICAgICAgICAgCiAgICAgICAgICAgIH0pOwoKICAgICAgICAgICAgCiAgICAgICAgICAgICAgICB2YXIgaHRtbF8zOTkwMDljYjBmZWQ0MTk5OTdkN2NkMmUxNTFiYTM4NSA9ICQoYDxkaXYgaWQ9Imh0bWxfMzk5MDA5Y2IwZmVkNDE5OTk3ZDdjZDJlMTUxYmEzODUiIHN0eWxlPSJ3aWR0aDogMTAwLjAlOyBoZWlnaHQ6IDEwMC4wJTsiPltBdmVyYWdlc19CZXN0I2ZpbGxtaXNtYXRjaF06IEJPU1RPTiAxNiBOTSBFYXN0IG9mIEJvc3RvbiwgTUE8L2Rpdj5gKVswXTsKICAgICAgICAgICAgICAgIHBvcHVwX2Q1ZjMyZDIzYzg1ZTQzYWQ4ZmViZjdjMGYxZjNkNjg1LnNldENvbnRlbnQoaHRtbF8zOTkwMDljYjBmZWQ0MTk5OTdkN2NkMmUxNTFiYTM4NSk7CiAgICAgICAgICAgIAoKICAgICAgICAgICAgbWFya2VyXzAzZWI2ODM3NDE4MTRkNTVhMTkzZWRhM2YzOTliNDNlLmJpbmRQb3B1cChwb3B1cF9kNWYzMmQyM2M4NWU0M2FkOGZlYmY3YzBmMWYzZDY4NSkKICAgICAgICAgICAgOwoKICAgICAgICAgICAgCiAgICAgICAgCiAgICAKICAgICAgICB2YXIgbWFya2VyX2M5M2JiYjJjOWRhYjQ0YzVhZWRjYjcxYmM5ZTJmZDE4ID0gTC5tYXJrZXIoCiAgICAgICAgICAgIFs0Mi4zNTM3MDI4MzA0NDk5NCwgLTcwLjY0MTQwMjcxNDg0NjExXSwKICAgICAgICAgICAgewogICAgICAgICAgICAgICAgaWNvbjogbmV3IEwuSWNvbi5EZWZhdWx0KCkKICAgICAgICAgICAgICAgIH0KICAgICAgICAgICAgKS5hZGRUbyhtYXJrZXJfY2x1c3Rlcl82ZDMxODJjOTM3MDI0NTMyOWFlM2MwNWFiNmM5ZDFmOCk7CiAgICAgICAgCiAgICAKICAgICAgICAgICAgdmFyIHBvcHVwX2I2MzlhOTM2NWU3NzRkNGVhYjk4ZjMxZTljMDE2NjE4ID0gTC5wb3B1cCh7bWF4V2lkdGg6ICczMDAnCiAgICAgICAgICAgIAogICAgICAgICAgICB9KTsKCiAgICAgICAgICAgIAogICAgICAgICAgICAgICAgdmFyIGh0bWxfYjUyNTQ0NWQ4NzcwNGY3MWEzZjY3MTJiOWFmOTc2YjQgPSAkKGA8ZGl2IGlkPSJodG1sX2I1MjU0NDVkODc3MDRmNzFhM2Y2NzEyYjlhZjk3NmI0IiBzdHlsZT0id2lkdGg6IDEwMC4wJTsgaGVpZ2h0OiAxMDAuMCU7Ij5bU0VDT09SQV9OQ1NVX0NOQVBTXTogQk9TVE9OIDE2IE5NIEVhc3Qgb2YgQm9zdG9uLCBNQTwvZGl2PmApWzBdOwogICAgICAgICAgICAgICAgcG9wdXBfYjYzOWE5MzY1ZTc3NGQ0ZWFiOThmMzFlOWMwMTY2MTguc2V0Q29udGVudChodG1sX2I1MjU0NDVkODc3MDRmNzFhM2Y2NzEyYjlhZjk3NmI0KTsKICAgICAgICAgICAgCgogICAgICAgICAgICBtYXJrZXJfYzkzYmJiMmM5ZGFiNDRjNWFlZGNiNzFiYzllMmZkMTguYmluZFBvcHVwKHBvcHVwX2I2MzlhOTM2NWU3NzRkNGVhYjk4ZjMxZTljMDE2NjE4KQogICAgICAgICAgICA7CgogICAgICAgICAgICAKICAgICAgICAKICAgIAogICAgICAgIHZhciBtYXJrZXJfNWJkMWU1M2U3ZTU0NDdiNjgzNThhM2VhNGFjYzUwNTUgPSBMLm1hcmtlcigKICAgICAgICAgICAgWzQyLjM0NjAwMDAwMDAwMDAwNCwgLTcwLjY1MTAwMDAwMDAwMDAxXSwKICAgICAgICAgICAgewogICAgICAgICAgICAgICAgaWNvbjogbmV3IEwuSWNvbi5EZWZhdWx0KCkKICAgICAgICAgICAgICAgIH0KICAgICAgICAgICAgKS5hZGRUbyhtYXBfODhkYjUwYmZiOWJhNGM3MTk2YTdjMTM4ZWYwY2NmODgpOwogICAgICAgIAogICAgCgogICAgICAgICAgICAgICAgdmFyIGljb25fOWUzNDRlMTEyNmI4NDllMzkzNzVkNGRjZDVmOTk4MGMgPSBMLkF3ZXNvbWVNYXJrZXJzLmljb24oewogICAgICAgICAgICAgICAgICAgIGljb246ICdzdGF0cycsCiAgICAgICAgICAgICAgICAgICAgaWNvbkNvbG9yOiAnd2hpdGUnLAogICAgICAgICAgICAgICAgICAgIG1hcmtlckNvbG9yOiAnZ3JlZW4nLAogICAgICAgICAgICAgICAgICAgIHByZWZpeDogJ2dseXBoaWNvbicsCiAgICAgICAgICAgICAgICAgICAgZXh0cmFDbGFzc2VzOiAnZmEtcm90YXRlLTAnCiAgICAgICAgICAgICAgICAgICAgfSk7CiAgICAgICAgICAgICAgICBtYXJrZXJfNWJkMWU1M2U3ZTU0NDdiNjgzNThhM2VhNGFjYzUwNTUuc2V0SWNvbihpY29uXzllMzQ0ZTExMjZiODQ5ZTM5Mzc1ZDRkY2Q1Zjk5ODBjKTsKICAgICAgICAgICAgCiAgICAKICAgICAgICAgICAgdmFyIHBvcHVwXzAzMGVlOWI4YTRmMzQ2YmRiNGYwOWRjOWY1NjQ3N2E3ID0gTC5wb3B1cCh7bWF4V2lkdGg6ICcyNjUwJwogICAgICAgICAgICAKICAgICAgICAgICAgfSk7CgogICAgICAgICAgICAKICAgICAgICAgICAgICAgIHZhciBpX2ZyYW1lXzYzYWRiMmU1NGVkNDQ4N2JiOWM2YWJkMmIyYmI1YjUzID0gJChgPGlmcmFtZSBzcmM9ImRhdGE6dGV4dC9odG1sO2NoYXJzZXQ9dXRmLTg7YmFzZTY0LENpQWdJQ0FLQ2dvS1BDRkVUME5VV1ZCRklHaDBiV3crQ2p4b2RHMXNJR3hoYm1jOUltVnVJajRLSUNBS0lDQThhR1ZoWkQ0S0lDQWdJQW9nSUNBZ0lDQThiV1YwWVNCamFHRnljMlYwUFNKMWRHWXRPQ0krQ2lBZ0lDQWdJRHgwYVhSc1pUNDBOREF4TXp3dmRHbDBiR1UrQ2lBZ0lDQWdJQW9nSUNBZ0lDQUtJQ0FnSUNBZ0lDQUtJQ0FnSUNBZ0lDQWdJQW9nSUNBZ0lDQWdJRHhzYVc1cklISmxiRDBpYzNSNWJHVnphR1ZsZENJZ2FISmxaajBpYUhSMGNITTZMeTlqWkc0dWNIbGtZWFJoTG05eVp5OWliMnRsYUM5eVpXeGxZWE5sTDJKdmEyVm9MVEV1TUM0MExtMXBiaTVqYzNNaUlIUjVjR1U5SW5SbGVIUXZZM056SWlBdlBnb2dJQ0FnSUNBZ0lBb2dJQ0FnSUNBZ0lBb2dJQ0FnSUNBZ0lDQWdDaUFnSUNBZ0lDQWdQSE5qY21sd2RDQjBlWEJsUFNKMFpYaDBMMnBoZG1GelkzSnBjSFFpSUhOeVl6MGlhSFIwY0hNNkx5OWpaRzR1Y0hsa1lYUmhMbTl5Wnk5aWIydGxhQzl5Wld4bFlYTmxMMkp2YTJWb0xURXVNQzQwTG0xcGJpNXFjeUkrUEM5elkzSnBjSFErQ2lBZ0lDQWdJQ0FnUEhOamNtbHdkQ0IwZVhCbFBTSjBaWGgwTDJwaGRtRnpZM0pwY0hRaVBnb2dJQ0FnSUNBZ0lDQWdJQ0JDYjJ0bGFDNXpaWFJmYkc5blgyeGxkbVZzS0NKcGJtWnZJaWs3Q2lBZ0lDQWdJQ0FnUEM5elkzSnBjSFErQ2lBZ0lDQWdJQ0FnQ2lBZ0lDQWdJQW9nSUNBZ0lDQUtJQ0FnSUFvZ0lEd3ZhR1ZoWkQ0S0lDQUtJQ0FLSUNBOFltOWtlVDRLSUNBZ0lBb2dJQ0FnSUNBS0lDQWdJQ0FnSUNBS0lDQWdJQ0FnSUNBZ0lBb2dJQ0FnSUNBZ0lDQWdDaUFnSUNBZ0lDQWdJQ0FnSUFvZ0lDQWdJQ0FnSUNBZ0lDQWdJRHhrYVhZZ1kyeGhjM005SW1KckxYSnZiM1FpSUdsa1BTSTJPRE5oWWpJelpDMDVOVGN4TFRSaU1tSXRPRFUwWWkwMU5qVXhOemsyWTJNeU16SWlJR1JoZEdFdGNtOXZkQzFwWkQwaU1UQXdNaUkrUEM5a2FYWStDaUFnSUNBZ0lDQWdJQ0FnSUFvZ0lDQWdJQ0FnSUNBZ0NpQWdJQ0FnSUNBZ0NpQWdJQ0FnSUFvZ0lDQWdJQ0FLSUNBZ0lDQWdJQ0E4YzJOeWFYQjBJSFI1Y0dVOUltRndjR3hwWTJGMGFXOXVMMnB6YjI0aUlHbGtQU0l4TkRBMklqNEtJQ0FnSUNBZ0lDQWdJSHNpTlRBNE1UTTFPVEl0WXpZMU5pMDBOVE5rTFRoa056a3RNak5qTUdNMllUQXlOV1JrSWpwN0luSnZiM1J6SWpwN0luSmxabVZ5Wlc1alpYTWlPbHQ3SW1GMGRISnBZblYwWlhNaU9uc2lZMkZzYkdKaFkyc2lPbTUxYkd3c0ltUmhkR0VpT25zaWVDSTZleUpmWDI1a1lYSnlZWGxmWHlJNklrRkJRVUZ5YlZkTFpHdEpRVUZQWjJOaFdYQXlVV2RCUVRCSmRITnBibHBEUVVGRE5DdHRLMHRrYTBsQlFVdENjR00wY0RKUlowRkJhVTVvTW1sdVdrTkJRVUozVWpOeFMyUnJTVUZCUm1reVpsbHdNbEZuUVVGUlExZENhVzVhUTBGQlFXOXNTVk5MWkd0SlFVRkNRVVJwU1hBeVVXZEJRU3RJUjB4cGJscERRVUZFWnpSSk5rdGthMGxCUVUxb1VHdHZjREpSWjBGQmMwdzJWbWx1V2tOQlFVTlpURnB0UzJSclNVRkJTVU5qYmtsd01sRm5RVUZoUVhWbmFXNWFRMEZCUWxGbGNVOUxaR3RKUVVGRWFuQndiM0F5VVdkQlFVbEdhWEZwYmxwRFFVRkJTWGcyTWt0a2EwbEJRVkJCTVhOWmNESlJaMEZCTWt0VE1HbHVXa05CUVVSQlJUZHBTMlJyU1VGQlMybERkVFJ3TWxGblFVRnJVRWNyYVc1YVEwRkJRalJaVFV0TFpHdEpRVUZIUkZCNFdYQXlVV2RCUVZORU4wcHBibHBEUVVGQmQzSmplVXRrYTBsQlFVSm5ZekJKY0RKUlowRkJRVWwyVkdsdVdrTkJRVVJ2SzJSaFMyUnJTVUZCVGtKdk1tOXdNbEZuUVVGMVRtWmthVzVhUTBGQlEyZFNkVWRMWkd0SlFVRkphVEUxU1hBeVVXZEJRV05EVkc5cGJscERRVUZDV1dzcmRVdGthMGxCUVVWQlF6YzBjREpSWjBGQlMwaEllV2x1V2tOQlFVRlJORkJYUzJSclNVRkJVR2hQSzFsd01sRm5RVUUwVERNNGFXNWFRMEZCUkVsTVFVTk1aR3RKUVVGTVEySkJOSFF5VVdkQlFXMUJiMGhwTTFwRFFVRkRRV1ZSY1V4a2EwbEJRVWRxYjBSWmRESlJaMEZCVlVaalVta3pXa05CUVVFMGVHaFRUR1JyU1VGQlEwRXhSMGwwTWxGblFVRkRTMUZpYVROYVEwRkJSSGRGYUN0TVpHdEpRVUZPYVVKSmIzUXlVV2RCUVhkUVFXeHBNMXBEUVVGRGIxaDViVXhrYTBsQlFVcEVUMHhKZERKUlowRkJaVVF3ZDJreldrTkJRVUpuY2tSUFRHUnJTVUZCUldkaVRqUjBNbEZuUVVGTlNXODJhVE5hUTBGQlFWa3JWREpNWkd0SlFVRkJRbTlSV1hReVVXZEJRVFpPV2tWcE0xcERRVUZFVVZKVmFVeGthMGxCUVV4cE1GTTBkREpSWjBGQmIwTk9VR2t6V2tOQlFVTkphMnhMVEdSclNVRkJTRUZDVm05ME1sRm5RVUZYU0VKYWFUTmFRMEZCUWtFek1YbE1aR3RKUVVGRGFFOVpTWFF5VVdkQlFVVk1NV3BwTTFwRFFVRkVORXN5WlV4a2EwbEJRVTlEWVdGdmRESlJaMEZCZVVGc2RXa3pXa05CUVVOM1pVaEhUR1JyU1VGQlNtcHVaRWwwTWxGblFVRm5SbG8wYVROYVEwRkJRbTk0V0hWTVpHdEpRVUZHUVRCbU5IUXlVV2RCUVU5TFQwTnBNMXBEUVVGQlowVnZZVXhrYTBsQlFVRnBRbWxaZERKUlowRkJPRThyVFdreldrTkJRVVJaV0hCRFRHUnJTVUZCVFVST2F6UjBNbEZuUVVGeFJIbFlhVE5hUTBGQlExRnhOWEZNWkd0SlFVRklaMkZ1YjNReVVXZEJRVmxKYldocE0xcERRVUZDU1N0TFUweGthMGxCUVVSQ2JuRkpkREpSWjBGQlIwNWhjbWt6V2tOQlFVRkJVbUVyVEdSclNVRkJUMmw2YzI5ME1sRm5RVUV3UTBzeWFUTmFRMEZCUXpSclltMU1aR3RKUVVGTFFVRjJXWFF5VVdkQlFXbEhMMEZwTTFwRFFVRkNkek56VDB4a2EwbEJRVVpvVG5nMGRESlJaMEZCVVV4NlMya3pXa05CUVVGdlN6ZzJUR1JyU1VGQlFrTmhNRmwwTWxGblFVRXJRV3BXYVROYVEwRkJSR2RrT1dsTVpHdEpRVUZOYW0weU5IUXlVV2RCUVhOR1dHWnBNMXBEUVVGRFdYaFBTMHhrYTBsQlFVbEJlalZ2ZERKUlowRkJZVXRNY0dreldrTkJRVUpSUldVeVRHUnJTVUZCUkdsQk9FbDBNbEZuUVVGSlR5OTZhVE5hUTBGQlFVbFlkbVZNWkd0SlFVRlFSRTByYjNReVVXZEJRVEpFZGl0cE0xcERRVUZFUVhGblIwMWthMGxCUVV0bldrSlplREpSWjBGQmEwbG5TV3BJV2tOQlFVSTBPWGQxVFdSclNVRkJSMEp0UkRSNE1sRm5RVUZUVGxWVGFraGFRMEZCUVhkU1FtRk5aR3RKUVVGQ2FYcEhXWGd5VVdkQlFVRkRTV1JxU0ZwRFFVRkViMnREUTAxa2EwbEJRVTVFTDBrMGVESlJaMEZCZFVjMGJtcElXa05CUVVObk0xTnhUV1JyU1VGQlNXaE5URzk0TWxGblFVRmpUSE40YWtoYVEwRkJRbGxMYWxkTlpHdEpRVUZGUTFwUFNYZ3lVV2RCUVV0Qlp6aHFTRnBEUVVGQlVXUjZLMDFrYTBsQlFWQnFiRkZ2ZURKUlp6MDlJaXdpWkhSNWNHVWlPaUptYkc5aGREWTBJaXdpYzJoaGNHVWlPbHN4TkRCZGZTd2llU0k2ZXlKZlgyNWtZWEp5WVhsZlh5STZJbHB0V20xYWJWcHRSV3RDYlZwdFdtMWFiVmxUVVVkYWJWcHRXbTFhYUVwQldtMWFiVnB0V20xRmEwSnRXbTFhYlZwdFdWTlJSMXB0V20xYWJWcG9Ta0ZhYlZwdFdtMWFiVVZyUW0xYWJWcHRXbTFaVTFGSFdtMWFiVnB0V21oS1FWcHRXbTFhYlZwdFJXdENiVnB0V20xYWJWbFRVVUZCUVVGQlFVRkJRa3BCUVVGQlFVRkJRVUZGYTBGQlFVRkJRVUZCUVZOUlFVRkJRVUZCUVVGQ1NrRkJRVUZCUVVGQlFVVnJRVUZCUVVGQlFVRkJVMUZCUVVGQlFVRkJRVUpLUVVGQlFVRkJRVUZCUld0QlFVRkJRVUZCUVVGVFVVRkJRVUZCUVVGQlFrcEJRVUZCUVVGQlFVRkZhMEZCUVVGQlFVRkJRVk5SUVVGQlFVRkJRVUZDU2tGQlFVRkJRVUZCUVVWclFVRkJRVUZCUVVGQlUxRkJRVUZCUVVGQlFVSktRVzF3YlZwdFdtMWFSVlZEWVcxYWJWcHRXbXRTVVVweFdtMWFiVnB0VWtaQmJYQnRXbTFhYlZwRlZVTmhiVnB0V20xYWExSlJTbkZhYlZwdFdtMVJNVUY2WTNwTmVrMTZUVVJGUkU1NlRYcE5lazEzVFZGTk0wMTZUWHBOZWtGNFFYcGplazE2VFhwTlJFVkVUbnBOZWsxNlRYZE5VVTB6VFhwTmVrMTZRWGhCZW1ONlRYcE5lazFFUlVST2VrMTZUWHBOZDAxUlNuRmFiVnB0V20xUk1VRk5lazE2VFhwTmVrUXdSRTU2VFhwTmVrMTNVVkZFVFhwTmVrMTZUWGhHUVUxNlRYcE5lazE2UlZWQmVrMTZUWHBOZWsxU1VVUk5lazE2VFhwTmVFWkJUWHBOZWsxNlRYcEZWVUY2VFhwTmVrMTZUVkpSUkUxNlRYcE5lazE0UmtGTmVrMTZUWHBOZWtWVlFYcE5lazE2VFhwTlVsRkVUWHBOZWsxNlRYaEdRVTE2VFhwTmVrMTZSVlZCZWsxNlRYcE5lazFTVVVSTmVrMTZUWHBOZUVaQlRYcE5lazE2VFhwRlZVRjZUWHBOZWsxNlRWSlJRVUZCUVVGQlFVRkNRa0ZOZWsxNlRYcE5la1F3UVVGQlFVRkJRVUZCVVZGTk0wMTZUWHBOZWtKQ1FVMTZUWHBOZWsxNlJWVkJlazE2VFhwTmVrMVNVVXB4V20xYWJWcHRVa1pCYlhCdFdtMWFiVnBGVlVGNlRYcE5lazE2VFZKUlNuRmFiVnB0V20xU1JrRnRjRzFhYlZwdFdrVlZRWHBOZWsxNlRYcE5VbEZLY1ZwdFdtMWFiVkpHUVcxd2JWcHRXbTFhUlZWQmVrMTZUWHBOZWsxU1VVUk5lazE2VFhwTmVFWkJUWHBOZWsxNlRYcEZWVUY2VFhwTmVrMTZUVkpSUkUxNlRYcE5lazE0UmtGTmVrMTZUWHBOZWtWVlFYcE5lazE2VFhwTlVsRkVUWHBOZWsxNlRYaEdRWHBqZWsxNlRYcE5SVVZFVG5wTmVrMTZUWGRSVVVkYWJWcHRXbTFhYUVKQldtMWFiVnB0V20xRlJVSnRXbTFhYlZwdFdWRlJRVUZCUVVGQlFVRkNRa0ZhYlZwdFdtMWFiVVZGUkU1NlRYcE5lazEzVVZGS2NWcHRXbTFhYlZKR1FVRkJRVUZCUVVGQlJXdENiVnB0V20xYWJWbFRVVUZCUVVGQlFVRkJRa3BCUVVGQlFVRkJRVUZGYTBOaGJWcHRXbTFhYTFKUlJFMTZUWHBOZWsxNFJrRk5lazE2VFhwTmVrVlZRWHBOZWsxNlRYcE5VbEZFVFhwTmVrMTZUWGhHUVUxNlRYcE5lazE2UlZWQmVrMTZUWHBOZWsxU1VVUk5lazE2VFhwTmVFWkJUWHBOZWsxNlRYcEZWVVJPZWsxNlRYcE5kMUZSUjFwdFdtMWFiVnBvUWtGQlFVRkJRVUZCUVVWRlFYcE5lazE2VFhwTlVGRkhXbTFhYlZwdFdtYzFRVnB0V20xYWJWcHRSR3RDYlZwdFdtMWFiVmxQVVVSTmVrMTZUWHBOZHpsQlFVRkJRVUZCUVVGRlJVSnRXbTFhYlZwdFdWRlJSRTE2VFhwTmVrMTRSa0ZCUVVGQlFVRkJRVVZyUW0xYWJWcHRXbTFaVTFGQlFVRkJRVUZCUVVKS1FVRkJRVUZCUVVGQlJXdERZVzFhYlZwdFdtdFNVVVJOZWsxNlRYcE5lRVpCVFhwTmVrMTZUWHBGVlVST2VrMTZUWHBOZDFGUlRUTk5lazE2VFhwQ1FrRk5lazE2VFhwTmVrVlZSRTU2VFhwTmVrMTNVVkZCUVVGQlFVRkJRVUpDUVVGQlFVRkJRVUZCUlVWQ2JWcHRXbTFhYlZsUlVVZGFiVnB0V20xYWFFSkJXbTFhYlZwdFdtMUZSVUp0V20xYWJWcHRXVkZSUjFwdFdtMWFiVnBvUWtGYWJWcHRXbTFhYlVWRlFtMWFiVnB0V20xWlVWRk5NMDE2VFhwTmVrSkNRVTE2VFhwTmVrMTZSVlZCZWsxNlRYcE5lazFTVVVSTmVrMTZUWHBOZUVaQmJYQnRXbTFhYlZwRlZVTmhiVnB0V20xYWExSlJRVDA5SWl3aVpIUjVjR1VpT2lKbWJHOWhkRFkwSWl3aWMyaGhjR1VpT2xzeE5EQmRmWDBzSW5ObGJHVmpkR1ZrSWpwN0ltbGtJam9pTVRFM055SXNJblI1Y0dVaU9pSlRaV3hsWTNScGIyNGlmU3dpYzJWc1pXTjBhVzl1WDNCdmJHbGplU0k2ZXlKcFpDSTZJakV4TnpZaUxDSjBlWEJsSWpvaVZXNXBiMjVTWlc1a1pYSmxjbk1pZlgwc0ltbGtJam9pTVRFeE5pSXNJblI1Y0dVaU9pSkRiMngxYlc1RVlYUmhVMjkxY21ObEluMHNleUpoZEhSeWFXSjFkR1Z6SWpwN0ltTmhiR3hpWVdOcklqcHVkV3hzTENKeVpXNWtaWEpsY25NaU9sdDdJbWxrSWpvaU1UQTJNU0lzSW5SNWNHVWlPaUpIYkhsd2FGSmxibVJsY21WeUluMWRMQ0owYjI5c2RHbHdjeUk2VzFzaVRtRnRaU0lzSWtocGMzUnZjbmxmUW1WemRDTm1hV3hzYldsemJXRjBZMmdpWFN4YklrSnBZWE1pTENJdE1DNDNOQ0pkTEZzaVUydHBiR3dpTENJd0xqTXlJbDFkZlN3aWFXUWlPaUl4TURnMElpd2lkSGx3WlNJNklraHZkbVZ5Vkc5dmJDSjlMSHNpWVhSMGNtbGlkWFJsY3lJNmUzMHNJbWxrSWpvaU1UQXdPQ0lzSW5SNWNHVWlPaUpNYVc1bFlYSlRZMkZzWlNKOUxIc2lZWFIwY21saWRYUmxjeUk2ZTMwc0ltbGtJam9pTVRFM055SXNJblI1Y0dVaU9pSlRaV3hsWTNScGIyNGlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2liR0ZpWld3aU9uc2lkbUZzZFdVaU9pSklhWE4wYjNKNVgwSmxjM1FqWm1sc2JHMXBjMjFoZEdOb0luMHNJbkpsYm1SbGNtVnljeUk2VzNzaWFXUWlPaUl4TURZeElpd2lkSGx3WlNJNklrZHNlWEJvVW1WdVpHVnlaWElpZlYxOUxDSnBaQ0k2SWpFd09ETWlMQ0owZVhCbElqb2lUR1ZuWlc1a1NYUmxiU0o5TEhzaVlYUjBjbWxpZFhSbGN5STZleUpzWVdKbGJDSTZleUoyWVd4MVpTSTZJbE5GUTA5UFVrRmZUa05UVlY5RFRrRlFVeTV1WXlObWFXeHNiV2x6YldGMFkyZ2lmU3dpY21WdVpHVnlaWEp6SWpwYmV5SnBaQ0k2SWpFeE5URWlMQ0owZVhCbElqb2lSMng1Y0doU1pXNWtaWEpsY2lKOVhYMHNJbWxrSWpvaU1URTNPU0lzSW5SNWNHVWlPaUpNWldkbGJtUkpkR1Z0SW4wc2V5SmhkSFJ5YVdKMWRHVnpJanA3ZlN3aWFXUWlPaUl4TURFd0lpd2lkSGx3WlNJNklreHBibVZoY2xOallXeGxJbjBzZXlKaGRIUnlhV0oxZEdWeklqcDdJbXhwYm1WZllXeHdhR0VpT2pBdU5qVXNJbXhwYm1WZlkyRndJam9pY205MWJtUWlMQ0pzYVc1bFgyTnZiRzl5SWpvaUkyRmxZemRsT0NJc0lteHBibVZmYW05cGJpSTZJbkp2ZFc1a0lpd2liR2x1WlY5M2FXUjBhQ0k2TlN3aWVDSTZleUptYVdWc1pDSTZJbmdpZlN3aWVTSTZleUptYVdWc1pDSTZJbmtpZlgwc0ltbGtJam9pTVRBMU9TSXNJblI1Y0dVaU9pSk1hVzVsSW4wc2V5SmhkSFJ5YVdKMWRHVnpJanA3ZlN3aWFXUWlPaUl4TWpBM0lpd2lkSGx3WlNJNklsVnVhVzl1VW1WdVpHVnlaWEp6SW4wc2V5SmhkSFJ5YVdKMWRHVnpJanA3SW14cGJtVmZZMkZ3SWpvaWNtOTFibVFpTENKc2FXNWxYMk52Ykc5eUlqb2lZM0pwYlhOdmJpSXNJbXhwYm1WZmFtOXBiaUk2SW5KdmRXNWtJaXdpYkdsdVpWOTNhV1IwYUNJNk5Td2llQ0k2ZXlKbWFXVnNaQ0k2SW5naWZTd2llU0k2ZXlKbWFXVnNaQ0k2SW5raWZYMHNJbWxrSWpvaU1URXhOeUlzSW5SNWNHVWlPaUpNYVc1bEluMHNleUpoZEhSeWFXSjFkR1Z6SWpwN0ltTmhiR3hpWVdOcklqcHVkV3hzTENKeVpXNWtaWEpsY25NaU9sdDdJbWxrSWpvaU1URTFNU0lzSW5SNWNHVWlPaUpIYkhsd2FGSmxibVJsY21WeUluMWRMQ0owYjI5c2RHbHdjeUk2VzFzaVRtRnRaU0lzSWxORlEwOVBVa0ZmVGtOVFZWOURUa0ZRVXk1dVl5Tm1hV3hzYldsemJXRjBZMmdpWFN4YklrSnBZWE1pTENJdE1DNHpPQ0pkTEZzaVUydHBiR3dpTENJd0xqSTFJbDFkZlN3aWFXUWlPaUl4TVRnd0lpd2lkSGx3WlNJNklraHZkbVZ5Vkc5dmJDSjlMSHNpWVhSMGNtbGlkWFJsY3lJNmV5SnNhVzVsWDJGc2NHaGhJam93TGpFc0lteHBibVZmWTJGd0lqb2ljbTkxYm1RaUxDSnNhVzVsWDJOdmJHOXlJam9pSXpGbU56ZGlOQ0lzSW14cGJtVmZhbTlwYmlJNkluSnZkVzVrSWl3aWJHbHVaVjkzYVdSMGFDSTZOU3dpZUNJNmV5Sm1hV1ZzWkNJNkluZ2lmU3dpZVNJNmV5Sm1hV1ZzWkNJNklua2lmWDBzSW1sa0lqb2lNVEV4T0NJc0luUjVjR1VpT2lKTWFXNWxJbjBzZXlKaGRIUnlhV0oxZEdWeklqcDdJbU5oYkd4aVlXTnJJanB1ZFd4c2ZTd2lhV1FpT2lJeE1EQTBJaXdpZEhsd1pTSTZJa1JoZEdGU1lXNW5aVEZrSW4wc2V5SmhkSFJ5YVdKMWRHVnpJanA3SW1OaGJHeGlZV05ySWpwdWRXeHNMQ0prWVhSaElqcDdJbmdpT25zaVgxOXVaR0Z5Y21GNVgxOGlPaUpCUVVGQmNtMVhTMlJyU1VGQlQyZGpZVmx3TWxGblFVRXdTWFJ6YVc1YVEwRkJRelFyYlN0TFpHdEpRVUZMUW5Cak5IQXlVV2RCUVdsT2FESnBibHBEUVVGQ2QxSXpjVXRrYTBsQlFVWnBNbVpaY0RKUlowRkJVVU5YUW1sdVdrTkJRVUZ2YkVsVFMyUnJTVUZCUWtGRWFVbHdNbEZuUVVFclNFZE1hVzVhUTBGQlJHYzBTVFpMWkd0SlFVRk5hRkJyYjNBeVVXZEJRWE5NTmxacGJscERRVUZEV1V4YWJVdGthMGxCUVVsRFkyNUpjREpSWjBGQllVRjFaMmx1V2tOQlFVSlJaWEZQUzJSclNVRkJSR3B3Y0c5d01sRm5RVUZKUm1seGFXNWFRMEZCUVVsNE5qSkxaR3RKUVVGUVFURnpXWEF5VVdkQlFUSkxVekJwYmxwRFFVRkVRVVUzYVV0a2EwbEJRVXRwUTNVMGNESlJaMEZCYTFCSEsybHVXa05CUVVJMFdVMUxTMlJyU1VGQlIwUlFlRmx3TWxGblFVRlRSRGRLYVc1YVEwRkJRWGR5WTNsTFpHdEpRVUZDWjJNd1NYQXlVV2RCUVVGSmRsUnBibHBEUVVGRWJ5dGtZVXRrYTBsQlFVNUNiekp2Y0RKUlowRkJkVTVtWkdsdVdrTkJRVU5uVW5WSFMyUnJTVUZCU1dreE5VbHdNbEZuUVVGalExUnZhVzVhUTBGQlFsbHJLM1ZMWkd0SlFVRkZRVU0zTkhBeVVXZEJRVXRJU0hscGJscERRVUZCVVRSUVYwdGthMGxCUVZCb1R5dFpjREpSWjBGQk5Fd3pPR2x1V2tOQlFVUkpURUZEVEdSclNVRkJURU5pUVRSME1sRm5RVUZ0UVc5SWFUTmFRMEZCUTBGbFVYRk1aR3RKUVVGSGFtOUVXWFF5VVdkQlFWVkdZMUpwTTFwRFFVRkJOSGhvVTB4a2EwbEJRVU5CTVVkSmRESlJaMEZCUTB0Ullta3pXa05CUVVSM1JXZ3JUR1JyU1VGQlRtbENTVzkwTWxGblFVRjNVRUZzYVROYVEwRkJRMjlZZVcxTVpHdEpRVUZLUkU5TVNYUXlVV2RCUVdWRU1IZHBNMXBEUVVGQ1ozSkVUMHhrYTBsQlFVVm5ZazQwZERKUlowRkJUVWx2Tm1reldrTkJRVUZaSzFReVRHUnJTVUZCUVVKdlVWbDBNbEZuUVVFMlRscEZhVE5hUTBGQlJGRlNWV2xNWkd0SlFVRk1hVEJUTkhReVVXZEJRVzlEVGxCcE0xcERRVUZEU1d0c1MweGthMGxCUVVoQlFsWnZkREpSWjBGQlYwaENXbWt6V2tOQlFVSkJNekY1VEdSclNVRkJRMmhQV1VsME1sRm5RVUZGVERGcWFUTmFRMEZCUkRSTE1tVk1aR3RKUVVGUFEyRmhiM1F5VVdkQlFYbEJiSFZwTTFwRFFVRkRkMlZJUjB4a2EwbEJRVXBxYm1SSmRESlJaMEZCWjBaYU5Ha3pXa05CUVVKdmVGaDFUR1JyU1VGQlJrRXdaalIwTWxGblFVRlBTMDlEYVROYVEwRkJRV2RGYjJGTVpHdEpRVUZCYVVKcFdYUXlVV2RCUVRoUEswMXBNMXBEUVVGRVdWaHdRMHhrYTBsQlFVMUVUbXMwZERKUlowRkJjVVI1V0dreldrTkJRVU5SY1RWeFRHUnJTVUZCU0dkaGJtOTBNbEZuUVVGWlNXMW9hVE5hUTBGQlFra3JTMU5NWkd0SlFVRkVRbTV4U1hReVVXZEJRVWRPWVhKcE0xcERRVUZCUVZKaEsweGthMGxCUVU5cGVuTnZkREpSWjBGQk1FTkxNbWt6V2tOQlFVTTBhMkp0VEdSclNVRkJTMEZCZGxsME1sRm5RVUZwUnk5QmFUTmFRMEZCUW5jemMwOU1aR3RKUVVGR2FFNTROSFF5VVdkQlFWRk1la3RwTTFwRFFVRkJiMHM0Tmt4a2EwbEJRVUpEWVRCWmRESlJaMEZCSzBGcVZta3pXa05CUVVSblpEbHBUR1JyU1VGQlRXcHRNalIwTWxGblFVRnpSbGhtYVROYVEwRkJRMWw0VDB0TVpHdEpRVUZKUVhvMWIzUXlVV2RCUVdGTFRIQnBNMXBEUVVGQ1VVVmxNa3hrYTBsQlFVUnBRVGhKZERKUlowRkJTVTh2ZW1reldrTkJRVUZKV0habFRHUnJTVUZCVUVSTksyOTBNbEZuUVVFeVJIWXJhVE5hUTBGQlJFRnhaMGROWkd0SlFVRkxaMXBDV1hneVVXZEJRV3RKWjBscVNGcERRVUZDTkRsM2RVMWthMGxCUVVkQ2JVUTBlREpSWjBGQlUwNVZVMnBJV2tOQlFVRjNVa0poVFdSclNVRkJRbWw2UjFsNE1sRm5RVUZCUTBsa2FraGFRMEZCUkc5clEwTk5aR3RKUVVGT1JDOUpOSGd5VVdkQlFYVkhORzVxU0ZwRFFVRkRaek5UY1Uxa2EwbEJRVWxvVFV4dmVESlJaMEZCWTB4emVHcElXa05CUVVKWlMycFhUV1JyU1VGQlJVTmFUMGw0TWxGblFVRkxRV2M0YWtoYVEwRkJRVkZrZWl0TlpHdEpRVUZRYW14UmIzZ3lVV2RCUVRSR1VrZHFTRnBEUVVGRVNYY3diVTFrYTBsQlFVeEJlVlJaZURKUlowRkJiVXRHVVdwSVdrTkJRVU5CUlVaVFRXUnJTVUZCUjJndlZqUjRNbEZuUVVGVlR6Vmhha2hhUTBGQlFUUllWalpOWkd0SlFVRkRSRTFaV1hneVVXZEJRVU5FZEd4cVNGcERRVUZFZDNGWGFVMWthMGxCUVU1bldXSkplREpSWjBGQmQwbGtkbXBJV2tOQlFVTnZPVzVMVFdSclNVRkJTa0pzWkc5NE1sRm5RVUZsVGxJMWFraGFRMEZCUW1kUk16Sk5aR3RKUVVGRmFYbG5TWGd5VVdkQlFVMURSMFZxU0ZwRFFVRkJXV3RKWlUxa2EwbEJRVUZFTDJsdmVESlJaMEZCTmtjeVQycElXa05CUVVSUk0wcEhUV1JyU1VGQlRHaE1iRmw0TWxGblFVRnZUSEZaYWtoYVEwRkJRMGxMV25sTlpHdEpRVUZJUTFsdU5IZ3lVV2RCUVZkQlpXcHFTRnBEUVVGQ1FXUnhZVTFrYTBsQlFVTnFiSEZaZURKUlowRkJSVVpUZEdwSVdrTkJRVVEwZDNKRFRXUnJTVUZCVDBGNGRFbDRNbEZuUVVGNVMwTXpha2hhUTBGQlEzZEVOM1ZOWkd0SlFVRkthQ3QyYjNneVVXZEJRV2RQTTBKcVNGcERRVUZDYjFoTlYwMWthMGxCUVVaRVRIbEplREpSWjBGQlQwUnlUV3BJV2tOQlFVRm5jV01yVFdSclNVRkJRV2RaTURSNE1sRm5QVDBpTENKa2RIbHdaU0k2SW1ac2IyRjBOalFpTENKemFHRndaU0k2V3pFNE1sMTlMQ0o1SWpwN0lsOWZibVJoY25KaGVWOWZJam9pUVVGQlFVbEJhR0pFVlVGQlFVRkNaMmxxWjA5UlFVRkJRVUZCUjBkUk5VRkJRVUZCWjBoWU1FUlZRVUZCUVVGQk1UaHpUbEZCUVVGQlEwSXljM2N4UVVGQlFVRTBUa05wUkZWQlFVRkJRMmRFWVdOT1VVRkJRVUZCUTNKMFp6RkJRVUZCUVZGUVZGTkVWVUZCUVVGQloxTm1hMDVSUVVGQlFVVkNZa2huTlVGQlFVRkJiMHBzUTBSclFVRkJRVVJCTTFkQlQxRkJRVUZCUzBKMFptYzFRVUZCUVVGdlJEWlZSR3RCUVVGQlFrRnNZVzlQVVVGQlFVRlBSSGg1VVRWQlFVRkJRVmxHY25WRWEwRkJRVUZDWjBaQ2ExQlJRVUZCUVVORU5GRlJPVUZCUVVGQmQwZFNhVVF3UVVGQlFVRkJaa2d3VUZGQlFVRkJUMEpUYkdjNVFVRkJRVUZaU0RaNFJEQkJRVUZCUkdkMmFXdFFVVUZCUVVGUFFscEpVVGxCUVVGQlFYZEVTVXRFTUVGQlFVRkJaMnhRVFU5UlFVRkJRVWRFVXpKbk5VRkJRVUZCU1VaTVFrUnJRVUZCUVVSbllXRk5UMUZCUVVGQlNVUkRhR2MxUVVGQlFVRjNUbXh4Ukd0QlFVRkJRMmROUmtGUFVVRkJRVUZKUVdaT2R6VkJRVUZCUVVGTE9HZEVhMEZCUVVGRVowZENaMDlSUVVGQlFVZERTRWhCTlVGQlFVRkJTVU5GYkVSclFVRkJRVU5CTDJsdlQxRkJRVUZCU1VOVFMzYzFRVUZCUVVGQlEyZHNSR3RCUVVGQlJFRXJhRkZQVVVGQlFVRk5SRVF2UVRGQlFVRkJRVFJRVUdORVZVRkJRVUZCWjJRM1VVNVJRVUZCUVU5RVMybDNNVUZCUVVGQmQwTm9hVVJWUVVGQlFVSkJNblJOVFZGQlFVRkJUVVJFYTBGNFFVRkJRVUUwUnpWVlJFVkJRVUZCUkdkcmVITk5VVUZCUVVGTFFuQTFkM1JCUVVGQlFWbEpia05ETUVGQlFVRkNRWFJMTkV4UlFVRkJRVWREZDNKQmRFRkJRVUZCTkVWSE9FTXdRVUZCUVVSQmJHUkpURkZCUVVGQlIwSnJPR2QwUVVGQlFVRkJTMUZXUkVWQlFVRkJSRUZrTUZsTlVVRkJRVUZOUWtOblFYaEJRVUZCUVRSTFR5dEVSVUZCUVVGQlFVWlJaMDVSUVVGQlFVZERaVmRuTVVGQlFVRkJRVVZsYUVSVlFVRkJRVUpuZGsxUlRsRkJRVUZCU1VSR2QyY3hRVUZCUVVFMFNUSnBSRlZCUVVGQlEwRndTR05PVVVGQlFVRkRRVWxWVVRGQlFVRkJRVFJOV1hkRVZVRkJRVUZDUVUxS1VVeFJRVUZCUVU5Q00yUm5kRUZCUVVGQlNVaG9VME13UVVGQlFVTkJiVU5yVEZGQlFVRkJTMEZtSzJkd1FVRkJRVUZKUjI1RlEydEJRVUZCUkdkd1dUQkxVVUZCUVVGRlEyaFhkM0JCUVVGQlFWbEJXVEZEYTBGQlFVRkVRV3REUVV0UlFVRkJRVWxDYjBSM2NFRkJRVUZCV1VVNFJFTnJRVUZCUVVKQlFtZFJTMUZCUVVGQlJVTkRSMEZ3UVVGQlFVRkJRVEF5UTJ0QlFVRkJRa0ZqTVdkTFVVRkJRVUZQUkdkblozQkJRVUZCUVdkTGVYTkRhMEZCUVVGRVoybzVXVXRSUVVGQlFVZEViRGQzY0VGQlFVRkJVVTFZZDBOclFVRkJRVU5CVkRsdlMxRkJRVUZCUlVRNGRsRndRVUZCUVVGQlIxTjJRMnRCUVVGQlFrRjNObFZMVVVGQlFVRkRRa2x0ZDNCQlFVRkJRVWxCY1VoRGEwRkJRVUZFWnpoc05FdFJRVUZCUVUxRFlVcEJjRUZCUVVGQlFVOXlXRU5WUVVGQlFVSkJVVWxuU2xGQlFVRkJUVVJ0VVhkc1FVRkJRVUZCUkdkUFExVkJRVUZCUkVGVFpYZEpVVUZCUVVGSFEya3lRV2hCUVVGQlFYZEpSRk5EUlVGQlFVRkJRVmxrTkVsUlFVRkJRVWREVkRoM2FFRkJRVUZCZDBsUlRrTlZRVUZCUVVGQk1GTkZTbEZCUVVGQlRVTllUWGRzUVVGQlFVRTBTVGxGUTFWQlFVRkJRVUZ6TVZWS1VVRkJRVUZIUkZsWloyeEJRVUZCUVRSQ1NtaERWVUZCUVVGQ1FVWldRVXBSUVVGQlFVVkRNbEZCYkVGQlFVRkJVVVpKTkVOVlFVRkJRVVJCYjNaM1RWRkJRVUZCUlVKck9IZDRRVUZCUVVFMFRWaHVSRVZCUVVGQlEwRnZkRVZOVVVGQlFVRkJRbVZ6WjNoQlFVRkJRVmxRSzB0RVJVRkJRVUZEWjNKR01FMVJRVUZCUVVORFlrMUJlRUZCUVVGQldVbGpTRVJGUVVGQlFVSkJaaXRqVEZGQlFVRkJTMEZ5TVhkMFFVRkJRVUZKVUdKU1F6QkJRVUZCUkVGaksyZE1VVUZCUVVGRlJHUkZaM2hCUVVGQlFWRkRjRXBFUlVGQlFVRkVaMkp2VFUxUlFVRkJRVUZCVkhaQmVFRkJRVUZCYjBoMmQwUkZRVUZCUVVGQlNsSjNUbEZCUVVGQlNVRldVRUV4UVVGQlFVRlJVRFZJUkZWQlFVRkJRMEZxZW10T1VVRkJRVUZKUTJ4TFVURkJRVUZCUVVGRWIxaEVWVUZCUVVGQlowdFJSVTVSUVVGQlFVZERRalpCZUVGQlFVRkJiMDFpVDBSRlFVRkJRVU5CVFdKWlRWRkJRVUZCU1VRcmJXZDRRVUZCUVVGblJVTkJSRVZCUVVGQlEyZFhiVmxOVVVGQlFVRkJRVzVWUVhoQlFVRkJRVWxKWkVGRVJVRkJRVUZFUVhCcVkwMVJRVUZCUVU5QlpFOUJlRUZCUVVGQk5FZE9Ra1JGUVVGQlFVUkJaa1pOVFZGQlFVRkJRMFExWWxGNFFVRkJRVUZSVG1WUFJFVkJRVUZCUTJjd1lsRk5VVUZCUVVGTFFVRXdVWGhCUVVGQlFXZEJibVZFUlVGQlFVRkJaMGM1ZDAxUlFVRkJRVTFCYmpCQmVFRkJRVUZCWjBWaVFVUkZRVUZCUVVKblUwczRUVkZCUVVGQlJVSkZjSGQ0UVVGQlFVRTBUbVZ4UkVWQlFVRkJSRUUzY2sxTlVVRkJRVUZMUWl0M1ozaEJRVUZCUVZsQ1JGSkVSVUZCUVVGRVFXdzVORTFSUVVGQlFVRkRjelYzZUVGQlFVRkJaMDl5Y2tSRlFVRkJRVU5CTUN0dlRWRkJRVUZCU1VOd05HZDRRVUZCUVVGQlRtSlZSRVZCUVVGQlJHZEtjMGxOVVVGQlFVRkhRbGh5VVhoQlFVRkJRVmxHWlhSRVJVRkJRVUZDWjFZMk1FMVJRVDA5SWl3aVpIUjVjR1VpT2lKbWJHOWhkRFkwSWl3aWMyaGhjR1VpT2xzeE9ESmRmWDBzSW5ObGJHVmpkR1ZrSWpwN0ltbGtJam9pTVRFeE1TSXNJblI1Y0dVaU9pSlRaV3hsWTNScGIyNGlmU3dpYzJWc1pXTjBhVzl1WDNCdmJHbGplU0k2ZXlKcFpDSTZJakV4TVRBaUxDSjBlWEJsSWpvaVZXNXBiMjVTWlc1a1pYSmxjbk1pZlgwc0ltbGtJam9pTVRBMU9DSXNJblI1Y0dVaU9pSkRiMngxYlc1RVlYUmhVMjkxY21ObEluMHNleUpoZEhSeWFXSjFkR1Z6SWpwN0lteGhZbVZzSWpwN0luWmhiSFZsSWpvaVFYWmxjbUZuWlhOZlFtVnpkQ05tYVd4c2JXbHpiV0YwWTJnaWZTd2ljbVZ1WkdWeVpYSnpJanBiZXlKcFpDSTZJakV3TXpRaUxDSjBlWEJsSWpvaVIyeDVjR2hTWlc1a1pYSmxjaUo5WFgwc0ltbGtJam9pTVRBMU5TSXNJblI1Y0dVaU9pSk1aV2RsYm1SSmRHVnRJbjBzZXlKaGRIUnlhV0oxZEdWeklqcDdmU3dpYVdRaU9pSXhNVEV4SWl3aWRIbHdaU0k2SWxObGJHVmpkR2x2YmlKOUxIc2lZWFIwY21saWRYUmxjeUk2ZXlKc1lXSmxiQ0k2ZXlKMllXeDFaU0k2SWs1RlEwOUdVMTlHVmtOUFRWOVBRMFZCVGw5TlFWTlRRa0ZaWDBaUFVrVkRRVk5VTG01akkyWnBiR3h0YVhOdFlYUmphQ0o5TENKeVpXNWtaWEpsY25NaU9sdDdJbWxrSWpvaU1UQTRPU0lzSW5SNWNHVWlPaUpIYkhsd2FGSmxibVJsY21WeUluMWRmU3dpYVdRaU9pSXhNVEV6SWl3aWRIbHdaU0k2SWt4bFoyVnVaRWwwWlcwaWZTeDdJbUYwZEhKcFluVjBaWE1pT25zaWFYUmxiWE1pT2x0N0ltbGtJam9pTVRBMU5TSXNJblI1Y0dVaU9pSk1aV2RsYm1SSmRHVnRJbjBzZXlKcFpDSTZJakV3T0RNaUxDSjBlWEJsSWpvaVRHVm5aVzVrU1hSbGJTSjlMSHNpYVdRaU9pSXhNVEV6SWl3aWRIbHdaU0k2SWt4bFoyVnVaRWwwWlcwaWZTeDdJbWxrSWpvaU1URTBOU0lzSW5SNWNHVWlPaUpNWldkbGJtUkpkR1Z0SW4wc2V5SnBaQ0k2SWpFeE56a2lMQ0owZVhCbElqb2lUR1ZuWlc1a1NYUmxiU0o5WFN3aWNHeHZkQ0k2ZXlKcFpDSTZJakV3TURJaUxDSnpkV0owZVhCbElqb2lSbWxuZFhKbElpd2lkSGx3WlNJNklsQnNiM1FpZlgwc0ltbGtJam9pTVRBMU5DSXNJblI1Y0dVaU9pSk1aV2RsYm1RaWZTeDdJbUYwZEhKcFluVjBaWE1pT250OUxDSnBaQ0k2SWpFeE1UQWlMQ0owZVhCbElqb2lWVzVwYjI1U1pXNWtaWEpsY25NaWZTeDdJbUYwZEhKcFluVjBaWE1pT250OUxDSnBaQ0k2SWpFd05USWlMQ0owZVhCbElqb2lXV1ZoY25OVWFXTnJaWElpZlN4N0ltRjBkSEpwWW5WMFpYTWlPbnNpWW1Wc2IzY2lPbHQ3SW1sa0lqb2lNVEF4TWlJc0luUjVjR1VpT2lKRVlYUmxkR2x0WlVGNGFYTWlmVjBzSW14bFpuUWlPbHQ3SW1sa0lqb2lNVEF4TnlJc0luUjVjR1VpT2lKTWFXNWxZWEpCZUdsekluMWRMQ0p3Ykc5MFgyaGxhV2RvZENJNk1qVXdMQ0p3Ykc5MFgzZHBaSFJvSWpvM05UQXNJbkpsYm1SbGNtVnljeUk2VzNzaWFXUWlPaUl4TURFeUlpd2lkSGx3WlNJNklrUmhkR1YwYVcxbFFYaHBjeUo5TEhzaWFXUWlPaUl4TURFMklpd2lkSGx3WlNJNklrZHlhV1FpZlN4N0ltbGtJam9pTVRBeE55SXNJblI1Y0dVaU9pSk1hVzVsWVhKQmVHbHpJbjBzZXlKcFpDSTZJakV3TWpFaUxDSjBlWEJsSWpvaVIzSnBaQ0o5TEhzaWFXUWlPaUl4TURJMklpd2lkSGx3WlNJNklrSnZlRUZ1Ym05MFlYUnBiMjRpZlN4N0ltbGtJam9pTVRBMU5DSXNJblI1Y0dVaU9pSk1aV2RsYm1RaWZTeDdJbWxrSWpvaU1UQXpOQ0lzSW5SNWNHVWlPaUpIYkhsd2FGSmxibVJsY21WeUluMHNleUpwWkNJNklqRXdOakVpTENKMGVYQmxJam9pUjJ4NWNHaFNaVzVrWlhKbGNpSjlMSHNpYVdRaU9pSXhNRGc1SWl3aWRIbHdaU0k2SWtkc2VYQm9VbVZ1WkdWeVpYSWlmU3g3SW1sa0lqb2lNVEV4T1NJc0luUjVjR1VpT2lKSGJIbHdhRkpsYm1SbGNtVnlJbjBzZXlKcFpDSTZJakV4TlRFaUxDSjBlWEJsSWpvaVIyeDVjR2hTWlc1a1pYSmxjaUo5WFN3aWRHbDBiR1VpT25zaWFXUWlPaUl4TURBeElpd2lkSGx3WlNJNklsUnBkR3hsSW4wc0luUnZiMnhpWVhJaU9uc2lhV1FpT2lJeE1ESTFJaXdpZEhsd1pTSTZJbFJ2YjJ4aVlYSWlmU3dpZEc5dmJHSmhjbDlzYjJOaGRHbHZiaUk2SW1GaWIzWmxJaXdpZUY5eVlXNW5aU0k2ZXlKcFpDSTZJakV3TURRaUxDSjBlWEJsSWpvaVJHRjBZVkpoYm1kbE1XUWlmU3dpZUY5elkyRnNaU0k2ZXlKcFpDSTZJakV3TURnaUxDSjBlWEJsSWpvaVRHbHVaV0Z5VTJOaGJHVWlmU3dpZVY5eVlXNW5aU0k2ZXlKcFpDSTZJakV3TURZaUxDSjBlWEJsSWpvaVJHRjBZVkpoYm1kbE1XUWlmU3dpZVY5elkyRnNaU0k2ZXlKcFpDSTZJakV3TVRBaUxDSjBlWEJsSWpvaVRHbHVaV0Z5VTJOaGJHVWlmWDBzSW1sa0lqb2lNVEF3TWlJc0luTjFZblI1Y0dVaU9pSkdhV2QxY21VaUxDSjBlWEJsSWpvaVVHeHZkQ0o5TEhzaVlYUjBjbWxpZFhSbGN5STZleUpqWVd4c1ltRmpheUk2Ym5Wc2JIMHNJbWxrSWpvaU1UQXdOaUlzSW5SNWNHVWlPaUpFWVhSaFVtRnVaMlV4WkNKOUxIc2lZWFIwY21saWRYUmxjeUk2ZXlKemIzVnlZMlVpT25zaWFXUWlPaUl4TURnMklpd2lkSGx3WlNJNklrTnZiSFZ0YmtSaGRHRlRiM1Z5WTJVaWZYMHNJbWxrSWpvaU1UQTVNQ0lzSW5SNWNHVWlPaUpEUkZOV2FXVjNJbjBzZXlKaGRIUnlhV0oxZEdWeklqcDdJbTF2Ym5Sb2N5STZXekFzTmwxOUxDSnBaQ0k2SWpFd05URWlMQ0owZVhCbElqb2lUVzl1ZEdoelZHbGphMlZ5SW4wc2V5SmhkSFJ5YVdKMWRHVnpJanA3SW1OaGJHeGlZV05ySWpwdWRXeHNMQ0prWVhSaElqcDdJbmdpT25zaVgxOXVaR0Z5Y21GNVgxOGlPaUpCUVVSbk5FazJTMlJyU1VGQlRXaFFhMjl3TWxGblFVRnpURFpXYVc1YVEwRkJRMmRTZFVkTFpHdEpRVUZKYVRFMVNYQXlVV2RCUVdORFZHOXBibHBEUVVGQ1ozSkVUMHhrYTBsQlFVVm5ZazQwZERKUlowRkJUVWx2Tm1reldrTkJRVUZuUlc5aFRHUnJTVUZCUVdsQ2FWbDBNbEZuUVVFNFR5dE5hVE5hUTBGQlJHZGtPV2xNWkd0SlFVRk5hbTB5TkhReVVXZEJRWE5HV0dacE0xcERRVUZEWnpOVGNVMWthMGxCUVVsb1RVeHZlREpSWjBGQlkweHplR3BJV2tNaUxDSmtkSGx3WlNJNkltWnNiMkYwTmpRaUxDSnphR0Z3WlNJNld6RTRYWDBzSW5raU9uc2lYMTl1WkdGeWNtRjVYMThpT2lKQlFVRkJORWhXZUVSclFVRkJRVUpuWTFoQlQxRkJRVUZCVDBKelluYzFRVUZCUVVFMFFXeGFSR3RCUVVGQlJFRnpWVmxQVVVGQlFVRk5RbHBPUVRWQlFVRkJRVzlOWldkRVJVRkJRVUZEWjNwdmMwMVJRVUZCUVUxRVZtUm5lRUZCUVVGQmQwaERjRU5yUVVGQlFVUkJNM0IzUzFGQlFVRkJUVUpOYTBGd1FVRkJRVUZCVFVJM1ExVkJRVUZCUVdkR1NqQktVVUZCUVVGSFFtOTJaMnhCUVVGQlFWRkxVMkpFUlVGQlFVRkNRWEJLYzAxUlFVRkJRVVZEYTIxM2VFRWlMQ0prZEhsd1pTSTZJbVpzYjJGME5qUWlMQ0p6YUdGd1pTSTZXekU0WFgxOUxDSnpaV3hsWTNSbFpDSTZleUpwWkNJNklqRXdPREVpTENKMGVYQmxJam9pVTJWc1pXTjBhVzl1SW4wc0luTmxiR1ZqZEdsdmJsOXdiMnhwWTNraU9uc2lhV1FpT2lJeE1EZ3dJaXdpZEhsd1pTSTZJbFZ1YVc5dVVtVnVaR1Z5WlhKekluMTlMQ0pwWkNJNklqRXdNekVpTENKMGVYQmxJam9pUTI5c2RXMXVSR0YwWVZOdmRYSmpaU0o5TEhzaVlYUjBjbWxpZFhSbGN5STZleUp0YjI1MGFITWlPbHN3TERRc09GMTlMQ0pwWkNJNklqRXdOVEFpTENKMGVYQmxJam9pVFc5dWRHaHpWR2xqYTJWeUluMHNleUpoZEhSeWFXSjFkR1Z6SWpwN2ZTd2lhV1FpT2lJeE1ESTBJaXdpZEhsd1pTSTZJbEpsYzJWMFZHOXZiQ0o5TEhzaVlYUjBjbWxpZFhSbGN5STZleUp0YjI1MGFITWlPbHN3TERJc05DdzJMRGdzTVRCZGZTd2lhV1FpT2lJeE1EUTVJaXdpZEhsd1pTSTZJazF2Ym5Sb2MxUnBZMnRsY2lKOUxIc2lZWFIwY21saWRYUmxjeUk2ZTMwc0ltbGtJam9pTVRJd09DSXNJblI1Y0dVaU9pSlRaV3hsWTNScGIyNGlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2liblZ0WDIxcGJtOXlYM1JwWTJ0eklqbzFMQ0owYVdOclpYSnpJanBiZXlKcFpDSTZJakV3TkRFaUxDSjBlWEJsSWpvaVFXUmhjSFJwZG1WVWFXTnJaWElpZlN4N0ltbGtJam9pTVRBME1pSXNJblI1Y0dVaU9pSkJaR0Z3ZEdsMlpWUnBZMnRsY2lKOUxIc2lhV1FpT2lJeE1EUXpJaXdpZEhsd1pTSTZJa0ZrWVhCMGFYWmxWR2xqYTJWeUluMHNleUpwWkNJNklqRXdORFFpTENKMGVYQmxJam9pUkdGNWMxUnBZMnRsY2lKOUxIc2lhV1FpT2lJeE1EUTFJaXdpZEhsd1pTSTZJa1JoZVhOVWFXTnJaWElpZlN4N0ltbGtJam9pTVRBME5pSXNJblI1Y0dVaU9pSkVZWGx6VkdsamEyVnlJbjBzZXlKcFpDSTZJakV3TkRjaUxDSjBlWEJsSWpvaVJHRjVjMVJwWTJ0bGNpSjlMSHNpYVdRaU9pSXhNRFE0SWl3aWRIbHdaU0k2SWsxdmJuUm9jMVJwWTJ0bGNpSjlMSHNpYVdRaU9pSXhNRFE1SWl3aWRIbHdaU0k2SWsxdmJuUm9jMVJwWTJ0bGNpSjlMSHNpYVdRaU9pSXhNRFV3SWl3aWRIbHdaU0k2SWsxdmJuUm9jMVJwWTJ0bGNpSjlMSHNpYVdRaU9pSXhNRFV4SWl3aWRIbHdaU0k2SWsxdmJuUm9jMVJwWTJ0bGNpSjlMSHNpYVdRaU9pSXhNRFV5SWl3aWRIbHdaU0k2SWxsbFlYSnpWR2xqYTJWeUluMWRmU3dpYVdRaU9pSXhNREV6SWl3aWRIbHdaU0k2SWtSaGRHVjBhVzFsVkdsamEyVnlJbjBzZXlKaGRIUnlhV0oxZEdWeklqcDdJbVJoZEdGZmMyOTFjbU5sSWpwN0ltbGtJam9pTVRBMU9DSXNJblI1Y0dVaU9pSkRiMngxYlc1RVlYUmhVMjkxY21ObEluMHNJbWRzZVhCb0lqcDdJbWxrSWpvaU1UQTFPU0lzSW5SNWNHVWlPaUpNYVc1bEluMHNJbWh2ZG1WeVgyZHNlWEJvSWpwdWRXeHNMQ0p0ZFhSbFpGOW5iSGx3YUNJNmJuVnNiQ3dpYm05dWMyVnNaV04wYVc5dVgyZHNlWEJvSWpwN0ltbGtJam9pTVRBMk1DSXNJblI1Y0dVaU9pSk1hVzVsSW4wc0luTmxiR1ZqZEdsdmJsOW5iSGx3YUNJNmJuVnNiQ3dpZG1sbGR5STZleUpwWkNJNklqRXdOaklpTENKMGVYQmxJam9pUTBSVFZtbGxkeUo5ZlN3aWFXUWlPaUl4TURZeElpd2lkSGx3WlNJNklrZHNlWEJvVW1WdVpHVnlaWElpZlN4N0ltRjBkSEpwWW5WMFpYTWlPbnNpWm05eWJXRjBkR1Z5SWpwN0ltbGtJam9pTVRBek9TSXNJblI1Y0dVaU9pSkVZWFJsZEdsdFpWUnBZMnRHYjNKdFlYUjBaWElpZlN3aWNHeHZkQ0k2ZXlKcFpDSTZJakV3TURJaUxDSnpkV0owZVhCbElqb2lSbWxuZFhKbElpd2lkSGx3WlNJNklsQnNiM1FpZlN3aWRHbGphMlZ5SWpwN0ltbGtJam9pTVRBeE15SXNJblI1Y0dVaU9pSkVZWFJsZEdsdFpWUnBZMnRsY2lKOWZTd2lhV1FpT2lJeE1ERXlJaXdpZEhsd1pTSTZJa1JoZEdWMGFXMWxRWGhwY3lKOUxIc2lZWFIwY21saWRYUmxjeUk2ZXlKd2JHOTBJanB1ZFd4c0xDSjBaWGgwSWpvaU5EUXdNVE1pZlN3aWFXUWlPaUl4TURBeElpd2lkSGx3WlNJNklsUnBkR3hsSW4wc2V5SmhkSFJ5YVdKMWRHVnpJanA3SW0xdmJuUm9jeUk2V3pBc01Td3lMRE1zTkN3MUxEWXNOeXc0TERrc01UQXNNVEZkZlN3aWFXUWlPaUl4TURRNElpd2lkSGx3WlNJNklrMXZiblJvYzFScFkydGxjaUo5TEhzaVlYUjBjbWxpZFhSbGN5STZleUpzYVc1bFgyRnNjR2hoSWpvd0xqWTFMQ0pzYVc1bFgyTmhjQ0k2SW5KdmRXNWtJaXdpYkdsdVpWOWpiMnh2Y2lJNklpTXhaamMzWWpRaUxDSnNhVzVsWDJwdmFXNGlPaUp5YjNWdVpDSXNJbXhwYm1WZmQybGtkR2dpT2pVc0luZ2lPbnNpWm1sbGJHUWlPaUo0SW4wc0lua2lPbnNpWm1sbGJHUWlPaUo1SW4xOUxDSnBaQ0k2SWpFd016SWlMQ0owZVhCbElqb2lUR2x1WlNKOUxIc2lZWFIwY21saWRYUmxjeUk2ZXlKa1lYbHpJanBiTVN3eE5WMTlMQ0pwWkNJNklqRXdORGNpTENKMGVYQmxJam9pUkdGNWMxUnBZMnRsY2lKOUxIc2lZWFIwY21saWRYUmxjeUk2ZXlKa1lYbHpJanBiTVN3NExERTFMREl5WFgwc0ltbGtJam9pTVRBME5pSXNJblI1Y0dVaU9pSkVZWGx6VkdsamEyVnlJbjBzZXlKaGRIUnlhV0oxZEdWeklqcDdJbU5oYkd4aVlXTnJJanB1ZFd4c0xDSnlaVzVrWlhKbGNuTWlPbHQ3SW1sa0lqb2lNVEE0T1NJc0luUjVjR1VpT2lKSGJIbHdhRkpsYm1SbGNtVnlJbjFkTENKMGIyOXNkR2x3Y3lJNlcxc2lUbUZ0WlNJc0lrNUZRMDlHVTE5R1ZrTlBUVjlQUTBWQlRsOU5RVk5UUWtGWlgwWlBVa1ZEUVZOVUxtNWpJMlpwYkd4dGFYTnRZWFJqYUNKZExGc2lRbWxoY3lJc0lpMHdMak16SWwwc1d5SlRhMmxzYkNJc0lqQXVNelVpWFYxOUxDSnBaQ0k2SWpFeE1UUWlMQ0owZVhCbElqb2lTRzkyWlhKVWIyOXNJbjBzZXlKaGRIUnlhV0oxZEdWeklqcDdJbVJoZVhNaU9sc3hMRFFzTnl3eE1Dd3hNeXd4Tml3eE9Td3lNaXd5TlN3eU9GMTlMQ0pwWkNJNklqRXdORFVpTENKMGVYQmxJam9pUkdGNWMxUnBZMnRsY2lKOUxIc2lZWFIwY21saWRYUmxjeUk2ZXlKa1lYUmhYM052ZFhKalpTSTZleUpwWkNJNklqRXdPRFlpTENKMGVYQmxJam9pUTI5c2RXMXVSR0YwWVZOdmRYSmpaU0o5TENKbmJIbHdhQ0k2ZXlKcFpDSTZJakV3T0RjaUxDSjBlWEJsSWpvaVRHbHVaU0o5TENKb2IzWmxjbDluYkhsd2FDSTZiblZzYkN3aWJYVjBaV1JmWjJ4NWNHZ2lPbTUxYkd3c0ltNXZibk5sYkdWamRHbHZibDluYkhsd2FDSTZleUpwWkNJNklqRXdPRGdpTENKMGVYQmxJam9pVEdsdVpTSjlMQ0p6Wld4bFkzUnBiMjVmWjJ4NWNHZ2lPbTUxYkd3c0luWnBaWGNpT25zaWFXUWlPaUl4TURrd0lpd2lkSGx3WlNJNklrTkVVMVpwWlhjaWZYMHNJbWxrSWpvaU1UQTRPU0lzSW5SNWNHVWlPaUpIYkhsd2FGSmxibVJsY21WeUluMHNleUpoZEhSeWFXSjFkR1Z6SWpwN0ltUmhlWE1pT2xzeExESXNNeXcwTERVc05pdzNMRGdzT1N3eE1Dd3hNU3d4TWl3eE15d3hOQ3d4TlN3eE5pd3hOeXd4T0N3eE9Td3lNQ3d5TVN3eU1pd3lNeXd5TkN3eU5Td3lOaXd5Tnl3eU9Dd3lPU3d6TUN3ek1WMTlMQ0pwWkNJNklqRXdORFFpTENKMGVYQmxJam9pUkdGNWMxUnBZMnRsY2lKOUxIc2lZWFIwY21saWRYUmxjeUk2ZXlKa1lYUmhYM052ZFhKalpTSTZleUpwWkNJNklqRXhNVFlpTENKMGVYQmxJam9pUTI5c2RXMXVSR0YwWVZOdmRYSmpaU0o5TENKbmJIbHdhQ0k2ZXlKcFpDSTZJakV4TVRjaUxDSjBlWEJsSWpvaVRHbHVaU0o5TENKb2IzWmxjbDluYkhsd2FDSTZiblZzYkN3aWJYVjBaV1JmWjJ4NWNHZ2lPbTUxYkd3c0ltNXZibk5sYkdWamRHbHZibDluYkhsd2FDSTZleUpwWkNJNklqRXhNVGdpTENKMGVYQmxJam9pVEdsdVpTSjlMQ0p6Wld4bFkzUnBiMjVmWjJ4NWNHZ2lPbTUxYkd3c0luWnBaWGNpT25zaWFXUWlPaUl4TVRJd0lpd2lkSGx3WlNJNklrTkVVMVpwWlhjaWZYMHNJbWxrSWpvaU1URXhPU0lzSW5SNWNHVWlPaUpIYkhsd2FGSmxibVJsY21WeUluMHNleUpoZEhSeWFXSjFkR1Z6SWpwN0luTnZkWEpqWlNJNmV5SnBaQ0k2SWpFeE1UWWlMQ0owZVhCbElqb2lRMjlzZFcxdVJHRjBZVk52ZFhKalpTSjlmU3dpYVdRaU9pSXhNVEl3SWl3aWRIbHdaU0k2SWtORVUxWnBaWGNpZlN4N0ltRjBkSEpwWW5WMFpYTWlPbnQ5TENKcFpDSTZJakV3T0RFaUxDSjBlWEJsSWpvaVUyVnNaV04wYVc5dUluMHNleUpoZEhSeWFXSjFkR1Z6SWpwN0ltRmpkR2wyWlY5a2NtRm5Jam9pWVhWMGJ5SXNJbUZqZEdsMlpWOXBibk53WldOMElqb2lZWFYwYnlJc0ltRmpkR2wyWlY5dGRXeDBhU0k2Ym5Wc2JDd2lZV04wYVhabFgzTmpjbTlzYkNJNkltRjFkRzhpTENKaFkzUnBkbVZmZEdGd0lqb2lZWFYwYnlJc0luUnZiMnh6SWpwYmV5SnBaQ0k2SWpFd01qSWlMQ0owZVhCbElqb2lVR0Z1Vkc5dmJDSjlMSHNpYVdRaU9pSXhNREl6SWl3aWRIbHdaU0k2SWtKdmVGcHZiMjFVYjI5c0luMHNleUpwWkNJNklqRXdNalFpTENKMGVYQmxJam9pVW1WelpYUlViMjlzSW4wc2V5SnBaQ0k2SWpFd05UWWlMQ0owZVhCbElqb2lTRzkyWlhKVWIyOXNJbjBzZXlKcFpDSTZJakV3T0RRaUxDSjBlWEJsSWpvaVNHOTJaWEpVYjI5c0luMHNleUpwWkNJNklqRXhNVFFpTENKMGVYQmxJam9pU0c5MlpYSlViMjlzSW4wc2V5SnBaQ0k2SWpFeE5EWWlMQ0owZVhCbElqb2lTRzkyWlhKVWIyOXNJbjBzZXlKcFpDSTZJakV4T0RBaUxDSjBlWEJsSWpvaVNHOTJaWEpVYjI5c0luMWRmU3dpYVdRaU9pSXhNREkxSWl3aWRIbHdaU0k2SWxSdmIyeGlZWElpZlN4N0ltRjBkSEpwWW5WMFpYTWlPbnNpWW1GelpTSTZNalFzSW0xaGJuUnBjM05oY3lJNld6RXNNaXcwTERZc09Dd3hNbDBzSW0xaGVGOXBiblJsY25aaGJDSTZORE15TURBd01EQXVNQ3dpYldsdVgybHVkR1Z5ZG1Gc0lqb3pOakF3TURBd0xqQXNJbTUxYlY5dGFXNXZjbDkwYVdOcmN5STZNSDBzSW1sa0lqb2lNVEEwTXlJc0luUjVjR1VpT2lKQlpHRndkR2wyWlZScFkydGxjaUo5TEhzaVlYUjBjbWxpZFhSbGN5STZlMzBzSW1sa0lqb2lNVEUwTWlJc0luUjVjR1VpT2lKVmJtbHZibEpsYm1SbGNtVnljeUo5TEhzaVlYUjBjbWxpZFhSbGN5STZlMzBzSW1sa0lqb2lNVEE0TUNJc0luUjVjR1VpT2lKVmJtbHZibEpsYm1SbGNtVnljeUo5TEhzaVlYUjBjbWxpZFhSbGN5STZleUpqWVd4c1ltRmpheUk2Ym5Wc2JDd2ljbVZ1WkdWeVpYSnpJanBiZXlKcFpDSTZJakV3TXpRaUxDSjBlWEJsSWpvaVIyeDVjR2hTWlc1a1pYSmxjaUo5WFN3aWRHOXZiSFJwY0hNaU9sdGJJazVoYldVaUxDSkJkbVZ5WVdkbGMxOUNaWE4wSTJacGJHeHRhWE50WVhSamFDSmRMRnNpUW1saGN5SXNJaTB3TGpVMUlsMHNXeUpUYTJsc2JDSXNJakF1TXpJaVhWMTlMQ0pwWkNJNklqRXdOVFlpTENKMGVYQmxJam9pU0c5MlpYSlViMjlzSW4wc2V5SmhkSFJ5YVdKMWRHVnpJanA3SW1KaGMyVWlPall3TENKdFlXNTBhWE56WVhNaU9sc3hMRElzTlN3eE1Dd3hOU3d5TUN3ek1GMHNJbTFoZUY5cGJuUmxjblpoYkNJNk1UZ3dNREF3TUM0d0xDSnRhVzVmYVc1MFpYSjJZV3dpT2pFd01EQXVNQ3dpYm5WdFgyMXBibTl5WDNScFkydHpJam93ZlN3aWFXUWlPaUl4TURReUlpd2lkSGx3WlNJNklrRmtZWEIwYVhabFZHbGphMlZ5SW4wc2V5SmhkSFJ5YVdKMWRHVnpJanA3SW05MlpYSnNZWGtpT25zaWFXUWlPaUl4TURJMklpd2lkSGx3WlNJNklrSnZlRUZ1Ym05MFlYUnBiMjRpZlgwc0ltbGtJam9pTVRBeU15SXNJblI1Y0dVaU9pSkNiM2hhYjI5dFZHOXZiQ0o5TEhzaVlYUjBjbWxpZFhSbGN5STZlMzBzSW1sa0lqb2lNVEUwTXlJc0luUjVjR1VpT2lKVFpXeGxZM1JwYjI0aWZTeDdJbUYwZEhKcFluVjBaWE1pT25zaWJHbHVaVjloYkhCb1lTSTZNQzQyTlN3aWJHbHVaVjlqWVhBaU9pSnliM1Z1WkNJc0lteHBibVZmWTI5c2IzSWlPaUlqWm1ZM1pqQmxJaXdpYkdsdVpWOXFiMmx1SWpvaWNtOTFibVFpTENKc2FXNWxYM2RwWkhSb0lqbzFMQ0o0SWpwN0ltWnBaV3hrSWpvaWVDSjlMQ0o1SWpwN0ltWnBaV3hrSWpvaWVTSjlmU3dpYVdRaU9pSXhNRGczSWl3aWRIbHdaU0k2SWt4cGJtVWlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2liR0ZpWld3aU9uc2lkbUZzZFdVaU9pSlBZbk5sY25aaGRHbHZibk1pZlN3aWNtVnVaR1Z5WlhKeklqcGJleUpwWkNJNklqRXhNVGtpTENKMGVYQmxJam9pUjJ4NWNHaFNaVzVrWlhKbGNpSjlYWDBzSW1sa0lqb2lNVEUwTlNJc0luUjVjR1VpT2lKTVpXZGxibVJKZEdWdEluMHNleUpoZEhSeWFXSjFkR1Z6SWpwN0lteHBibVZmWVd4d2FHRWlPakF1TVN3aWJHbHVaVjlqWVhBaU9pSnliM1Z1WkNJc0lteHBibVZmWTI5c2IzSWlPaUlqTVdZM04ySTBJaXdpYkdsdVpWOXFiMmx1SWpvaWNtOTFibVFpTENKc2FXNWxYM2RwWkhSb0lqbzFMQ0o0SWpwN0ltWnBaV3hrSWpvaWVDSjlMQ0o1SWpwN0ltWnBaV3hrSWpvaWVTSjlmU3dpYVdRaU9pSXhNRGc0SWl3aWRIbHdaU0k2SWt4cGJtVWlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2ljR3h2ZENJNmV5SnBaQ0k2SWpFd01ESWlMQ0p6ZFdKMGVYQmxJam9pUm1sbmRYSmxJaXdpZEhsd1pTSTZJbEJzYjNRaWZTd2lkR2xqYTJWeUlqcDdJbWxrSWpvaU1UQXhNeUlzSW5SNWNHVWlPaUpFWVhSbGRHbHRaVlJwWTJ0bGNpSjlmU3dpYVdRaU9pSXhNREUySWl3aWRIbHdaU0k2SWtkeWFXUWlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2laR2x0Wlc1emFXOXVJam94TENKd2JHOTBJanA3SW1sa0lqb2lNVEF3TWlJc0luTjFZblI1Y0dVaU9pSkdhV2QxY21VaUxDSjBlWEJsSWpvaVVHeHZkQ0o5TENKMGFXTnJaWElpT25zaWFXUWlPaUl4TURFNElpd2lkSGx3WlNJNklrSmhjMmxqVkdsamEyVnlJbjE5TENKcFpDSTZJakV3TWpFaUxDSjBlWEJsSWpvaVIzSnBaQ0o5TEhzaVlYUjBjbWxpZFhSbGN5STZleUpqWVd4c1ltRmpheUk2Ym5Wc2JDd2ljbVZ1WkdWeVpYSnpJanBiZXlKcFpDSTZJakV4TVRraUxDSjBlWEJsSWpvaVIyeDVjR2hTWlc1a1pYSmxjaUo5WFN3aWRHOXZiSFJwY0hNaU9sdGJJazVoYldVaUxDSlBZbk5sY25aaGRHbHZibk1pWFN4YklrSnBZWE1pTENKT1FTSmRMRnNpVTJ0cGJHd2lMQ0pPUVNKZFhYMHNJbWxrSWpvaU1URTBOaUlzSW5SNWNHVWlPaUpJYjNabGNsUnZiMndpZlN4N0ltRjBkSEpwWW5WMFpYTWlPbnNpYldGdWRHbHpjMkZ6SWpwYk1Td3lMRFZkTENKdFlYaGZhVzUwWlhKMllXd2lPalV3TUM0d0xDSnVkVzFmYldsdWIzSmZkR2xqYTNNaU9qQjlMQ0pwWkNJNklqRXdOREVpTENKMGVYQmxJam9pUVdSaGNIUnBkbVZVYVdOclpYSWlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2liR2x1WlY5aGJIQm9ZU0k2TUM0eExDSnNhVzVsWDJOaGNDSTZJbkp2ZFc1a0lpd2liR2x1WlY5amIyeHZjaUk2SWlNeFpqYzNZalFpTENKc2FXNWxYMnB2YVc0aU9pSnliM1Z1WkNJc0lteHBibVZmZDJsa2RHZ2lPalVzSW5naU9uc2labWxsYkdRaU9pSjRJbjBzSW5raU9uc2labWxsYkdRaU9pSjVJbjE5TENKcFpDSTZJakV3TmpBaUxDSjBlWEJsSWpvaVRHbHVaU0o5TEhzaVlYUjBjbWxpZFhSbGN5STZlMzBzSW1sa0lqb2lNVEF6T1NJc0luUjVjR1VpT2lKRVlYUmxkR2x0WlZScFkydEdiM0p0WVhSMFpYSWlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2lZMkZzYkdKaFkyc2lPbTUxYkd3c0ltUmhkR0VpT25zaWVDSTZleUpmWDI1a1lYSnlZWGxmWHlJNklrRkJRVUZ5YlZkTFpHdEpRVUZQWjJOaFdYQXlVV2RCUVRCSmRITnBibHBEUVVGRE5DdHRLMHRrYTBsQlFVdENjR00wY0RKUlowRkJhVTVvTW1sdVdrTkJRVUozVWpOeFMyUnJTVUZCUm1reVpsbHdNbEZuUVVGUlExZENhVzVhUTBGQlFXOXNTVk5MWkd0SlFVRkNRVVJwU1hBeVVXZEJRU3RJUjB4cGJscERRVUZFWnpSSk5rdGthMGxCUVUxb1VHdHZjREpSWjBGQmMwdzJWbWx1V2tOQlFVTlpURnB0UzJSclNVRkJTVU5qYmtsd01sRm5RVUZoUVhWbmFXNWFRMEZCUWxGbGNVOUxaR3RKUVVGRWFuQndiM0F5VVdkQlFVbEdhWEZwYmxwRFFVRkJTWGcyTWt0a2EwbEJRVkJCTVhOWmNESlJaMEZCTWt0VE1HbHVXa05CUVVSQlJUZHBTMlJyU1VGQlMybERkVFJ3TWxGblFVRnJVRWNyYVc1YVEwRkJRalJaVFV0TFpHdEpRVUZIUkZCNFdYQXlVV2RCUVZORU4wcHBibHBEUVVGQmQzSmplVXRrYTBsQlFVSm5ZekJKY0RKUlowRkJRVWwyVkdsdVdrTkJRVVJ2SzJSaFMyUnJTVUZCVGtKdk1tOXdNbEZuUVVGMVRtWmthVzVhUTBGQlEyZFNkVWRMWkd0SlFVRkphVEUxU1hBeVVXZEJRV05EVkc5cGJscERRVUZDV1dzcmRVdGthMGxCUVVWQlF6YzBjREpSWjBGQlMwaEllV2x1V2tOQlFVRlJORkJYUzJSclNVRkJVR2hQSzFsd01sRm5RVUUwVERNNGFXNWFRMEZCUkVsTVFVTk1aR3RKUVVGTVEySkJOSFF5VVdkQlFXMUJiMGhwTTFwRFFVRkRRV1ZSY1V4a2EwbEJRVWRxYjBSWmRESlJaMEZCVlVaalVta3pXa05CUVVFMGVHaFRUR1JyU1VGQlEwRXhSMGwwTWxGblFVRkRTMUZpYVROYVEwRkJSSGRGYUN0TVpHdEpRVUZPYVVKSmIzUXlVV2RCUVhkUVFXeHBNMXBEUVVGRGIxaDViVXhrYTBsQlFVcEVUMHhKZERKUlowRkJaVVF3ZDJreldrTkJRVUpuY2tSUFRHUnJTVUZCUldkaVRqUjBNbEZuUVVGTlNXODJhVE5hUTBGQlFWa3JWREpNWkd0SlFVRkJRbTlSV1hReVVXZEJRVFpPV2tWcE0xcERRVUZFVVZKVmFVeGthMGxCUVV4cE1GTTBkREpSWjBGQmIwTk9VR2t6V2tOQlFVTkphMnhMVEdSclNVRkJTRUZDVm05ME1sRm5RVUZYU0VKYWFUTmFRMEZCUWtFek1YbE1aR3RKUVVGRGFFOVpTWFF5VVdkQlFVVk1NV3BwTTFwRFFVRkVORXN5WlV4a2EwbEJRVTlEWVdGdmRESlJaMEZCZVVGc2RXa3pXa05CUVVOM1pVaEhUR1JyU1VGQlNtcHVaRWwwTWxGblFVRm5SbG8wYVROYVEwRkJRbTk0V0hWTVpHdEpRVUZHUVRCbU5IUXlVV2RCUVU5TFQwTnBNMXBEUVVGQlowVnZZVXhrYTBsQlFVRnBRbWxaZERKUlowRkJPRThyVFdreldrTkJRVVJaV0hCRFRHUnJTVUZCVFVST2F6UjBNbEZuUVVGeFJIbFlhVE5hUTBGQlExRnhOWEZNWkd0SlFVRklaMkZ1YjNReVVXZEJRVmxKYldocE0xcERRVUZDU1N0TFUweGthMGxCUVVSQ2JuRkpkREpSWjBGQlIwNWhjbWt6V2tOQlFVRkJVbUVyVEdSclNVRkJUMmw2YzI5ME1sRm5RVUV3UTBzeWFUTmFRMEZCUXpSclltMU1aR3RKUVVGTFFVRjJXWFF5VVdkQlFXbEhMMEZwTTFwRFFVRkNkek56VDB4a2EwbEJRVVpvVG5nMGRESlJaMEZCVVV4NlMya3pXa05CUVVGdlN6ZzJUR1JyU1VGQlFrTmhNRmwwTWxGblFVRXJRV3BXYVROYVEwRkJSR2RrT1dsTVpHdEpRVUZOYW0weU5IUXlVV2RCUVhOR1dHWnBNMXBEUVVGRFdYaFBTMHhrYTBsQlFVbEJlalZ2ZERKUlowRkJZVXRNY0dreldrTkJRVUpSUldVeVRHUnJTVUZCUkdsQk9FbDBNbEZuUVVGSlR5OTZhVE5hUTBGQlFVbFlkbVZNWkd0SlFVRlFSRTByYjNReVVXZEJRVEpFZGl0cE0xcERRVUZFUVhGblIwMWthMGxCUVV0bldrSlplREpSWjBGQmEwbG5TV3BJV2tOQlFVSTBPWGQxVFdSclNVRkJSMEp0UkRSNE1sRm5RVUZUVGxWVGFraGFRMEZCUVhkU1FtRk5aR3RKUVVGQ2FYcEhXWGd5VVdkQlFVRkRTV1JxU0ZwRFFVRkViMnREUTAxa2EwbEJRVTVFTDBrMGVESlJaMEZCZFVjMGJtcElXa05CUVVObk0xTnhUV1JyU1VGQlNXaE5URzk0TWxGblFVRmpUSE40YWtoYVEwRkJRbGxMYWxkTlpHdEpRVUZGUTFwUFNYZ3lVV2RCUVV0Qlp6aHFTRnBEUVVGQlVXUjZLMDFrYTBsQlFWQnFiRkZ2ZURKUlowRkJORVpTUjJwSVdrTkJRVVJKZHpCdFRXUnJTVUZCVEVGNVZGbDRNbEZuUVVGdFMwWlJha2hhUTBGQlEwRkZSbE5OWkd0SlFVRkhhQzlXTkhneVVXZEJRVlZQTldGcVNGcERRVUZCTkZoV05rMWthMGxCUVVORVRWbFplREpSWjBGQlEwUjBiR3BJV2tOQlFVUjNjVmRwVFdSclNVRkJUbWRaWWtsNE1sRm5RVUYzU1dSMmFraGFRMEZCUTI4NWJrdE5aR3RKUVVGS1FteGtiM2d5VVdkQlFXVk9ValZxU0ZwRFFVRkNaMUV6TWsxa2EwbEJRVVZwZVdkSmVESlJaMEZCVFVOSFJXcElXa05CUVVGWmEwbGxUV1JyU1VGQlFVUXZhVzk0TWxGblFVRTJSekpQYWtoYVEwRkJSRkV6U2tkTlpHdEpRVUZNYUV4c1dYZ3lVV2RCUVc5TWNWbHFTRnBEUVVGRFNVdGFlVTFrYTBsQlFVaERXVzQwZURKUlowRkJWMEZsYW1wSVdrTWlMQ0prZEhsd1pTSTZJbVpzYjJGME5qUWlMQ0p6YUdGd1pTSTZXekUyT0YxOUxDSjVJanA3SWw5ZmJtUmhjbkpoZVY5Zklqb2lRVUZCUVZsUFRtMUZSVUZCUVVGRVFVNUdSVkZSUVVGQlFVVkRSMDk0UWtGQlFVRkJiMDVqYkVWRlFVRkJRVU5uTWxFMFVWRkJRVUZCUTBNek4zYzVRVUZCUVVGSlRIWkNSREJCUVVGQlFtZGtPRTFRVVVGQlFVRkxRWHA0VVRsQlFVRkJRVFJQTDBkRU1FRkJRVUZFUVRsek5GQlJRVUZCUVVsRU9URm5PVUZCUVVGQldVRlVaa1F3UVVGQlFVUkJRbVEwVUZGQlFVRkJSVUZJTTFFNVFVRkJRVUZ2UVdwalJEQkJRVUZCUVVFeU9FMVFVVUZCUVVGRlEzUnhkemxCUVVGQlFXOUlLMVJFTUVGQlFVRkVaMUJaUlZCUlFVRkJRVU5FT0dKbk9VRkJRVUZCV1V4d1kwUXdRVUZCUVVSQlIyMWpVRkZCUVVGQlJVSTNZMUU1UVVGQlFVRnZUblEzUkRCQlFVRkJSR2RCUjFsUVVVRkJRVUZCUVcxVlFUbEJRVUZCUVZGRmN6WkVNRUZCUVVGQlp6bEJTVkJSUVVGQlFVTkRaSGwzTlVGQlFVRkJRVVZoVlVSclFVRkJRVUpCVWxrd1QxRkJRVUZCUjBKRmFHYzFRVUZCUVVGdlJVNHZSR3RCUVVGQlFXZEhTMFZQVVVGQlFVRkxSSE4zWnpWQlFVRkJRVWxOU0d0RWEwRkJRVUZDWjJOMk9FOVJRVUZCUVVsQmFrZG5PVUZCUVVGQmQwNVJNRVF3UVVGQlFVUm5ibmRaVUZGQlFVRkJRVUp5TWtFMVFVRkJRVUZKUkdGeFJHdEJRVUZCUVVGRFdrRlBVVUZCUVVGTlJHSmtVVFZCUVVGQlFXOUxOV0pFYTBGQlFVRkJRWG8wVVU5UlFVRkJRVWRFZG5KUk5VRkJRVUZCZDBFdldFUnJRVUZCUVVKQmVXVlZUMUZCUVVGQlRVTkRPVUUxUVVGQlFVRlJSSGRFUkRCQlFVRkJRVUZPZERCUFVVRkJRVUZOUVhaMGR6VkJRVUZCUVdkRGJWSkVhMEZCUVVGRFFWaFlkMDlSUVVGQlFVbERVbHAzTlVGQlFVRkJaMDFXVTBSclFVRkJRVUZCU2tnNFQxRkJRVUZCU1VORGNYYzFRVUZCUVVGQlQwaFlSR3RCUVVGQlEyZG9RakJRVVVGQlFVRkhRVzlaZHpsQlFVRkJRVUZOZVc5RU1FRkJRVUZDWjJwbk1GRlJRVUZCUVU5RE1sSm9Ra0ZCUVVGQlVVNDVMMFZGUVVGQlFVSkJiVFJqVVZGQlFVRkJSVUpZYW5oQ1FVRkJRVUZSUWs5WVJVVkJRVUZCUTJkMEx6UlFVVUZCUVVGTlFrbDZkelZCUVVGQlFUUk9iV1pFVlVGQlFVRkVaM2hLTUU1UlFVRkJRVTFEZG0xM01VRkJRVUZCZDBweFdrUlZRVUZCUVVSbmN6SjNUbEZCUVVGQlQwUk5VSGN4UVVGQlFVRkJUMWxUUkZWQlFVRkJRVUV6SzFWTlVVRkJRVUZCUkZsMVFYaEJRVUZCUVVGT1IweEVSVUZCUVVGQlozZGhhMDFSUVVGQlFVTkRlSGgzZUVGQlFVRkJVVXRJYkVSRlFVRkJRVUpCUTBSQlRsRkJRVUZCUTBKMlpXY3hRVUZCUVVGSlRtSkZSRlZCUVVGQlFVRnVWVTFQVVVGQlFVRlBRbXAzWnpWQlFVRkJRWGREY0VKRU1FRkJRVUZFWjIxcE9GQlJRVUZCUVVGQlRFaG5PVUZCUVVGQlNVaHpUVVF3UVVGQlFVUm5UalEwVUZGQlFVRkJSMFEyUW5oQ1FVRkJRVUYzVG1oSlJVVkJRVUZCUTBGNlZVRlJVVUZCUVVGSFJFTlBRa0pCUVVGQlFVbE1ZM2RGUlVGQlFVRkNRV1JTU1ZGUlFVRkJRVWxDYlRaQk9VRkJRVUZCYjA5TGNrUXdRVUZCUVVOQlJHMUZVRkZCUVVGQlIwRTJSbWM1UVVGQlFVRlJSMkpNUkd0QlFVRkJRbWRUVGpSUFVVRkJRVUZKUVhFNFVUVkJRVUZCUVc5QmQwVkVNRUZCUVVGRFozZHVORkJSUVVGQlFVMUNOQ3RST1VGQlFVRkJXVUpqTmtWRlFVRkJRVU5uVGtoVlVWRkJRVUZCVFVKU2MwSkNRVUZCUVVGQlJ5OXlSVVZCUVVGQlFtY3daRFJSVVVGQlFVRkxRWG93YUVKQlFVRkJRVUZLWWtaRlJVRkJRVUZEUVVzMU5GRlJRVUZCUVU5RVFXUm9Ra0ZCUVVGQldVWmFVRVZGUVVGQlFVSkJNMFZ6VVZGQlFVRkJSVUpwVTBKQ1FVRkJRVUZKVDJoRlJVVkJRVUZCUkVGeVUzTlJVVUZCUVVGSlFucEZhRUpCUVVGQlFWRklUSGxFTUVGQlFVRkVaMnB5WTFCUlFVRkJRVXREY21aQk9VRkJRVUZCVVUxb1FrUXdRVUZCUVVOQmJVVkZVRkZCUVVGQlQwSnZVVkU1UVVGQlFVRkpSR3hDUkRCQlFVRkJRbWQwU2xsUVVVRkJRVUZKUVhZM1FUbEJRVUZCUVZsT1ZXZEZSVUZCUVVGQlFXcHFZMUZSUVVGQlFVbENSMVJvUWtGQlFVRkJTVkE1YTBWRlFVRkJRVVJCVEhwclVWRkJRVUZCU1VKblJGSkNRVUZCUVVGUlEweEVSREJCUVVGQlFXZGtTemhRVVVGQlFVRkJSRWR0ZHpsQlFVRkJRVFJDWlVsRU1FRkJRVUZDUVZGTVVWQlJRVUZCUVV0Q2J6UkJPVUZCUVVGQlowVm5SMFZGUVVGQlFVTm5SM2RqVVZGQlFVRkJTMFIxUW5oQ1FVRkJRVUYzVFVWSlJVVkJRVUZCUVdjMFpFRlFVVUZCUVVGTFFTdHJRVGxCUVVGQlFWRktlRkJFTUVGQlFVRkNRVzlxYjFCUlFVRkJRVWREYjBwUk9VRkJRVUZCV1VzMFVVUXdRVUZCUVVGbmMydEJVRkZCUVVGQlRVTXhZMEU1UVVGQlFVRm5URzFuUkRCQlFVRkJRMmRUTjI5UVVVRkJRVUZQUkdRd2R6bEJRVUZCUVVGSVJIUkVNRUZCUVVGQ1p6bGlORkJSUVVGQlFVMUNObXRCT1VGQlFVRkJTVUZDYVVRd1FVRkJRVUZuUVVkSlVGRkJRVUZCUTBGQldXYzVRU0lzSW1SMGVYQmxJam9pWm14dllYUTJOQ0lzSW5Ob1lYQmxJanBiTVRZNFhYMTlMQ0p6Wld4bFkzUmxaQ0k2ZXlKcFpDSTZJakV5TURnaUxDSjBlWEJsSWpvaVUyVnNaV04wYVc5dUluMHNJbk5sYkdWamRHbHZibDl3YjJ4cFkza2lPbnNpYVdRaU9pSXhNakEzSWl3aWRIbHdaU0k2SWxWdWFXOXVVbVZ1WkdWeVpYSnpJbjE5TENKcFpDSTZJakV4TkRnaUxDSjBlWEJsSWpvaVEyOXNkVzF1UkdGMFlWTnZkWEpqWlNKOUxIc2lZWFIwY21saWRYUmxjeUk2ZTMwc0ltbGtJam9pTVRBek55SXNJblI1Y0dVaU9pSkNZWE5wWTFScFkydEdiM0p0WVhSMFpYSWlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2liR2x1WlY5aGJIQm9ZU0k2TUM0Mk5Td2liR2x1WlY5allYQWlPaUp5YjNWdVpDSXNJbXhwYm1WZlkyOXNiM0lpT2lJalptWmlZamM0SWl3aWJHbHVaVjlxYjJsdUlqb2ljbTkxYm1RaUxDSnNhVzVsWDNkcFpIUm9Jam8xTENKNElqcDdJbVpwWld4a0lqb2llQ0o5TENKNUlqcDdJbVpwWld4a0lqb2llU0o5ZlN3aWFXUWlPaUl4TVRRNUlpd2lkSGx3WlNJNklreHBibVVpZlN4N0ltRjBkSEpwWW5WMFpYTWlPbnNpYzI5MWNtTmxJanA3SW1sa0lqb2lNVEExT0NJc0luUjVjR1VpT2lKRGIyeDFiVzVFWVhSaFUyOTFjbU5sSW4xOUxDSnBaQ0k2SWpFd05qSWlMQ0owZVhCbElqb2lRMFJUVm1sbGR5SjlMSHNpWVhSMGNtbGlkWFJsY3lJNmV5SmliM1IwYjIxZmRXNXBkSE1pT2lKelkzSmxaVzRpTENKbWFXeHNYMkZzY0doaElqcDdJblpoYkhWbElqb3dMalY5TENKbWFXeHNYMk52Ykc5eUlqcDdJblpoYkhWbElqb2liR2xuYUhSbmNtVjVJbjBzSW14bFpuUmZkVzVwZEhNaU9pSnpZM0psWlc0aUxDSnNaWFpsYkNJNkltOTJaWEpzWVhraUxDSnNhVzVsWDJGc2NHaGhJanA3SW5aaGJIVmxJam94TGpCOUxDSnNhVzVsWDJOdmJHOXlJanA3SW5aaGJIVmxJam9pWW14aFkyc2lmU3dpYkdsdVpWOWtZWE5vSWpwYk5DdzBYU3dpYkdsdVpWOTNhV1IwYUNJNmV5SjJZV3gxWlNJNk1uMHNJbkJzYjNRaU9tNTFiR3dzSW5KbGJtUmxjbDl0YjJSbElqb2lZM056SWl3aWNtbG5hSFJmZFc1cGRITWlPaUp6WTNKbFpXNGlMQ0owYjNCZmRXNXBkSE1pT2lKelkzSmxaVzRpZlN3aWFXUWlPaUl4TURJMklpd2lkSGx3WlNJNklrSnZlRUZ1Ym05MFlYUnBiMjRpZlN4N0ltRjBkSEpwWW5WMFpYTWlPbnNpYkdsdVpWOWhiSEJvWVNJNk1DNHhMQ0pzYVc1bFgyTmhjQ0k2SW5KdmRXNWtJaXdpYkdsdVpWOWpiMnh2Y2lJNklpTXhaamMzWWpRaUxDSnNhVzVsWDJwdmFXNGlPaUp5YjNWdVpDSXNJbXhwYm1WZmQybGtkR2dpT2pVc0luZ2lPbnNpWm1sbGJHUWlPaUo0SW4wc0lua2lPbnNpWm1sbGJHUWlPaUo1SW4xOUxDSnBaQ0k2SWpFeE5UQWlMQ0owZVhCbElqb2lUR2x1WlNKOUxIc2lZWFIwY21saWRYUmxjeUk2ZTMwc0ltbGtJam9pTVRBeE9DSXNJblI1Y0dVaU9pSkNZWE5wWTFScFkydGxjaUo5TEhzaVlYUjBjbWxpZFhSbGN5STZleUp6YjNWeVkyVWlPbnNpYVdRaU9pSXhNRE14SWl3aWRIbHdaU0k2SWtOdmJIVnRia1JoZEdGVGIzVnlZMlVpZlgwc0ltbGtJam9pTVRBek5TSXNJblI1Y0dVaU9pSkRSRk5XYVdWM0luMHNleUpoZEhSeWFXSjFkR1Z6SWpwN0ltUmhkR0ZmYzI5MWNtTmxJanA3SW1sa0lqb2lNVEUwT0NJc0luUjVjR1VpT2lKRGIyeDFiVzVFWVhSaFUyOTFjbU5sSW4wc0ltZHNlWEJvSWpwN0ltbGtJam9pTVRFME9TSXNJblI1Y0dVaU9pSk1hVzVsSW4wc0ltaHZkbVZ5WDJkc2VYQm9JanB1ZFd4c0xDSnRkWFJsWkY5bmJIbHdhQ0k2Ym5Wc2JDd2libTl1YzJWc1pXTjBhVzl1WDJkc2VYQm9JanA3SW1sa0lqb2lNVEUxTUNJc0luUjVjR1VpT2lKTWFXNWxJbjBzSW5ObGJHVmpkR2x2Ymw5bmJIbHdhQ0k2Ym5Wc2JDd2lkbWxsZHlJNmV5SnBaQ0k2SWpFeE5USWlMQ0owZVhCbElqb2lRMFJUVm1sbGR5SjlmU3dpYVdRaU9pSXhNVFV4SWl3aWRIbHdaU0k2SWtkc2VYQm9VbVZ1WkdWeVpYSWlmU3g3SW1GMGRISnBZblYwWlhNaU9udDlMQ0pwWkNJNklqRXhOellpTENKMGVYQmxJam9pVlc1cGIyNVNaVzVrWlhKbGNuTWlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2lZMkZzYkdKaFkyc2lPbTUxYkd3c0ltUmhkR0VpT25zaWVDSTZleUpmWDI1a1lYSnlZWGxmWHlJNklrRkJRMEZsVVhGTVpHdEpRVUZIYW05RVdYUXlVV2RCUVZWR1kxSnBNMXBEUVVGQk5IaG9VMHhrYTBsQlFVTkJNVWRKZERKUlowRkJRMHRSWW1reldrTkJRVVIzUldnclRHUnJTVUZCVG1sQ1NXOTBNbEZuUVVGM1VFRnNhVE5hUTBGQlEyOVllVzFNWkd0SlFVRktSRTlNU1hReVVXZEJRV1ZFTUhkcE0xcERRVUZDWjNKRVQweGthMGxCUVVWbllrNDBkREpSWjBGQlRVbHZObWt6V2tOQlFVRlpLMVF5VEdSclNVRkJRVUp2VVZsME1sRm5RVUUyVGxwRmFUTmFRMEZCUkZGU1ZXbE1aR3RKUVVGTWFUQlROSFF5VVdkQlFXOURUbEJwTTFwRFFVRkRTV3RzUzB4a2EwbEJRVWhCUWxadmRESlJaMEZCVjBoQ1dta3pXa05CUVVKQk16RjVUR1JyU1VGQlEyaFBXVWwwTWxGblFVRkZUREZxYVROYVEwRkJSRFJMTW1WTVpHdEpRVUZQUTJGaGIzUXlVV2RCUVhsQmJIVnBNMXBEUVVGRGQyVklSMHhrYTBsQlFVcHFibVJKZERKUlowRkJaMFphTkdreldrTkJRVUp2ZUZoMVRHUnJTVUZCUmtFd1pqUjBNbEZuUVVGUFMwOURhVE5hUTBGQlFXZEZiMkZNWkd0SlFVRkJhVUpwV1hReVVXZEJRVGhQSzAxcE0xcERRVUZFV1Zod1EweGthMGxCUVUxRVRtczBkREpSWjBGQmNVUjVXR2t6V2tOQlFVTlJjVFZ4VEdSclNVRkJTR2RoYm05ME1sRm5RVUZaU1cxb2FUTmFRMEZCUWtrclMxTk1aR3RKUVVGRVFtNXhTWFF5VVdkQlFVZE9ZWEpwTTFwRFFVRkJRVkpoSzB4a2EwbEJRVTlwZW5OdmRESlJaMEZCTUVOTE1ta3pXa05CUVVNMGEySnRUR1JyU1VGQlMwRkJkbGwwTWxGblFVRnBSeTlCYVROYVEwRkJRbmN6YzA5TVpHdEpRVUZHYUU1NE5IUXlVV2RCUVZGTWVrdHBNMXBEUVVGQmIwczROa3hrYTBsQlFVSkRZVEJaZERKUlowRkJLMEZxVm1reldrTkJRVVJuWkRscFRHUnJTVUZCVFdwdE1qUjBNbEZuUVVGelJsaG1hVE5hUTBGQlExbDRUMHRNWkd0SlFVRkpRWG8xYjNReVVXZEJRV0ZMVEhCcE0xcERRVUZDVVVWbE1reGthMGxCUVVScFFUaEpkREpSWjBGQlNVOHZlbWt6V2tOQlFVRkpXSFpsVEdSclNVRkJVRVJOSzI5ME1sRm5RVUV5UkhZcmFUTmFRMEZCUkVGeFowZE5aR3RKUVVGTFoxcENXWGd5VVdkQlFXdEpaMGxxU0ZwRFFVRkNORGwzZFUxa2EwbEJRVWRDYlVRMGVESlJaMEZCVTA1VlUycElXa05CUVVGM1VrSmhUV1JyU1VGQlFtbDZSMWw0TWxGblFVRkJRMGxrYWtoYVEwRkJSRzlyUTBOTlpHdEpRVUZPUkM5Sk5IZ3lVV2RCUVhWSE5HNXFTRnBEUVVGRFp6TlRjVTFrYTBsQlFVbG9UVXh2ZURKUlowRkJZMHh6ZUdwSVdrTkJRVUpaUzJwWFRXUnJTVUZCUlVOYVQwbDRNbEZuUVVGTFFXYzRha2hhUTBGQlFWRmtlaXROWkd0SlFVRlFhbXhSYjNneVVXZEJRVFJHVWtkcVNGcERRVUZFU1hjd2JVMWthMGxCUVV4QmVWUlplREpSWjBGQmJVdEdVV3BJV2tOQlFVTkJSVVpUVFdSclNVRkJSMmd2VmpSNE1sRm5RVUZWVHpWaGFraGFRMEZCUVRSWVZqWk5aR3RKUVVGRFJFMVpXWGd5VVdkQlFVTkVkR3hxU0ZwRFFVRkVkM0ZYYVUxa2EwbEJRVTVuV1dKSmVESlJaMEZCZDBsa2RtcElXa05CUVVOdk9XNUxUV1JyU1VGQlNrSnNaRzk0TWxGblFVRmxUbEkxYWtoYVEwRkJRbWRSTXpKTlpHdEpRVUZGYVhsblNYZ3lVV2RCUVUxRFIwVnFTRnBEUVVGQldXdEpaVTFrYTBsQlFVRkVMMmx2ZURKUlowRkJOa2N5VDJwSVdrTkJRVVJSTTBwSFRXUnJTVUZCVEdoTWJGbDRNbEZuUVVGdlRIRlpha2hhUTBGQlEwbExXbmxOWkd0SlFVRklRMWx1TkhneVVXZEJRVmRCWldwcVNGcERRVUZDUVdSeFlVMWthMGxCUVVOcWJIRlplREpSWjBGQlJVWlRkR3BJV2tOQlFVUTBkM0pEVFdSclNVRkJUMEY0ZEVsNE1sRm5RVUY1UzBNemFraGFRMEZCUTNkRU4zVk5aR3RKUVVGS2FDdDJiM2d5VVdkQlFXZFBNMEpxU0ZwRFFVRkNiMWhOVjAxa2EwbEJRVVpFVEhsSmVESlJaMEZCVDBSeVRXcElXa05CUVVGbmNXTXJUV1JyU1VGQlFXZFpNRFI0TWxGblFVRTRTV0pYYWtoYVEwRkJSRms1WkcxTlpHdEpRVUZOUW1zeldYZ3lVV2RCUVhGT1VHZHFTRnBEUVVGRFVWRjFVMDFrYTBsQlFVaHBlRFUwZURKUlowRkJXVU5FY21wSVdrTkJRVUpKYWlzMlRXUnJTVUZCUkVRck9GbDRNbEZuUVVGSFJ6TXhha2hhUXlJc0ltUjBlWEJsSWpvaVpteHZZWFEyTkNJc0luTm9ZWEJsSWpwYk1UUTBYWDBzSW5raU9uc2lYMTl1WkdGeWNtRjVYMThpT2lKQlFVRkJaMHBSYjBRd1FVRkJRVVJCUWpoVlQxRkJRVUZCUTBJM1dWRTFRVUZCUVVGWlR6YzVSRlZCUVVGQlFXZFdOM2RPVVVGQlFVRk5ReTlsWnpGQlFVRkJRV2REWnpWRVZVRkJRVUZCUVRabU5FMVJRVUZCUVVsRGNIaEJlRUZCUVVGQlFVZHhTMFJGUVVGQlFVTkJOREJSVFZGQlFVRkJRMEprTDNkMFFVRkJRVUZ2VG1FMVF6QkJRVUZCUVVGWlUyTk5VVUZCUVVGSlJISnNRWGhCUVVGQlFUUklWVU5FVlVGQlFVRkRaME5FUVU5UlFVRkJRVVZEWWxoUk9VRkJRVUZCUVVwa1JrVkZRVUZCUVVKblNWVTBVVkZCUVVGQlMwTnlWbWhDUVVGQlFVRkJSRnBtUlVWQlFVRkJRV2RUVW5OUlVVRkJRVUZMUXpSeVp6bEJRVUZCUVRST05HMUVNRUZCUVVGRFp6bHdiMDlSUVVGQlFVbEJUMFIzTlVGQlFVRkJVVU5oUkVSVlFVRkJRVVJuWnpSRlRsRkJRVUZCUjBSb1puY3hRVUZCUVVGQlJEa3JSRlZCUVVGQlFVRkxXa1ZPVVVGQlFVRkRRVlJ3UVRGQlFVRkJRVWxRTWpKRVZVRkJRVUZFUVVwMFJVNVJRVUZCUVVkQ1VUWjNNVUZCUVVGQlFVaHZSa1JyUVVGQlFVSm5OMnhyVDFGQlFVRkJTMEpwY21jMVFVRkJRVUZCVG1ORFJEQkJRVUZCUW1kc04zTlFVVUZCUVVGUFFYSlBhRUpCUVVGQlFVRkplVmRGUlVGQlFVRkRRWGhpYTFGUlFVRkJRVTlFS3pOQ1FrRkJRVUZCV1VSblFVVlZRVUZCUVVGblJuVk5VVkZCUVVGQlQwUjZlRkpDUVVGQlFVRnZUa2R2UlVWQlFVRkJRVUZFWTBWUlVVRkJRVUZGUWtreVVrSkJRVUZCUVc5SlVIaEZSVUZCUVVGQlp6RXJWVkZSUVVGQlFVbEJjVEpvUWtGQlFVRkJRVWczVDBWRlFVRkJRVUpuUnpoalVWRkJRVUZCVDBNMGRuaENRVUZCUVVGUlJtRTBSVVZCUVVGQlEwRnljMjlSVVVGQlFVRkxRVWN6VWtKQlFVRkJRVFJHTjNaRlJVRkJRVUZDUVZkbmExSlJRVUZCUVVsQ1ZrbDRSa0ZCUVVGQk5FWkJPVVZWUVVGQlFVSm5jMGxOVWxGQlFVRkJUMEZRZVdoR1FVRkJRVUZaUnpoUlJXdEJRVUZCUkdjNWFsRlRVVUZCUVVGSlFpdFhVa3BCUVVGQlFVRkJXaXRGYTBGQlFVRkVaMWRGVFZOUlFVRkJRVXREY2tOQ1NrRkJRVUZCWjFBM1RrVlZRVUZCUVVKbmNEbEpVbEZCUVVGQlIwSlJNWGhHUVVGQlFVRlJVRzVpUlZWQlFVRkJRMmRWWW10U1VVRkJRVUZEUTNGc2FFWkJRVUZCUVdkQlNqQkZWVUZCUVVGQ1oxZHJkMUpSUVVGQlFVTkRlVXBDUmtGQlFVRkJRVUZ5T1VWRlFVRkJRVUZCVUU5elVWRkJRVUZCVDBKME1sSkNRVUZCUVVFMFNpOUlSVVZCUVVGQlJHY3JVSE5SVVVGQlFVRk5RbEpOUWtaQlFVRkJRWGRMY0d0RlZVRkJRVUZFUVc1MVRWSlJRVUZCUVV0RFUxbG9Ta0ZCUVVGQmIwbGlhRVZyUVVGQlFVRkJVM2RSVkZGQlFVRkJSMEZRU25oT1FVRkJRVUYzVGs1S1JUQkJRVUZCUTJkb2RqQlRVVUZCUVVGTFFUVnpVa3BCUVVGQlFXZFBlR3RGYTBGQlFVRkNaMjB3UlZOUlFVRkJRVWRDUzBob1NrRkJRVUZCVVZCdU5rVlZRVUZCUVVOblYzVkJVbEZCUVVGQlFVTTRlRkpHUVVGQlFVRlpRakp5UlZWQlFVRkJSRUY2TTJ0U1VVRkJRVUZCUTBOVFFrWkJRVUZCUVZsRVVWaEZWVUZCUVVGQlp6ZDFRVkZSUVVGQlFVRkRiM0ZvUWtGQlFVRkJkMGRHTUVWRlFVRkJRVUZuTlc0d1VWRkJRVUZCU1VKeGFIaENRVUZCUVVFMFR6WlJSVVZCUVVGQlJHZFhjbWRSVVVGQlFVRlBSRWN6ZUVKQlFVRkJRVFJFU1VoRlZVRkJRVUZCWjA5bVdWRlJRVUZCUVVWQkx6VlNRa0ZCUVVGQlowVllWVVZGUVVGQlFVRkJaVTFCVVZGQlFVRkJTME54Y2tKQ1FVRkJRVUZKVGpKWlJVVkJRVUZCUWtFMWIydFJVVUZCUVVGSFJIWmxhRUpCUVVGQlFXZFFhSEpGUlVGQlFVRkNaMjFHYTFGUlFVRkJRVVZCTkZKNFFrRkJRVUZCU1U1bk1FVkZRVUZCUVVOQlVuaHJVVkZCUVVGQlQwSjBLM2M1UVVGQlFVRjNSWHBGUkRCQlFVRkJRVUZUYVZWUlVVRkJRVUZMUW5SaFFrSkJRVUZCUVZGS1IzSkZSVUZCUVVGQlFYRjNSVkpSUVVGQlFVOUVSVlo0UmtGQlFVRkJiMDQyZEVWVlFVRkJRVUZuYWxGM1UxRkJRVUZCUzBFM1lYaEtRVUZCUVVGSlQzSktSV3RCUVVGQlEwRjBaV3RUVVVGQlFVRlBRMEZEVWs1QlFVRkJRVkZGZDNCRk1FRkJRVUZDUVZSRGExUlJRVUZCUVVWQ1RVdFNUa0VpTENKa2RIbHdaU0k2SW1ac2IyRjBOalFpTENKemFHRndaU0k2V3pFME5GMTlmU3dpYzJWc1pXTjBaV1FpT25zaWFXUWlPaUl4TVRReklpd2lkSGx3WlNJNklsTmxiR1ZqZEdsdmJpSjlMQ0p6Wld4bFkzUnBiMjVmY0c5c2FXTjVJanA3SW1sa0lqb2lNVEUwTWlJc0luUjVjR1VpT2lKVmJtbHZibEpsYm1SbGNtVnljeUo5ZlN3aWFXUWlPaUl4TURnMklpd2lkSGx3WlNJNklrTnZiSFZ0YmtSaGRHRlRiM1Z5WTJVaWZTeDdJbUYwZEhKcFluVjBaWE1pT250OUxDSnBaQ0k2SWpFd01qSWlMQ0owZVhCbElqb2lVR0Z1Vkc5dmJDSjlMSHNpWVhSMGNtbGlkWFJsY3lJNmV5SmtZWFJoWDNOdmRYSmpaU0k2ZXlKcFpDSTZJakV3TXpFaUxDSjBlWEJsSWpvaVEyOXNkVzF1UkdGMFlWTnZkWEpqWlNKOUxDSm5iSGx3YUNJNmV5SnBaQ0k2SWpFd016SWlMQ0owZVhCbElqb2lUR2x1WlNKOUxDSm9iM1psY2w5bmJIbHdhQ0k2Ym5Wc2JDd2liWFYwWldSZloyeDVjR2dpT201MWJHd3NJbTV2Ym5ObGJHVmpkR2x2Ymw5bmJIbHdhQ0k2ZXlKcFpDSTZJakV3TXpNaUxDSjBlWEJsSWpvaVRHbHVaU0o5TENKelpXeGxZM1JwYjI1ZloyeDVjR2dpT201MWJHd3NJblpwWlhjaU9uc2lhV1FpT2lJeE1ETTFJaXdpZEhsd1pTSTZJa05FVTFacFpYY2lmWDBzSW1sa0lqb2lNVEF6TkNJc0luUjVjR1VpT2lKSGJIbHdhRkpsYm1SbGNtVnlJbjBzZXlKaGRIUnlhV0oxZEdWeklqcDdJbk52ZFhKalpTSTZleUpwWkNJNklqRXhORGdpTENKMGVYQmxJam9pUTI5c2RXMXVSR0YwWVZOdmRYSmpaU0o5ZlN3aWFXUWlPaUl4TVRVeUlpd2lkSGx3WlNJNklrTkVVMVpwWlhjaWZTeDdJbUYwZEhKcFluVjBaWE1pT25zaVptOXliV0YwZEdWeUlqcDdJbWxrSWpvaU1UQXpOeUlzSW5SNWNHVWlPaUpDWVhOcFkxUnBZMnRHYjNKdFlYUjBaWElpZlN3aWNHeHZkQ0k2ZXlKcFpDSTZJakV3TURJaUxDSnpkV0owZVhCbElqb2lSbWxuZFhKbElpd2lkSGx3WlNJNklsQnNiM1FpZlN3aWRHbGphMlZ5SWpwN0ltbGtJam9pTVRBeE9DSXNJblI1Y0dVaU9pSkNZWE5wWTFScFkydGxjaUo5ZlN3aWFXUWlPaUl4TURFM0lpd2lkSGx3WlNJNklreHBibVZoY2tGNGFYTWlmU3g3SW1GMGRISnBZblYwWlhNaU9uc2liR2x1WlY5aGJIQm9ZU0k2TUM0eExDSnNhVzVsWDJOaGNDSTZJbkp2ZFc1a0lpd2liR2x1WlY5amIyeHZjaUk2SWlNeFpqYzNZalFpTENKc2FXNWxYMnB2YVc0aU9pSnliM1Z1WkNJc0lteHBibVZmZDJsa2RHZ2lPalVzSW5naU9uc2labWxsYkdRaU9pSjRJbjBzSW5raU9uc2labWxsYkdRaU9pSjVJbjE5TENKcFpDSTZJakV3TXpNaUxDSjBlWEJsSWpvaVRHbHVaU0o5WFN3aWNtOXZkRjlwWkhNaU9sc2lNVEF3TWlKZGZTd2lkR2wwYkdVaU9pSkNiMnRsYUNCQmNIQnNhV05oZEdsdmJpSXNJblpsY25OcGIyNGlPaUl4TGpBdU5DSjlmUW9nSUNBZ0lDQWdJRHd2YzJOeWFYQjBQZ29nSUNBZ0lDQWdJRHh6WTNKcGNIUWdkSGx3WlQwaWRHVjRkQzlxWVhaaGMyTnlhWEIwSWo0S0lDQWdJQ0FnSUNBZ0lDaG1kVzVqZEdsdmJpZ3BJSHNLSUNBZ0lDQWdJQ0FnSUNBZ2RtRnlJR1p1SUQwZ1puVnVZM1JwYjI0b0tTQjdDaUFnSUNBZ0lDQWdJQ0FnSUNBZ1FtOXJaV2d1YzJGbVpXeDVLR1oxYm1OMGFXOXVLQ2tnZXdvZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnS0daMWJtTjBhVzl1S0hKdmIzUXBJSHNLSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnWm5WdVkzUnBiMjRnWlcxaVpXUmZaRzlqZFcxbGJuUW9jbTl2ZENrZ2V3b2dJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQW9nSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0IyWVhJZ1pHOWpjMTlxYzI5dUlEMGdaRzlqZFcxbGJuUXVaMlYwUld4bGJXVnVkRUo1U1dRb0p6RTBNRFluS1M1MFpYaDBRMjl1ZEdWdWREc0tJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdkbUZ5SUhKbGJtUmxjbDlwZEdWdGN5QTlJRnQ3SW1SdlkybGtJam9pTlRBNE1UTTFPVEl0WXpZMU5pMDBOVE5rTFRoa056a3RNak5qTUdNMllUQXlOV1JrSWl3aWNtOXZkSE1pT25zaU1UQXdNaUk2SWpZNE0yRmlNak5rTFRrMU56RXROR0l5WWkwNE5UUmlMVFUyTlRFM09UWmpZekl6TWlKOWZWMDdDaUFnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJSEp2YjNRdVFtOXJaV2d1WlcxaVpXUXVaVzFpWldSZmFYUmxiWE1vWkc5amMxOXFjMjl1TENCeVpXNWtaWEpmYVhSbGJYTXBPd29nSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdDaUFnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJSDBLSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnYVdZZ0tISnZiM1F1UW05clpXZ2dJVDA5SUhWdVpHVm1hVzVsWkNrZ2V3b2dJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJR1Z0WW1Wa1gyUnZZM1Z0Wlc1MEtISnZiM1FwT3dvZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNCOUlHVnNjMlVnZXdvZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lIWmhjaUJoZEhSbGJYQjBjeUE5SURBN0NpQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdkbUZ5SUhScGJXVnlJRDBnYzJWMFNXNTBaWEoyWVd3b1puVnVZM1JwYjI0b2NtOXZkQ2tnZXdvZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdhV1lnS0hKdmIzUXVRbTlyWldnZ0lUMDlJSFZ1WkdWbWFXNWxaQ2tnZXdvZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0JsYldKbFpGOWtiMk4xYldWdWRDaHliMjkwS1RzS0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnWTJ4bFlYSkpiblJsY25aaGJDaDBhVzFsY2lrN0NpQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0I5Q2lBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQmhkSFJsYlhCMGN5c3JPd29nSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ2FXWWdLR0YwZEdWdGNIUnpJRDRnTVRBd0tTQjdDaUFnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lHTnZibk52YkdVdWJHOW5LQ0pDYjJ0bGFEb2dSVkpTVDFJNklGVnVZV0pzWlNCMGJ5QnlkVzRnUW05clpXaEtVeUJqYjJSbElHSmxZMkYxYzJVZ1FtOXJaV2hLVXlCc2FXSnlZWEo1SUdseklHMXBjM05wYm1jaUtUc0tJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ1kyeGxZWEpKYm5SbGNuWmhiQ2gwYVcxbGNpazdDaUFnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNCOUNpQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdmU3dnTVRBc0lISnZiM1FwQ2lBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUgwS0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUgwcEtIZHBibVJ2ZHlrN0NpQWdJQ0FnSUNBZ0lDQWdJQ0FnZlNrN0NpQWdJQ0FnSUNBZ0lDQWdJSDA3Q2lBZ0lDQWdJQ0FnSUNBZ0lHbG1JQ2hrYjJOMWJXVnVkQzV5WldGa2VWTjBZWFJsSUNFOUlDSnNiMkZrYVc1bklpa2dabTRvS1RzS0lDQWdJQ0FnSUNBZ0lDQWdaV3h6WlNCa2IyTjFiV1Z1ZEM1aFpHUkZkbVZ1ZEV4cGMzUmxibVZ5S0NKRVQwMURiMjUwWlc1MFRHOWhaR1ZrSWl3Z1ptNHBPd29nSUNBZ0lDQWdJQ0FnZlNrb0tUc0tJQ0FnSUNBZ0lDQThMM05qY21sd2RENEtJQ0FnSUFvZ0lEd3ZZbTlrZVQ0S0lDQUtQQzlvZEcxc1BnPT0iIHdpZHRoPSI3OTAiIHN0eWxlPSJib3JkZXI6bm9uZSAhaW1wb3J0YW50OyIgaGVpZ2h0PSIzMzAiPjwvaWZyYW1lPmApWzBdOwogICAgICAgICAgICAgICAgcG9wdXBfMDMwZWU5YjhhNGYzNDZiZGI0ZjA5ZGM5ZjU2NDc3YTcuc2V0Q29udGVudChpX2ZyYW1lXzYzYWRiMmU1NGVkNDQ4N2JiOWM2YWJkMmIyYmI1YjUzKTsKICAgICAgICAgICAgCgogICAgICAgICAgICBtYXJrZXJfNWJkMWU1M2U3ZTU0NDdiNjgzNThhM2VhNGFjYzUwNTUuYmluZFBvcHVwKHBvcHVwXzAzMGVlOWI4YTRmMzQ2YmRiNGYwOWRjOWY1NjQ3N2E3KQogICAgICAgICAgICA7CgogICAgICAgICAgICAKICAgICAgICAKICAgIAogICAgICAgICAgICB2YXIgbGF5ZXJfY29udHJvbF9mMGQyOGY3MTZkOWM0ODQ3OWFjY2IyMjk2M2VkMTA1YSA9IHsKICAgICAgICAgICAgICAgIGJhc2VfbGF5ZXJzIDogeyAib3BlbnN0cmVldG1hcCIgOiB0aWxlX2xheWVyX2RkM2I4ZTRkN2E4MDQ4OTE5ZDcxY2U3ZDRlNGMxNjA1LCB9LAogICAgICAgICAgICAgICAgb3ZlcmxheXMgOiB7ICJTZWEgU3VyZmFjZSBUZW1wZXJhdHVyZSIgOiBtYWNyb19lbGVtZW50XzY1YzcwNDc2MTcxZjQzMTY5Yjc2ZThlYmZkMjY1OWZiLCJDbHVzdGVyIiA6IG1hcmtlcl9jbHVzdGVyXzZkMzE4MmM5MzcwMjQ1MzI5YWUzYzA1YWI2YzlkMWY4LCB9CiAgICAgICAgICAgICAgICB9OwogICAgICAgICAgICBMLmNvbnRyb2wubGF5ZXJzKAogICAgICAgICAgICAgICAgbGF5ZXJfY29udHJvbF9mMGQyOGY3MTZkOWM0ODQ3OWFjY2IyMjk2M2VkMTA1YS5iYXNlX2xheWVycywKICAgICAgICAgICAgICAgIGxheWVyX2NvbnRyb2xfZjBkMjhmNzE2ZDljNDg0NzlhY2NiMjI5NjNlZDEwNWEub3ZlcmxheXMsCiAgICAgICAgICAgICAgICB7cG9zaXRpb246ICd0b3ByaWdodCcsCiAgICAgICAgICAgICAgICAgY29sbGFwc2VkOiB0cnVlLAogICAgICAgICAgICAgICAgIGF1dG9aSW5kZXg6IHRydWUKICAgICAgICAgICAgICAgIH0pLmFkZFRvKG1hcF84OGRiNTBiZmI5YmE0YzcxOTZhN2MxMzhlZjBjY2Y4OCk7CiAgICAgICAgICAgIAogICAgICAgIAo8L3NjcmlwdD4=" style="position:absolute;width:100%;height:100%;left:0;top:0;border:none !important;" allowfullscreen webkitallowfullscreen mozallowfullscreen></iframe></div></div>



Now we can navigate the map and click on the markers to explorer our findings.

The green markers locate the observations locations. They pop-up an interactive plot with the time-series and scores for the models (hover over the lines to se the scores). The blue markers indicate the nearest model grid point found for the comparison.
<br>
Right click and choose Save link as... to
[download](https://raw.githubusercontent.com/ioos/notebooks_demos/master/notebooks/2016-12-22-boston_light_swim.ipynb)
this notebook, or click [here](https://mybinder.org/v2/gh/ioos/notebooks_demos/master?filepath=notebooks/2016-12-22-boston_light_swim.ipynb) to run a live instance of this notebook.