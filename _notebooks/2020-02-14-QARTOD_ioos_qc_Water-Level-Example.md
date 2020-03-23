---
title: ""
layout: notebook

---
This post will demonstrate how to [run ``ioos_qc``](https://github.com/ioos/ioos_qc) on a time-series dataset. ``ioos_qc`` implements the [Quality Assurance / Quality Control of Real Time Oceanographic Data (QARTOD)](https://ioos.noaa.gov/project/qartod/).

We will [use `bokeh`](https://docs.bokeh.org/en/latest/) for interactive plots, so let's start by loading the interactive notebook output.

<div class="prompt input_prompt">
In&nbsp;[1]:
</div>

```python
from bokeh.plotting import output_notebook
output_notebook()
```



<div class="bk-root">
    <a href="https://bokeh.org" target="_blank" class="bk-logo bk-logo-small bk-logo-notebook"></a>
    <span id="1001">Loading BokehJS ...</span>
</div>




We will be using the water level data from a [fixed station in Kotzebue, AK](https://www.google.com/maps?q=66.895035,-162.566752).

Below we create a simple Quality Assurance/Quality Control (QA/QC) configuration that will be used as input for ``ioos_qc``. All the interval values are in the same units as the data.

For more information on the tests and recommended values for QA/QC check the documentation of each test and its inputs: 
https://ioos.github.io/ioos_qc/api/ioos_qc.html#module-ioos_qc.qartod

<div class="prompt input_prompt">
In&nbsp;[2]:
</div>

```python
variable_name = "sea_surface_height_above_sea_level_geoid_mhhw"


qc_config = {
    "qartod": {
      "gross_range_test": {
        "fail_span": [-10, 10],
        "suspect_span": [-2, 3]
      },
      "flat_line_test": {
        "tolerance": 0.001,
        "suspect_threshold": 10800,
        "fail_threshold": 21600
      },
      "spike_test": {
        "suspect_threshold": 0.8,
        "fail_threshold": 3,
      }
    }
}
```

Now we are ready to load the data, run tests and plot results!

We will get the data from the [AOOS ERDDAP server](http://erddap.aoos.org/erddap/). Note that the data may change in the future. For reproducibility's sake we will save the data downloaded into a CSV file.

<div class="prompt input_prompt">
In&nbsp;[3]:
</div>

```python
from pathlib import Path
import pandas as pd
from erddapy import ERDDAP


path = Path().absolute()
fname = path.joinpath("data", "water_level_example.csv")

if fname.is_file():
    data = pd.read_csv(fname, parse_dates=["time (UTC)"])
else:
    e = ERDDAP(
        server="http://erddap.aoos.org/erddap/",
        protocol="tabledap"
    )
    e.dataset_id = "kotzebue-alaska-water-level"
    e.constraints = {
        "time>=": "2018-09-05T21:00:00Z",
        "time<=": "2019-07-10T19:00:00Z",
    }
    e.variables = [
        variable_name,
        "time",
        "z",
    ]
    data = e.to_pandas(
        index_col="time (UTC)",
        parse_dates=True,
    )
    data["timestamp"] = data.index.astype("int64") // 1e9
    data.to_csv(fname)

data.head()
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
      <th>time (UTC)</th>
      <th>sea_surface_height_above_sea_level_geoid_mhhw (m)</th>
      <th>z (m)</th>
      <th>timestamp</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2018-09-05 21:00:00+00:00</td>
      <td>0.4785</td>
      <td>0.0</td>
      <td>1.536181e+09</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2018-09-05 22:00:00+00:00</td>
      <td>0.4420</td>
      <td>0.0</td>
      <td>1.536185e+09</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2018-09-05 23:00:00+00:00</td>
      <td>0.4968</td>
      <td>0.0</td>
      <td>1.536188e+09</td>
    </tr>
    <tr>
      <th>3</th>
      <td>2018-09-06 01:00:00+00:00</td>
      <td>0.5456</td>
      <td>0.0</td>
      <td>1.536196e+09</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2018-09-06 02:00:00+00:00</td>
      <td>0.5761</td>
      <td>0.0</td>
      <td>1.536199e+09</td>
    </tr>
  </tbody>
</table>
</div>



<div class="prompt input_prompt">
In&nbsp;[4]:
</div>

```python
from ioos_qc.config import QcConfig


qc = QcConfig(qc_config)

qc_results =  qc.run(
    inp=data["sea_surface_height_above_sea_level_geoid_mhhw (m)"],
    tinp=data["timestamp"],
    zinp=data["z (m)"],
    gen_agg=True
)

qc_results
```




    OrderedDict([('qartod',
                  OrderedDict([('gross_range_test',
                                masked_array(data=[1, 1, 1, ..., 1, 1, 1],
                                             mask=False,
                                       fill_value=999999,
                                            dtype=uint8)),
                               ('flat_line_test', array([1, 1, 1, ..., 1, 1, 1])),
                               ('spike_test',
                                masked_array(data=[1, 1, 1, ..., 1, 1, 1],
                                             mask=False,
                                       fill_value=999999,
                                            dtype=uint8))]))])



The results are returned in a dictionary format, similar to the input configuration, with a mask for each test. While the mask **is** a masked array it should not be applied as such. The results range from 1 to 4 meaning:

1. data passed the QA/QC
2. did not run on this data point
3. flag as suspect
4. flag as failed

Now we can write a plotting function that will read these results and flag the data.

<div class="prompt input_prompt">
In&nbsp;[5]:
</div>

```python
import numpy as np
from datetime import datetime

from bokeh.layouts import gridplot
from bokeh.plotting import figure, show


def plot_results(data, var_name, results, title, test_name):
    time = data["time (UTC)"]
    obs = data[var_name]
    qc_test = results["qartod"][test_name]

    qc_pass = np.ma.masked_where(qc_test != 1, obs)
    qc_suspect = np.ma.masked_where(qc_test != 3, obs)
    qc_fail = np.ma.masked_where(qc_test != 4, obs)
    qc_notrun = np.ma.masked_where(qc_test != 2, obs)

    p1 = figure(x_axis_type="datetime", title=test_name + " : " + title)
    p1.grid.grid_line_alpha=0.3
    p1.xaxis.axis_label = "Time"
    p1.yaxis.axis_label = "Observation Value"

    p1.line(time, obs,  legend_label="obs", color="#A6CEE3")
    p1.circle(time, qc_notrun, size=2, legend_label="qc not run", color="gray", alpha=0.2)
    p1.circle(time, qc_pass, size=4, legend_label="qc pass", color="green", alpha=0.5)
    p1.circle(time, qc_suspect, size=4, legend_label="qc suspect", color="orange", alpha=0.7)
    p1.circle(time, qc_fail, size=6, legend_label="qc fail", color="red", alpha=1.0)

    show(gridplot([[p1]], plot_width=800, plot_height=400))


title = "Water Level [MHHW] [m] : Kotzebue, AK"
```

The gross range test test should fail data outside the $\pm$ 10 range and suspect data below -2, and greater than 3. As one can easily see all the major spikes are flagged as expected.

<div class="prompt input_prompt">
In&nbsp;[6]:
</div>

```python
plot_results(
    data,
    "sea_surface_height_above_sea_level_geoid_mhhw (m)",
    qc_results,
    title,
    "gross_range_test"
)
```








<div class="bk-root" id="61664200-7302-48df-8440-c0eddaef53e6" data-root-id="1209"></div>





An actual spike test, based on a data increase threshold, flags similar spikes to the gross range test but also indetifies other suspect unsual increases in the series.

<div class="prompt input_prompt">
In&nbsp;[7]:
</div>

```python
plot_results(
    data,
    "sea_surface_height_above_sea_level_geoid_mhhw (m)",
    qc_results,
    title,
    "spike_test"
)
```








<div class="bk-root" id="8461cb95-e270-4a82-98b3-7fd023b036ea" data-root-id="1613"></div>





The flat line test identifies issues with the data where values are "stuck."

`ioos_qc` succefully identified a huge portion of the data where that happens and flagged a smaller one as suspect. (Zoom in the red point to the left to see this one.)

<div class="prompt input_prompt">
In&nbsp;[8]:
</div>

```python
plot_results(
    data,
    "sea_surface_height_above_sea_level_geoid_mhhw (m)",
    qc_results,
    title,
    "flat_line_test"
)
```








<div class="bk-root" id="0def68ea-de5a-4140-8207-8bfe2b6bff84" data-root-id="2045"></div>





This notebook was adapt from Jessica Austin and Kyle Wilcox's [original ioos_qc examples](https://github.com/ioos/ioos_qc/blob/b34b3762d659362fb3af11f52d8905d18cd6ec7b/docs/source/examples/QartodTestExample_WaterLevel.ipynb). Please [see the ``ioos_qc`` documentation](https://ioos.github.io/ioos_qc/) for more examples.
<br>
Right click and choose Save link as... to
[download](https://raw.githubusercontent.com/ioos/notebooks_demos/master/notebooks/2020-02-14-QARTOD_ioos_qc_Water-Level-Example.ipynb)
this notebook, or click [here](https://binder.pangeo.io/v2/gh/ioos/notebooks_demos/master?filepath=notebooks/2020-02-14-QARTOD_ioos_qc_Water-Level-Example.ipynb) to run a live instance of this notebook.