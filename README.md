# `pysplot`: Spatiotemporal Plotting in Python

The `pysplot` package hosts tools and visualization methods for making spatial and temporal plots with heliophysics timeseries data.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

To install the package:
```
pip install pysplot
```
Install additional dependencies if running through the examples:
```
pip install pysplot[examples]
```
or using the WHPI-specific tools:
```
pip install pysplot[whpi]
```

## Usage

After installation, you can use the package as follows:

```python
from pysplot.io.data import SpatialData, ScienceData, SpatialTimeData
from pysplot.plots.plottypes import orbit_plot
import matplotlib.pyplot as plt

# input data as dictionaries
spatial_data_dictionary = {
    'x': ['2020-01-01 01:00', '2020-01-01 02:00', '2020-01-01 03:00'], 
    'y': [[0,1,2],[3,4,5],[6,7,8]]
}
science_data_dictionary = {
    'x': ['2020-01-01 02:03', '2020-01-01 02:30', '2020-01-01 02:45'], 
    'y': [[1.08, 14.14, -16.24], [1.39, 14.08, -17.80], [1.86, 13.56, -18.79]]
}

# Initiate data objects
spatial_data = SpatialData(spatial_data_dictionary, spatial_columns_names=['x','y','z'])
science_data = ScienceData(science_data_dictionary, science_columns_names=['Bx', 'By', 'Bz'])
combined_data = SpatialTimeData(spatial_data, science=science_data)

# Make plots
fig, ax = plt.subplots()
orbit_plot(combined_data.data, ax, 'x', 'y')
```

For a full example of how to use `pysplot`, refer to our [Getting Started notebook](examples/getting_started.ipynb).


## Contributing

We welcome contributions to improve this package! To contribute, please follow these steps:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Make your changes.
4. Commit your changes (`git commit -am 'Add new feature'`).
5. Push to the branch (`git push origin feature-branch`).
6. Create a pull request.

Please make sure your code adheres to the existing style guidelines and includes tests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
